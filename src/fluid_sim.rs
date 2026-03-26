//! MAC (Marker-And-Cell) staggered-grid incompressible fluid simulator.
//!
//! Grid layout for a cell (x, y) that occupies [x, x+1] × [y, y+1]:
//!   • pressure  p    stored at cell centre  (x+0.5, y+0.5)
//!   • x-velocity u   stored at left/right faces  (x,     y+0.5)  → `velocities_x`
//!   • y-velocity v   stored at bottom/top  faces  (x+0.5, y    )  → `velocities_y`
//!
//! Simulation step (per frame):
//!   1. Caller sets inflow / boundary values.
//!   2. `advect_velocities`  – semi-Lagrangian advection of u and v.
//!   3. `pressure_projection`– Gauss-Seidel pressure solve + velocity correction.
//!   4. `enforce_walls`      – hard-set no-penetration on top/bottom walls.

use crate::Direction::{self, Down, Left, Right, Up};
use crate::Matrix;
use crate::matrix::bilinear;

pub struct FluidGrid {
    pub time_step: f32,
    pub density: f32,
    /// (width, height) in cells.
    pub cell_count: (u32, u32),
    pub cell_size: f32,

    /// Horizontal velocity at x-faces.  Size: (width+1) × height.
    pub velocities_x: Matrix<f32>,
    /// Vertical velocity at y-faces.    Size: width × (height+1).
    pub velocities_y: Matrix<f32>,
    /// Pressure at cell centres.        Size: width × height.
    pub pressure_map: Matrix<f32>,
}

impl FluidGrid {
    pub fn new(cell_count: (u32, u32), cell_size: f32) -> FluidGrid {
        let (w, h) = cell_count;
        FluidGrid {
            time_step: 0.1,
            density: 1.0,
            cell_count,
            cell_size,
            velocities_x: Matrix::new(w + 1, h, 0.5),
            velocities_y: Matrix::new(w, h + 1, 0.5),
            pressure_map: Matrix::new(w, h, 0.0),
        }
    }

    // ── Helpers ──────────────────────────────────────────────────────────────

    fn is_solid(&self, x: u32, y: u32) -> bool {
        x >= self.cell_count.0 || y >= self.cell_count.1
    }

    fn is_fluid_edge(&self, x: u32, y: u32, dir: Direction) -> bool {
        match (dir, x, y) {
            (Left, 0, _) => false,
            (Down, _, 0) => false,
            (Left, _, _) => !self.is_solid(x - 1, y),
            (Down, _, _) => !self.is_solid(x, y - 1),
            (Up, _, _) => !self.is_solid(x, y + 1),
            (Right, _, _) => !self.is_solid(x + 1, y),
        }
    }

    fn fluid_edge_f(&self, x: u32, y: u32, dir: Direction) -> f32 {
        if self.is_fluid_edge(x, y, dir) {
            1.0
        } else {
            0.0
        }
    }

    fn flow_around_cell(&self, x: u32, y: u32) -> (f32, f32, f32, f32) {
        (
            self.fluid_edge_f(x, y, Up),
            self.fluid_edge_f(x, y, Left),
            self.fluid_edge_f(x, y, Right),
            self.fluid_edge_f(x, y, Down),
        )
    }

    /// Pressure at (x,y), returning 0 for out-of-domain cells (ghost pressure = 0).
    fn pressure_at(&self, x: i32, y: i32) -> f32 {
        if x < 0 || y < 0 || x >= self.cell_count.0 as i32 || y >= self.cell_count.1 as i32 {
            0.0
        } else {
            self.pressure_map[(x as u32, y as u32)]
        }
    }

    // ── Velocity sampling ────────────────────────────────────────────────────

    /// Sample u (x-velocity) at arbitrary grid position (px, py).
    ///
    /// u is stored at face (i, j+0.5), so in storage coordinates the sample
    /// point becomes (px,  py − 0.5).
    pub fn sample_vx(&self, px: f32, py: f32) -> f32 {
        let gx = px.clamp(0.0, self.cell_count.0 as f32);
        let gy = (py - 0.5).clamp(0.0, (self.cell_count.1 - 1) as f32);
        bilinear(&self.velocities_x, gx, gy)
    }

    /// Sample v (y-velocity) at arbitrary grid position (px, py).
    ///
    /// v is stored at face (i+0.5, j), so storage coordinates become (px − 0.5, py).
    pub fn sample_vy(&self, px: f32, py: f32) -> f32 {
        let gx = (px - 0.5).clamp(0.0, (self.cell_count.0 - 1) as f32);
        let gy = py.clamp(0.0, self.cell_count.1 as f32);
        bilinear(&self.velocities_y, gx, gy)
    }

    /// Velocity vector at cell-centre position.
    pub fn velocity_at_center(&self, x: u32, y: u32) -> (f32, f32) {
        let px = x as f32 + 0.5;
        let py = y as f32 + 0.5;
        (self.sample_vx(px, py), self.sample_vy(px, py))
    }

    // ── Advection ────────────────────────────────────────────────────────────

    /// Semi-Lagrangian advection of both velocity components.
    ///
    /// Fixed bugs vs. original:
    ///   1. Uses a clone so earlier cells don't corrupt later ones.
    ///   2. Face positions are correct: u-face at (x, y+0.5), v-face at (x+0.5, y).
    ///   3. Bilinear interpolation in both axes.
    pub fn advect_velocities(&mut self) {
        let vx_old = self.velocities_x.clone();
        let vy_old = self.velocities_y.clone();
        let (w, h) = self.cell_count;
        let dt = self.time_step;

        // Closures sampling from the *old* (pre-step) velocity fields.
        let old_vx = |px: f32, py: f32| -> f32 {
            let gx = px.clamp(0.0, w as f32);
            let gy = (py - 0.5).clamp(0.0, (h - 1) as f32);
            bilinear(&vx_old, gx, gy)
        };
        let old_vy = |px: f32, py: f32| -> f32 {
            let gx = (px - 0.5).clamp(0.0, (w - 1) as f32);
            let gy = py.clamp(0.0, h as f32);
            bilinear(&vy_old, gx, gy)
        };

        // ---- Advect u (x-faces) -------------------------------------------------
        // u[(x,y)] lives at continuous position (x, y+0.5).
        // Skip x=0 (inflow, set by caller each step).
        for x in 1..=(w) {
            for y in 0..h {
                let px = x as f32;
                let py = y as f32 + 0.5;
                let vel_x = old_vx(px, py);
                let vel_y = old_vy(px, py);
                let prev_x = (px - dt * vel_x).clamp(0.0, w as f32);
                let prev_y = (py - dt * vel_y).clamp(0.5, h as f32 - 0.5);
                self.velocities_x[(x, y)] = old_vx(prev_x, prev_y);
            }
        }

        // ---- Advect v (y-faces) -------------------------------------------------
        // v[(x,y)] lives at continuous position (x+0.5, y).
        // Skip y=0 and y=h (top/bottom walls, enforced later).
        for x in 0..w {
            for y in 1..h {
                let px = x as f32 + 0.5;
                let py = y as f32;
                let vel_x = old_vx(px, py);
                let vel_y = old_vy(px, py);
                let prev_x = (px - dt * vel_x).clamp(0.5, w as f32 - 0.5);
                let prev_y = (py - dt * vel_y).clamp(0.0, h as f32);
                self.velocities_y[(x, y)] = old_vy(prev_x, prev_y);
            }
        }
    }

    // ── Pressure projection ───────────────────────────────────────────────────

    /// One Gauss-Seidel sweep over all cells.
    fn pressure_sweep(&mut self) {
        for x in 0..self.cell_count.0 {
            for y in 0..self.cell_count.1 {
                self.pressure_solve_cell(x, y);
            }
        }
    }

    fn pressure_solve_cell(&mut self, x: u32, y: u32) {
        if self.is_solid(x, y) {
            return;
        }
        let flow = self.flow_around_cell(x, y);
        let fluid_count = flow.0 + flow.1 + flow.2 + flow.3;
        if fluid_count == 0.0 {
            return;
        }

        let xi = x as i32;
        let yi = y as i32;

        // Neighbour pressures (ghost = 0 outside domain).
        let p_sum = self.pressure_at(xi, yi + 1) * flow.0   // top
            + self.pressure_at(xi - 1, yi) * flow.1          // left
            + self.pressure_at(xi + 1, yi) * flow.2          // right
            + self.pressure_at(xi, yi - 1) * flow.3; // bottom

        // Divergence (scaled by cell_size to get velocity sum).
        let div = (self.velocities_x[(x + 1, y)] - self.velocities_x[(x, y)]) * flow.2
            + (self.velocities_x[(x + 1, y)] - self.velocities_x[(x, y)]) * 0.0  // handled above
            + self.velocities_x[(x + 1, y)] * flow.2
            - self.velocities_x[(x, y)] * flow.1
            + self.velocities_y[(x, y + 1)] * flow.0
            - self.velocities_y[(x, y)] * flow.3;

        // Standard Gauss-Seidel pressure update derived from the discrete
        // incompressibility constraint: ∇·u = 0.
        //   p[i,j] = ( Σ p_neighbour − ρ·h·div / dt ) / num_fluid_neighbours
        self.pressure_map[(x, y)] =
            (p_sum - self.density * self.cell_size * div / self.time_step) / fluid_count;
    }

    /// Subtract pressure gradient from velocity field (makes flow divergence-free).
    fn update_velocities(&mut self) {
        let k = self.time_step / (self.cell_size * self.density);
        let (w, h) = self.cell_count;

        // u-faces: skip x=0 (inflow set by caller).
        for x in 1..=(w) {
            for y in 0..h {
                let p_right = self.pressure_at(x as i32, y as i32);
                let p_left = self.pressure_at(x as i32 - 1, y as i32);
                self.velocities_x[(x, y)] -= k * (p_right - p_left);
            }
        }

        // v-faces: skip y=0 and y=h (walls enforced separately).
        for x in 0..w {
            for y in 1..h {
                let p_top = self.pressure_at(x as i32, y as i32);
                let p_bottom = self.pressure_at(x as i32, y as i32 - 1);
                self.velocities_y[(x, y)] -= k * (p_top - p_bottom);
            }
        }
    }

    /// Subtract the mean pressure so it doesn't drift.
    fn anchor_pressure(&mut self) {
        let n = (self.cell_count.0 * self.cell_count.1) as f32;
        let mean = self.pressure_map.sum() / n;
        for (x, y) in self.pressure_map.indices() {
            self.pressure_map[(x, y)] -= mean;
        }
    }

    // ── Boundary conditions ───────────────────────────────────────────────────

    /// No-penetration on top and bottom walls; open outflow on the right.
    /// The left wall (inflow) is set by the caller in main every frame.
    pub fn enforce_walls(&mut self) {
        let (w, h) = self.cell_count;
        for x in 0..w {
            self.velocities_y[(x, 0)] = 0.0;
            self.velocities_y[(x, h)] = 0.0;
        }
    }

    // ── Public simulation step ────────────────────────────────────────────────

    /// Full simulation step.  Call `inject_inflow` before and after.
    pub fn step(&mut self, pressure_iters: usize) {
        self.advect_velocities();
        for _ in 0..pressure_iters {
            self.pressure_sweep();
        }
        self.anchor_pressure();
        self.update_velocities();
        self.enforce_walls();
    }

    // ── Diagnostics ──────────────────────────────────────────────────────────

    pub fn divergence_at(&self, x: u32, y: u32) -> f32 {
        (self.velocities_x[(x + 1, y)] - self.velocities_x[(x, y)] + self.velocities_y[(x, y + 1)]
            - self.velocities_y[(x, y)])
            / self.cell_size
    }

    pub fn max_divergence(&self) -> f32 {
        let (w, h) = self.cell_count;
        let mut max = 0.0f32;
        for x in 0..w {
            for y in 0..h {
                max = max.max(self.divergence_at(x, y).abs());
            }
        }
        max
    }
}
