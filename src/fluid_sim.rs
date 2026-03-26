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
//!   4. `enforce_walls`      – hard-set no-penetration on top/bottom walls and solid cells.

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
    /// Solid mask: true = solid obstacle (no flow inside).  Size: width × height.
    pub solid: Matrix<bool>,
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
            solid: Matrix::new(w, h, false),
        }
    }
    /// Calcule la force exercée par le fluide sur le solide.
    /// Retourne (force_horizontale, force_verticale) = (traînée, portance)
    /// (forces par unité de profondeur, car 2D)
    pub fn compute_forces(&self) -> (f32, f32) {
        let mut fx = 0.0;
        let mut fy = 0.0;
        let (w, h) = self.cell_count;
        let dx = self.cell_size;
        let mu = 0.0; // viscosité dynamique (non utilisée actuellement)
        // Si vous voulez inclure la viscosité, il faudrait définir mu = densité * nu
        // Mais votre simulateur n'a pas de viscosité explicite, donc on met 0.

        // On parcourt toutes les cellules solides
        for x in 0..w {
            for y in 0..h {
                if !self.solid[(x, y)] {
                    continue;
                }

                // Pour chaque face de la cellule, on regarde si le voisin est du fluide
                // Si oui, on ajoute la contribution de pression et de cisaillement.

                // Face gauche : x, y (normale sortante du solide : (-1, 0))
                if !self.is_solid(x.saturating_sub(1), y) {
                    // Pression sur la face (on prend la pression de la cellule voisine)
                    let p = self.pressure_at(x as i32 - 1, y as i32);
                    fx += p * dx; // contribution horizontale : p * n_x * surface (dx * 1)
                }

                // Face droite : x+1, y (normale (1, 0))
                if !self.is_solid(x + 1, y) {
                    let p = self.pressure_at(x as i32 + 1, y as i32);
                    fx -= p * dx; // car n_x = +1, la force est -p*dx (car p agit dans la direction opposée)
                }

                // Face basse : x, y (normale (0, -1))
                if !self.is_solid(x, y.saturating_sub(1)) {
                    let p = self.pressure_at(x as i32, y as i32 - 1);
                    fy += p * dx; // n_y = -1, donc force = p * dx (vers le bas)
                }

                // Face haute : x, y+1 (normale (0, 1))
                if !self.is_solid(x, y + 1) {
                    let p = self.pressure_at(x as i32, y as i32 + 1);
                    fy -= p * dx; // n_y = +1, force = -p*dx (vers le haut)
                }
            }
        }

        // portance positive vers le haut
        (fx, fy)
    }

    // ── Helpers ──────────────────────────────────────────────────────────────

    /// Returns true if the cell is outside the domain OR is a solid obstacle.
    fn is_solid(&self, x: u32, y: u32) -> bool {
        x >= self.cell_count.0 || y >= self.cell_count.1 || self.solid[(x, y)]
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

    /// Pressure at (x,y), returning 0 for out-of-domain or solid cells.
    fn pressure_at(&self, x: i32, y: i32) -> f32 {
        if x < 0 || y < 0 || x >= self.cell_count.0 as i32 || y >= self.cell_count.1 as i32 {
            return 0.0;
        }
        if self.solid[(x as u32, y as u32)] {
            return 0.0;
        }
        self.pressure_map[(x as u32, y as u32)]
    }

    // ── Velocity sampling ────────────────────────────────────────────────────

    pub fn sample_vx(&self, px: f32, py: f32) -> f32 {
        let gx = px.clamp(0.0, self.cell_count.0 as f32);
        let gy = (py - 0.5).clamp(0.0, (self.cell_count.1 - 1) as f32);
        bilinear(&self.velocities_x, gx, gy)
    }

    pub fn sample_vy(&self, px: f32, py: f32) -> f32 {
        let gx = (px - 0.5).clamp(0.0, (self.cell_count.0 - 1) as f32);
        let gy = py.clamp(0.0, self.cell_count.1 as f32);
        bilinear(&self.velocities_y, gx, gy)
    }

    pub fn velocity_at_center(&self, x: u32, y: u32) -> (f32, f32) {
        let px = x as f32 + 0.5;
        let py = y as f32 + 0.5;
        (self.sample_vx(px, py), self.sample_vy(px, py))
    }

    // ── Advection ────────────────────────────────────────────────────────────

    pub fn advect_velocities(&mut self) {
        let vx_old = self.velocities_x.clone();
        let vy_old = self.velocities_y.clone();
        let (w, h) = self.cell_count;
        let dt = self.time_step;

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
        for x in 1..=(w) {
            for y in 0..h {
                // Skip faces adjacent to solid cells (will be zeroed by enforce_solid).
                let left_solid = x > 0 && x <= w && self.solid[(x - 1, y)];
                let right_solid = x < w && self.solid[(x, y)];
                if left_solid || right_solid {
                    continue;
                }
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
        for x in 0..w {
            for y in 1..h {
                let below_solid = self.solid[(x, y - 1)];
                let above_solid = self.solid[(x, y)];
                if below_solid || above_solid {
                    continue;
                }
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

        let p_sum = self.pressure_at(xi, yi + 1) * flow.0   // top
            + self.pressure_at(xi - 1, yi) * flow.1          // left
            + self.pressure_at(xi + 1, yi) * flow.2          // right
            + self.pressure_at(xi, yi - 1) * flow.3; // bottom

        let div = self.velocities_x[(x + 1, y)] * flow.2 - self.velocities_x[(x, y)] * flow.1
            + self.velocities_y[(x, y + 1)] * flow.0
            - self.velocities_y[(x, y)] * flow.3;

        self.pressure_map[(x, y)] =
            (p_sum - self.density * self.cell_size * div / self.time_step) / fluid_count;
    }

    fn update_velocities(&mut self) {
        let k = self.time_step / (self.cell_size * self.density);
        let (w, h) = self.cell_count;

        for x in 1..=(w) {
            for y in 0..h {
                let p_right = self.pressure_at(x as i32, y as i32);
                let p_left = self.pressure_at(x as i32 - 1, y as i32);
                self.velocities_x[(x, y)] -= k * (p_right - p_left);
            }
        }

        for x in 0..w {
            for y in 1..h {
                let p_top = self.pressure_at(x as i32, y as i32);
                let p_bottom = self.pressure_at(x as i32, y as i32 - 1);
                self.velocities_y[(x, y)] -= k * (p_top - p_bottom);
            }
        }
    }

    fn anchor_pressure(&mut self) {
        let n = (self.cell_count.0 * self.cell_count.1) as f32;
        let mean = self.pressure_map.sum() / n;
        for (x, y) in self.pressure_map.indices() {
            self.pressure_map[(x, y)] -= mean;
        }
    }

    // ── Boundary conditions ───────────────────────────────────────────────────

    /// No-penetration on domain walls AND solid obstacle surfaces.
    pub fn enforce_walls(&mut self) {
        let (w, h) = self.cell_count;

        // Top / bottom domain walls (y-velocity only).
        for x in 0..w {
            self.velocities_y[(x, 0)] = 0.0;
            self.velocities_y[(x, h)] = 0.0;
        }

        // Solid obstacle: zero all velocity components that touch a solid cell.
        for x in 0..w {
            for y in 0..h {
                if !self.solid[(x, y)] {
                    continue;
                }
                // Zero the four surrounding face velocities.
                self.velocities_x[(x, y)] = 0.0;
                self.velocities_x[(x + 1, y)] = 0.0;
                self.velocities_y[(x, y)] = 0.0;
                self.velocities_y[(x, y + 1)] = 0.0;
                self.pressure_map[(x, y)] = 0.0;
            }
        }
    }

    // ── Public simulation step ────────────────────────────────────────────────

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
                if !self.solid[(x, y)] {
                    max = max.max(self.divergence_at(x, y).abs());
                }
            }
        }
        max
    }
}
