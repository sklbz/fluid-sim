use crate::Direction;
use crate::Matrix;
use Direction::*;

pub struct FluidGrid {
    pub time_step: f32,
    pub density: f32,
    pub cell_count: (u32, u32), /*(x,y)*/
    pub cell_size: f32,
    pub velocities_x: Matrix<f32>,
    pub velocities_y: Matrix<f32>,
    pub pressure_map: Matrix<f32>,
}

impl FluidGrid {
    pub fn new(cell_count: (u32, u32), cell_size: f32) -> FluidGrid {
        FluidGrid {
            time_step: 0.1,
            density: 1.0,
            cell_count,
            cell_size,
            velocities_x: Matrix::new(cell_count.0 + 1, cell_count.1, 0.0),
            velocities_y: Matrix::new(cell_count.0, cell_count.1 + 1, 0.0),
            pressure_map: Matrix::new(cell_count.0, cell_count.1, 0.0),
        }
    }

    pub fn divergence_at_cell(&self, x: u32, y: u32) -> f32 {
        let flow = self.flow_around_cell(x, y);
        let (top_velocity, left_velocity, right_velocity, bottom_velocity) =
            self.get_neighbor_velocity(x, y, flow);

        // rate of change of fluid velocity in either axis
        let gradient_x = (right_velocity - left_velocity) / self.cell_size;
        let gradient_y = (top_velocity - bottom_velocity) / self.cell_size;

        gradient_x + gradient_y
    }

    pub fn get_pressure(&self, x: u32, y: u32) -> f32 {
        if x >= self.cell_count.0 || y >= self.cell_count.1 {
            0.0
        } else {
            self.pressure_map[(x, y)]
        }
    }

    fn is_solid(&self, x: u32, y: u32) -> bool {
        /* x == 0 || y == 0 ||*/
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

    fn is_fluid_edge_float(&self, x: u32, y: u32, dir: Direction) -> f32 {
        if self.is_fluid_edge(x, y, dir) {
            1.0
        } else {
            0.0
        }
    }

    fn get_neighbor_pressure(
        &self,
        x: u32,
        y: u32,
        flow: (f32, f32, f32, f32),
    ) -> (f32, f32, f32, f32) {
        let top_pressure = self.get_pressure(x, y + 1);
        let left_pressure = if x == 0 {
            0.0
        } else {
            self.get_pressure(x - 1, y)
        };
        let right_pressure = self.get_pressure(x + 1, y);
        let bottom_pressure = if y == 0 {
            0.0
        } else {
            self.get_pressure(x, y - 1)
        };

        (
            top_pressure * flow.0,
            left_pressure * flow.1,
            right_pressure * flow.2,
            bottom_pressure * flow.3,
        )
    }

    fn get_neighbor_velocity(
        &self,
        x: u32,
        y: u32,
        flow: (f32, f32, f32, f32),
    ) -> (f32, f32, f32, f32) {
        let top_velocity = self.velocities_y[(x, y + 1)] * flow.0;
        let left_velocity = self.velocities_x[(x, y)] * flow.1;
        let right_velocity = self.velocities_x[(x + 1, y)] * flow.2;
        let bottom_velocity = self.velocities_y[(x, y)] * flow.3;

        (top_velocity, left_velocity, right_velocity, bottom_velocity)
    }

    fn flow_around_cell(&self, x: u32, y: u32) -> (f32, f32, f32, f32) {
        let flow_top = self.is_fluid_edge_float(x, y, Up);
        let flow_left = self.is_fluid_edge_float(x, y, Left);
        let flow_right = self.is_fluid_edge_float(x, y, Right);
        let flow_bottom = self.is_fluid_edge_float(x, y, Down);

        (flow_top, flow_left, flow_right, flow_bottom)
    }

    pub fn pressure_solve_cell(&mut self, x: u32, y: u32) {
        let flow = self.flow_around_cell(x, y);

        let fluid_edge_count = flow.0 + flow.1 + flow.2 + flow.3;
        if self.is_solid(x, y) || fluid_edge_count == 0.0 {
            self.pressure_map[(x, y)] = 0.0;
            return;
        }

        let (top_pressure, left_pressure, right_pressure, bottom_pressure) =
            self.get_neighbor_pressure(x, y, flow);
        let (top_velocity, left_velocity, right_velocity, bottom_velocity) =
            self.get_neighbor_velocity(x, y, flow);

        let pressure_sum = top_pressure + left_pressure + right_pressure + bottom_pressure;
        let delta_velocity_sum = right_velocity - left_velocity + top_velocity - bottom_velocity;

        self.pressure_map[(x, y)] = (pressure_sum
            - self.density * self.cell_size * delta_velocity_sum / self.time_step)
            / fluid_edge_count;
    }

    pub fn pressure_solve(&mut self) {
        for x in 0..self.cell_count.0 {
            for y in 0..self.cell_count.1 {
                self.pressure_solve_cell(x, y);
            }
        }

        self.anchor_pressure();
    }

    fn anchor_pressure(&mut self) {
        // Ancrer la pression : fixer P[0,0] = 0 comme référence
        // et soustraire la moyenne pour stabiliser
        let mean = self.pressure_map.sum() / (self.cell_count.0 * self.cell_count.1) as f32;
        for (x, y) in self.pressure_map.indices() {
            self.pressure_map[(x, y)] -= mean;
        }
    }

    pub fn update_velocities(&mut self) {
        let k: f32 = self.time_step / (self.cell_size * self.density);

        // ---- Horizontal -----------
        for (x, y) in self.velocities_x.indices() {
            let edge_is_solid = self.is_solid(x, y)
                || x.checked_sub(1)
                    .map(|px| self.is_solid(px, y))
                    .unwrap_or(true);
            if edge_is_solid {
                self.velocities_x[(x, y)] = 0.0;
                continue;
            }

            let pressure_right = self.get_pressure(x, y);
            let pressure_left = x
                .checked_sub(1)
                .map(|px| self.get_pressure(px, y))
                .unwrap_or(0.0);
            self.velocities_x[(x, y)] -= k * (pressure_right - pressure_left);
        }

        // ---- Vertical -------------
        for (x, y) in self.velocities_y.indices() {
            let edge_is_solid = self.is_solid(x, y)
                || y.checked_sub(1)
                    .map(|py| self.is_solid(x, py))
                    .unwrap_or(true);
            if edge_is_solid {
                self.velocities_y[(x, y)] = 0.0;
                continue;
            }

            let pressure_top = self.get_pressure(x, y);
            let pressure_bottom = y
                .checked_sub(1)
                .map(|py| self.get_pressure(x, py))
                .unwrap_or(0.0);
            self.velocities_y[(x, y)] -= k * (pressure_top - pressure_bottom);
        }
    }

    pub fn gauss_seidel(&mut self) {
        for _ in 0..10 {
            self.pressure_solve();
        }

        self.update_velocities();
    }
}
