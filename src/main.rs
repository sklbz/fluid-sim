pub mod matrix;

use matrix::Matrix;

fn main() {
    let mut grid = FluidGrid::new((10, 10), 1.0);
    grid.velocities_x.randomize(-1.0..=1.0);
    grid.velocities_y.randomize(-1.0..=1.0);

    println!("Divergence:");
    for x in 0..grid.cell_count.0 {
        for y in 0..grid.cell_count.1 {
            print!("{:^5.1}", grid.divergence_at_cell(x, y));
        }
        println!();
    }
    println!("Pressure:");
    for x in 0..grid.cell_count.0 {
        for y in 0..grid.cell_count.1 {
            print!("{:^5.1}", grid.get_pressure(x, y));
        }
        println!();
    }

    grid.pressure_solve();

    println!("Divergence:");
    for x in 0..grid.cell_count.0 {
        for y in 0..grid.cell_count.1 {
            print!("{:^5.1}", grid.divergence_at_cell(x, y));
        }
        println!();
    }
    println!("Pressure:");
    for x in 0..grid.cell_count.0 {
        for y in 0..grid.cell_count.1 {
            print!("{:^5.1}", grid.get_pressure(x, y));
        }
        println!();
    }

    grid.pressure_solve();

    println!("Divergence:");
    for x in 0..grid.cell_count.0 {
        for y in 0..grid.cell_count.1 {
            let div = grid.divergence_at_cell(x, y);
            if div > 0.0 {
                print!("\x1b[31mRouge\x1b[0m");
            } else {
                print!("\x1b[31mBleu\x1b[0m");
            }
            print!("{:^5.1}", div);
        }
            print!("\x1b[0m");
        }
        println!();
    }
    println!("Pressure:");
    for x in 0..grid.cell_count.0 {
        for y in 0..grid.cell_count.1 {
            print!("{:^5.1}", grid.get_pressure(x, y));
        }
        println!();
    }

    grid.pressure_solve();

    println!("Divergence:");
    for x in 0..grid.cell_count.0 {
        for y in 0..grid.cell_count.1 {
            print!("{:^5.1}", grid.divergence_at_cell(x, y));
        }
        println!();
    }
    println!("Pressure:");
    for x in 0..grid.cell_count.0 {
        for y in 0..grid.cell_count.1 {
            print!("{:^5.1}", grid.get_pressure(x, y));
        }
        println!();
    }
}

struct FluidGrid {
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
        let top_velocity = self.velocities_y[(x, y + 1)];
        let left_velocity = self.velocities_x[(x, y)];
        let right_velocity = self.velocities_x[(x + 1, y)];
        let bottom_velocity = self.velocities_y[(x, y)];

        // rate of change of fluid velocity in either axis
        let gradient_x = (right_velocity - left_velocity) / self.cell_size;
        let gradient_y = (top_velocity - bottom_velocity) / self.cell_size;

        gradient_x + gradient_y
    }

    fn get_pressure(&self, x: u32, y: u32) -> f32 {
        if x < 0 || x >= self.cell_count.0 || y < 0 || y >= self.cell_count.1 {
            0.0
        } else {
            self.pressure_map[(x, y)]
        }
    }

    pub fn pressure_solve_cell(&mut self, x: u32, y: u32) {
        let top_pressure = self.get_pressure(x, y + 1);
        let left_pressure = self.get_pressure(x, y);
        let right_pressure = self.get_pressure(x + 1, y);
        let bottom_pressure = self.get_pressure(x, y);

        let top_velocity = self.velocities_y[(x, y + 1)];
        let left_velocity = self.velocities_x[(x, y)];
        let right_velocity = self.velocities_x[(x + 1, y)];
        let bottom_velocity = self.velocities_y[(x, y)];

        let pressure_sum = top_pressure + left_pressure + right_pressure + bottom_pressure;
        let delta_velocity_sum = right_velocity - left_velocity + top_velocity - bottom_velocity;

        self.pressure_map[(x, y)] = (pressure_sum
            - self.density * self.cell_size * delta_velocity_sum / self.time_step)
            / 4.0;
    }

    pub fn pressure_solve(&mut self) {
        for x in 0..self.cell_count.0 {
            for y in 0..self.cell_count.1 {
                self.pressure_solve_cell(x, y);
            }
        }
    }

    pub fn update_velocities(&mut self) {
        let k: f32 = self.time_step / (self.cell_size * self.density);

        // ---- Horizontal -----------
        for x in 0..self.velocities_x.width() {
            for y in 0..self.velocities_x.height() {
                let pressure_right = self.get_pressure(x, y);
                let pressure_left = self.get_pressure(x - 1, y);
                self.velocities_x[(x, y)] -= k * (pressure_right - pressure_left);
            }
        }

        // ---- Vertical -------------
        for x in 0..self.velocities_y.width() {
            for y in 0..self.velocities_y.height() {
                let pressure_top = self.get_pressure(x, y);
                let pressure_bottom = self.get_pressure(x, y - 1);
                self.velocities_y[(x, y)] -= k * (pressure_top - pressure_bottom);
            }
        }
    }
}
