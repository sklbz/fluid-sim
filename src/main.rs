pub mod matrix;

use std::{thread::sleep, time::Duration};

use matrix::Matrix;

fn main() {
    let mut grid = FluidGrid::new((10, 10), 1.0);
    grid.velocities_x.randomize(-10.0..=10.0);
    grid.velocities_y.randomize(-10.0..=10.0);

    display_divergence(&grid);
    display_pressure(&grid);

    for _ in 0..10 {
        for _ in 0..50 {
            println!("\x1b[2J\x1b[H");
            display_divergence(&grid);
            display_pressure(&grid);

            grid.pressure_solve();
            sleep(Duration::from_millis(100));
        }
        grid.update_velocities();
        sleep(Duration::from_millis(200));
    }
}
fn display_divergence(grid: &FluidGrid) {
    println!("Divergence:");
    for x in 0..grid.cell_count.0 {
        for y in 0..grid.cell_count.1 {
            let div = grid.divergence_at_cell(x, y);
            if div.abs() > 1.0 {
                if div > 0.0 {
                    print!("\x1b[31m");
                } else {
                    print!("\x1b[34m");
                }
            } else {
                print!("\x1b[32m");
            }
            print!("{:^5.1}", div);
            print!("\x1b[0m");
        }
        println!();
    }
}

fn display_pressure(grid: &FluidGrid) {
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
                let pressure_left = x
                    .checked_sub(1)
                    .map(|px| self.get_pressure(px, y))
                    .unwrap_or(0.0);
                self.velocities_x[(x, y)] -= k * (pressure_right - pressure_left);
            }
        }

        // ---- Vertical -------------
        for x in 0..self.velocities_y.width() - 1 {
            // Pourquoi ? Je ne sais pas
            for y in 0..self.velocities_y.height() {
                let pressure_top = self.get_pressure(x, y);
                let pressure_bottom = y
                    .checked_sub(1)
                    .map(|py| self.get_pressure(x, py))
                    .unwrap_or(0.0);
                self.velocities_y[(x, y)] -= k * (pressure_top - pressure_bottom);
            }
        }
    }
}
