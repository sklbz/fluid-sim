use std::ops::{Index, IndexMut};

fn main() {
    let mut grid = FluidGrid::new((10, 10), 1.0);
    grid.velocities_x.set(0, 0, 1.0);
    grid.velocities_y[(0, 0)] = 1.0;
    for x in 0..grid.cell_count.0 {
        for y in 0..grid.cell_count.1 {
            print!("{:^3}", grid.divergence_at_cell(x, y));
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
        let top_velocity = self.velocities_y.get(x, y + 1);
        let left_velocity = self.velocities_x.get(x, y);
        let right_velocity = self.velocities_x.get(x + 1, y);
        let bottom_velocity = self.velocities_y.get(x, y);

        // rate of change of fluid velocity in either axis
        let gradient_x = (right_velocity - left_velocity) / self.cell_size;
        let gradient_y = (top_velocity - bottom_velocity) / self.cell_size;

        gradient_x + gradient_y
    }

    pub fn pressure_solve_cell(&mut self, x: u32, y: u32) {
        let top_pressure = self.velocities_y.get(x, y + 1);
        let left_pressure = self.velocities_x.get(x, y);
        let right_pressure = self.velocities_x.get(x + 1, y);
        let bottom_pressure = self.velocities_y.get(x, y);

        let top_velocity = self.velocities_y.get(x, y + 1);
        let left_velocity = self.velocities_x.get(x, y);
        let right_velocity = self.velocities_x.get(x + 1, y);
        let bottom_velocity = self.velocities_y.get(x, y);

        let pressure_sum = top_pressure + left_pressure + right_pressure + bottom_pressure;
        let delta_velocity_sum = right_velocity - left_velocity + top_velocity - bottom_velocity;

        self.pressure_map[(x, y)] = (pressure_sum
            - self.density * self.cell_size * delta_velocity_sum / self.time_step)
            / 4.0;
    }
}

struct Matrix<T> {
    rows: u32,
    cols: u32,
    data: Vec<T>,
}

impl<T: Clone> Matrix<T> {
    fn new(rows: u32, cols: u32, default: T) -> Matrix<T> {
        Matrix {
            rows,
            cols,
            data: vec![default; rows as usize * cols as usize],
        }
    }
}

impl<T> Matrix<T> {
    fn get(&self, row: u32, col: u32) -> &T {
        &self.data[(row * self.cols + col) as usize]
    }

    fn set(&mut self, row: u32, col: u32, value: T) {
        self.data[(row * self.cols + col) as usize] = value;
    }
}

impl<T> Index<(u32, u32)> for Matrix<T> {
    type Output = T;

    fn index(&self, (row, col): (u32, u32)) -> &T {
        self.get(row, col)
    }
}

impl<T> IndexMut<(u32, u32)> for Matrix<T> {
    fn index_mut(&mut self, (row, col): (u32, u32)) -> &mut T {
        &mut self.data[(row * self.cols + col) as usize]
    }
}
