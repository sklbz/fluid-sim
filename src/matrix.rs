use rand::RngExt;
use rand::distr::uniform::SampleUniform;
use std::{
    iter::Sum,
    ops::{Index, IndexMut, RangeInclusive},
};

pub struct Matrix<T> {
    rows: u32,
    cols: u32,
    data: Vec<T>,
}

pub fn bilinear(m: &Matrix<f32>, x: f32, y: f32) -> f32 {
    let x0 = (x.floor() as i32).clamp(0, m.width() as i32 - 1) as u32;
    let y0 = (y.floor() as i32).clamp(0, m.height() as i32 - 1) as u32;
    let x1 = (x0 + 1).min(m.width() - 1);
    let y1 = (y0 + 1).min(m.height() - 1);

    let tx = x - x.floor();
    let ty = y - y.floor();

    let v00 = m[(x0, y0)];
    let v10 = m[(x1, y0)];
    let v01 = m[(x0, y1)];
    let v11 = m[(x1, y1)];

    v00 * (1.0 - tx) * (1.0 - ty) + v10 * tx * (1.0 - ty) + v01 * (1.0 - tx) * ty + v11 * tx * ty
}

impl<T: Clone> Clone for Matrix<T> {
    fn clone(&self) -> Self {
        Matrix {
            rows: self.rows,
            cols: self.cols,
            data: self.data.clone(),
        }
    }
}

impl<T> Matrix<T>
where
    T: SampleUniform + PartialOrd,
    RangeInclusive<T>: Clone,
{
    pub fn random(rows: u32, cols: u32, range: RangeInclusive<T>) -> Matrix<T> {
        let mut rng = rand::rng();
        let mut data = Vec::with_capacity(rows as usize * cols as usize);
        for _ in 0..rows * cols {
            data.push(rng.random_range(range.clone()));
        }
        Matrix { rows, cols, data }
    }

    pub fn randomize(&mut self, range: RangeInclusive<T>) {
        let mut rng = rand::rng();
        for i in 0..self.rows * self.cols {
            self.data[i as usize] = rng.random_range(range.clone());
        }
    }
}

impl<T> Matrix<T>
where
    T: Sum + std::marker::Copy,
{
    pub fn sum(&self) -> T {
        self.data.iter().copied().sum()
    }
}

impl<T: Clone> Matrix<T> {
    pub fn new(rows: u32, cols: u32, default: T) -> Matrix<T> {
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

    pub fn width(&self) -> u32 {
        self.cols
    }

    pub fn height(&self) -> u32 {
        self.rows
    }

    pub fn iter(&'_ self) -> MatrixIter<'_, T> {
        MatrixIter {
            matrix: self,
            index: 0,
        }
    }
    pub fn indices(&self) -> MatrixIndices {
        MatrixIndices {
            cols: self.cols,
            total: self.rows * self.cols,
            index: 0,
        }
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

// The iterator struct holds a reference to the matrix and tracks position
pub struct MatrixIter<'a, T> {
    matrix: &'a Matrix<T>,
    index: u32,
}

impl<'a, T> Iterator for MatrixIter<'a, T> {
    // Each item yields the (row, col) coordinate alongside a reference to the value
    type Item = (u32, u32, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        let total = self.matrix.rows * self.matrix.cols;
        if self.index >= total {
            return None;
        }
        let row = self.index / self.matrix.cols;
        let col = self.index % self.matrix.cols;
        let value = self.matrix.get(row, col);
        self.index += 1;
        Some((row, col, value))
    }
}

pub struct MatrixIndices {
    cols: u32,
    total: u32,
    index: u32,
}

impl Iterator for MatrixIndices {
    type Item = (u32, u32);

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.total {
            return None;
        }
        let row = self.index / self.cols;
        let col = self.index % self.cols;
        self.index += 1;
        Some((row, col))
    }
}
