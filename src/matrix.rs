use rand::RngExt;
use rand::distr::uniform::SampleUniform;
use std::ops::{Index, IndexMut, RangeInclusive};

pub struct Matrix<T> {
    rows: u32,
    cols: u32,
    data: Vec<T>,
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
