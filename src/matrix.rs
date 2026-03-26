use rand::RngExt;
use rand::distr::uniform::SampleUniform;
use std::{
    iter::Sum,
    ops::{Index, IndexMut, RangeInclusive},
};

/// A 2-D grid stored in column-major order.
/// First index is x (horizontal, 0..width), second is y (vertical, 0..height).
/// `data[x * height + y]`
pub struct Matrix<T> {
    width: u32,
    height: u32,
    data: Vec<T>,
}

// ── Bilinear interpolation ────────────────────────────────────────────────────

/// Sample the matrix at continuous position (x, y) using bilinear interpolation.
/// Coordinates are clamped to the valid range so out-of-bounds queries are safe.
pub fn bilinear(m: &Matrix<f32>, x: f32, y: f32) -> f32 {
    let x0 = (x.floor() as i32).clamp(0, m.width as i32 - 1) as u32;
    let y0 = (y.floor() as i32).clamp(0, m.height as i32 - 1) as u32;
    let x1 = (x0 + 1).min(m.width - 1);
    let y1 = (y0 + 1).min(m.height - 1);

    let tx = (x - x.floor()).clamp(0.0, 1.0);
    let ty = (y - y.floor()).clamp(0.0, 1.0);

    let v00 = m[(x0, y0)];
    let v10 = m[(x1, y0)];
    let v01 = m[(x0, y1)];
    let v11 = m[(x1, y1)];

    v00 * (1.0 - tx) * (1.0 - ty) + v10 * tx * (1.0 - ty) + v01 * (1.0 - tx) * ty + v11 * tx * ty
}

// ── Clone ─────────────────────────────────────────────────────────────────────

impl<T: Clone> Clone for Matrix<T> {
    fn clone(&self) -> Self {
        Matrix {
            width: self.width,
            height: self.height,
            data: self.data.clone(),
        }
    }
}

// ── Constructors ──────────────────────────────────────────────────────────────

impl<T: Clone> Matrix<T> {
    pub fn new(width: u32, height: u32, default: T) -> Matrix<T> {
        Matrix {
            width,
            height,
            data: vec![default; (width * height) as usize],
        }
    }
}

// ── Random fill ───────────────────────────────────────────────────────────────

impl<T> Matrix<T>
where
    T: SampleUniform + PartialOrd,
    RangeInclusive<T>: Clone,
{
    pub fn randomize(&mut self, range: RangeInclusive<T>) {
        let mut rng = rand::rng();
        for val in &mut self.data {
            *val = rng.random_range(range.clone());
        }
    }
}

// ── Sum ───────────────────────────────────────────────────────────────────────

impl<T: Copy + Sum> Matrix<T> {
    pub fn sum(&self) -> T {
        self.data.iter().copied().sum()
    }
}

// ── Accessors ─────────────────────────────────────────────────────────────────

impl<T> Matrix<T> {
    #[inline]
    pub fn width(&self) -> u32 {
        self.width
    }
    #[inline]
    pub fn height(&self) -> u32 {
        self.height
    }

    /// Iterate over all (x, y) index pairs.
    pub fn indices(&self) -> MatrixIndices {
        MatrixIndices {
            height: self.height,
            total: self.width * self.height,
            index: 0,
        }
    }
}

// ── Index / IndexMut ──────────────────────────────────────────────────────────

impl<T> Index<(u32, u32)> for Matrix<T> {
    type Output = T;
    #[inline]
    fn index(&self, (x, y): (u32, u32)) -> &T {
        debug_assert!(
            x < self.width && y < self.height,
            "Matrix index ({x},{y}) out of bounds ({},{})",
            self.width,
            self.height
        );
        &self.data[(x * self.height + y) as usize]
    }
}

impl<T> IndexMut<(u32, u32)> for Matrix<T> {
    #[inline]
    fn index_mut(&mut self, (x, y): (u32, u32)) -> &mut T {
        debug_assert!(x < self.width && y < self.height);
        &mut self.data[(x * self.height + y) as usize]
    }
}

// ── Iterator ──────────────────────────────────────────────────────────────────

pub struct MatrixIndices {
    height: u32,
    total: u32,
    index: u32,
}

impl Iterator for MatrixIndices {
    type Item = (u32, u32);
    fn next(&mut self) -> Option<(u32, u32)> {
        if self.index >= self.total {
            return None;
        }
        let x = self.index / self.height;
        let y = self.index % self.height;
        self.index += 1;
        Some((x, y))
    }
}
