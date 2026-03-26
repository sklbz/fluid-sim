use core::f32;
use std::ops::{Add, Sub};
pub struct Vector2<T> {
    pub x: T,
    pub y: T,
}

impl<T> Vector2<T> {
    pub fn new(x: T, y: T) -> Vector2<T> {
        Vector2 { x, y }
    }
}

impl Add for Vector2<f32> {
    type Output = Vector2<f32>;

    fn add(self, other: Vector2<f32>) -> Vector2<f32> {
        Vector2 {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl Sub for Vector2<f32> {
    type Output = Vector2<f32>;

    fn sub(self, other: Vector2<f32>) -> Vector2<f32> {
        Vector2 {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}
