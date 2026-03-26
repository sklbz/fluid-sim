use std::ops::{Add, Mul, Sub};

#[derive(Clone, Copy)]
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
    fn add(self, o: Vector2<f32>) -> Vector2<f32> {
        Vector2::new(self.x + o.x, self.y + o.y)
    }
}

impl Sub for Vector2<f32> {
    type Output = Vector2<f32>;
    fn sub(self, o: Vector2<f32>) -> Vector2<f32> {
        Vector2::new(self.x - o.x, self.y - o.y)
    }
}

impl Mul<f32> for Vector2<f32> {
    type Output = Vector2<f32>;
    fn mul(self, s: f32) -> Vector2<f32> {
        Vector2::new(self.x * s, self.y * s)
    }
}
