#[derive(Debug, Clone, Copy)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

impl Complex {
    pub fn new(re: f64, im: f64) -> Self {
        Complex { re, im }
    }
    pub fn expi(theta: f64) -> Self {
        Complex::new(theta.cos(), theta.sin())
    }
    pub fn conj(&self) -> Self {
        Complex::new(self.re, -self.im)
    }
    pub fn abs(&self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }
}

impl From<f64> for Complex {
    fn from(re: f64) -> Self {
        Complex::new(re, 0.0)
    }
}

pub fn to_complex(vec: &Vec<f64>) -> Vec<Complex> {
    vec.iter().map(|&x| Complex::from(x)).collect()
}
pub fn to_re(vec: &Vec<Complex>) -> Vec<f64> {
    vec.iter().map(|x| x.re).collect()
}
pub fn to_im(vec: &Vec<Complex>) -> Vec<f64> {
    vec.iter().map(|x| x.im).collect()
}
pub fn to_abs(vec: &Vec<Complex>) -> Vec<f64> {
    vec.iter()
        .map(|x| (x.re.powi(2) + x.im.powi(2)).sqrt())
        .collect()
}

use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};
impl Add<Self> for Complex {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Complex::new(self.re + rhs.re, self.im + rhs.im)
    }
}
impl AddAssign<Self> for Complex {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}
impl Sub<Self> for Complex {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Complex::new(self.re - rhs.re, self.im - rhs.im)
    }
}
impl SubAssign<Self> for Complex {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}
impl Mul<Self> for Complex {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Complex::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}
impl MulAssign<f64> for Complex {
    fn mul_assign(&mut self, rhs: f64) {
        *self = *self * rhs;
    }
}
impl Mul<f64> for Complex {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Complex::new(self.re * rhs, self.im * rhs)
    }
}
impl MulAssign<Self> for Complex {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}
impl Div<Complex> for Complex {
    type Output = Self;
    fn div(self, rhs: Complex) -> Self {
        let denom = rhs.re * rhs.re + rhs.im * rhs.im;
        Complex::new(
            (self.re * rhs.re + self.im * rhs.im) / denom,
            (self.im * rhs.re - self.re * rhs.im) / denom,
        )
    }
}
impl DivAssign<Complex> for Complex {
    fn div_assign(&mut self, rhs: Complex) {
        *self = *self / rhs;
    }
}
impl Div<f64> for Complex {
    type Output = Self;
    fn div(self, rhs: f64) -> Self {
        Complex::new(self.re / rhs, self.im / rhs)
    }
}
impl DivAssign<f64> for Complex {
    fn div_assign(&mut self, rhs: f64) {
        *self = *self / rhs;
    }
}
