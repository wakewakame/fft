use math::complex::*;
use math::fft::*;
use math::*;

fn main() {
    let mut x = vec![0.0; 8];
    add_sin(&mut x, 1.0, 0.0, 0.25);
    add_sin(&mut x, 3.0, 1.1, 0.125);
    let x = to_complex(&x);
    //let t = dft(&x, false);
    //let t = fft_cooley_tukey(&x, false);
    let t = fft_stockham(&x, false);
    let t_abs = to_abs(&t);
    plot(&t_abs);
    println!("{:?}", t_abs);
}

// 便利関数
pub fn add_sin(vec: &mut Vec<f64>, freq: f64, phase: f64, amp: f64) {
    let len = vec.len() as f64;
    vec.iter_mut().enumerate().for_each(|(i, x)| {
        *x += amp * (std::f64::consts::TAU * freq * i as f64 / len + phase).sin();
    });
}
