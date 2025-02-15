#[derive(Debug, Clone, Copy)]
pub struct Complex(pub f64, pub f64);

pub fn dft(xf: &Vec<f64>) -> Vec<Complex> {
    // TODO
}

// 便利関数
pub fn add_sin(vec: &mut Vec<f64>, freq: f64, phase: f64, amp: f64) {
    let len = vec.len() as f64;
    vec.iter_mut().enumerate().for_each(|(i, x)| {
        *x += amp * (std::f64::consts::TAU * freq * i as f64 / len + phase).sin();
    });
}

pub fn to_complex(vec: &Vec<f64>) -> Vec<Complex> {
    vec.iter().map(|&x| Complex(x, 0.0)).collect()
}

pub fn to_re(vec: &Vec<Complex>) -> Vec<f64> {
    vec.iter().map(|x| x.0).collect()
}

pub fn to_im(vec: &Vec<Complex>) -> Vec<f64> {
    vec.iter().map(|x| x.1).collect()
}

pub fn to_abs(vec: &Vec<Complex>) -> Vec<f64> {
    vec.iter()
        .map(|x| (x.0.powi(2) + x.1.powi(2)).sqrt())
        .collect()
}

use plotters::prelude::*;
pub fn plot(vec: &Vec<f64>) {
    let root = BitMapBackend::new("./img.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .x_label_area_size(10)
        .y_label_area_size(10)
        .build_cartesian_2d(0f64..vec.len() as f64, -2f64..2f64)
        .unwrap();
    chart.configure_mesh().draw().unwrap();
    chart
        .draw_series(LineSeries::new(
            vec.iter().enumerate().map(|(x, y)| (x as f64, *y)),
            &RED,
        ))
        .unwrap();
    root.present().unwrap();
}
