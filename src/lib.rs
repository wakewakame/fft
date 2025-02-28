pub mod fft;
pub mod rand;

use plotters::prelude::*;
pub fn plot(vec: &Vec<f64>) {
    let root = BitMapBackend::new("./img.png", (640, 480)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .x_label_area_size(10)
        .y_label_area_size(10)
        .build_cartesian_2d(0f64..(vec.len() - 1) as f64, -2f64..2f64)
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
