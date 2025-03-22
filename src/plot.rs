use super::complex::*;

pub fn plot(v: &Vec<Complex>, width: usize, height: usize) {
    let len = (v.len() - 1).max(1);
    let y: Vec<f64> = (0..width)
        .into_iter()
        .map(|x| {
            let i = len * x / width;
            let p = (len * x) as f64 / width as f64 - i as f64;
            let y = v[i] * (1.0 - p) as f64 + v[i + 1] * p as f64;
            return y.abs();
        })
        .collect();
    let y_min = y.iter().fold(0.0, |a: f64, b: &f64| a.min(*b));
    let y_max = y.iter().fold(0.0, |a: f64, b: &f64| a.max(*b));
    let y_range = y_max - y_min;
    let y = y
        .iter()
        .map(|y| (height as f64) * (1.0 - (y - y_min) / y_range))
        .collect::<Vec<f64>>();
    for h in 0..height {
        for w in 0..width {
            if y[w] <= h as f64 {
                print!(".");
            } else {
                print!(" ");
            }
        }
        println!();
    }
}
