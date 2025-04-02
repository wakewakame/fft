use math::complex::*;
use math::fft::*;
use math::plot::*;

fn main() {
    let x = (0..8)
        .map(|i| {
            let i = i as f64;
            Complex::expi(1.0 * i * std::f64::consts::TAU / 8.0)
                + Complex::expi(3.0 * i * std::f64::consts::TAU / 8.0)
                + Complex::expi(5.0 * i * std::f64::consts::TAU / 8.0)
                + Complex::expi(7.0 * i * std::f64::consts::TAU / 8.0)
        })
        .collect::<Vec<Complex>>();
    let t = fft(&x, false);
    plot(&t, 40, 10, true);

    math::fft::bench();
}
