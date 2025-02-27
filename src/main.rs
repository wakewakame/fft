use math::*;

fn main() {
    let mut x = vec![0.0; 8];
    fft::add_sin(&mut x, 1.0, 0.0, 0.25);
    fft::add_sin(&mut x, 3.0, 1.1, 0.125);
    let x = fft::to_complex(&x);
    //let t = fft::dft(&x, false);
    //let t = fft::fft_cooley_tukey(&x, false);
    let t = fft::fft_stockham(&x, false);
    let t_abs = fft::to_abs(&t);
    plot(&t_abs);
    println!("{:?}", t_abs);
}
