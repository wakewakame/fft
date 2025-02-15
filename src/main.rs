use fft::*;

fn main() {
    let mut x = vec![0.0; 8];
    add_sin(&mut x, 1.0, 0.0, 0.25);
    add_sin(&mut x, 3.0, 1.1, 0.125);
    let x = to_complex(&x);
    //let t = dft(&x, false);
    //let t = fft_cooley_tukey(&x, false);
    let t = fft_stockham(&x, false);
    let t_abs = to_abs(&t);
    fft::plot(&t_abs);
    println!("{:?}", t_abs);
}
