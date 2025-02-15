use fft::*;

fn main() {
    let mut xf = vec![0.0; 128];
    add_sin(&mut xf, 1.0, 0.0, 1.0);
    add_sin(&mut xf, 2.0, 0.0, 1.0);
    println!("{:?}", xf);
    fft::plot(&xf);
}
