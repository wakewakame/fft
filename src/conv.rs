use super::complex::Complex;
use super::fft::fft;

pub fn convolution(f: &Vec<Complex>, g: &Vec<Complex>) -> Vec<Complex> {
    let (f, g) = if f.len() > g.len() { (f, g) } else { (g, f) };
    if g.len() == 0 {
        return f.clone();
    }
    let mut m = vec![Complex::new(0.0, 0.0); f.len() + g.len() - 1];
    for fi in 0..f.len() {
        for gi in 0..g.len() {
            m[fi + gi] += f[fi] * g[gi];
        }
    }
    m
}

pub fn fft_convolution(f: &Vec<Complex>, g: &Vec<Complex>) -> Vec<Complex> {
    if f.len() == 0 || g.len() == 0 {
        return vec![];
    }
    let len = f.len() + g.len() - 1;
    let len_bit_ceil = len.next_power_of_two();
    let mut f = f.clone();
    f.extend(vec![Complex::new(0.0, 0.0); len_bit_ceil - f.len()]);
    let mut g = g.clone();
    g.extend(vec![Complex::new(0.0, 0.0); len_bit_ceil - g.len()]);
    let f_fft = fft(&f, false);
    let g_fft = fft(&g, false);
    let m_fft: Vec<Complex> = f_fft
        .iter()
        .zip(g_fft.iter())
        .map(|(a, b)| *a * *b)
        .collect();
    let mut m = fft(&m_fft, true);
    m.truncate(len);
    m
}

#[cfg(test)]
mod tests {
    use super::super::complex::*;
    use super::*;

    #[test]
    fn test_convolution() {
        let f: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0];
        let g: Vec<f64> = vec![5.0, 6.0, 7.0];
        let expect: Vec<f64> = vec![5.0, 16.0, 34.0, 52.0, 45.0, 28.0];
        let f = to_complex(&f);
        let g = to_complex(&g);
        let expect = to_complex(&expect);
        let actual1 = convolution(&f, &g);
        let actual2 = convolution(&g, &f);
        assert_eq!(actual1.len(), expect.len());
        for (a, b) in expect.iter().zip(actual1.iter()) {
            assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
            assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
        }
        assert_eq!(actual2.len(), expect.len());
        for (a, b) in expect.iter().zip(actual2.iter()) {
            assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
            assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
        }
    }

    #[test]
    fn test_fft_convolution() {
        let f: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0];
        let g: Vec<f64> = vec![5.0, 6.0, 7.0];
        let expect: Vec<f64> = vec![5.0, 16.0, 34.0, 52.0, 45.0, 28.0];
        let f = to_complex(&f);
        let g = to_complex(&g);
        let expect = to_complex(&expect);
        let actual1 = fft_convolution(&f, &g);
        let actual2 = fft_convolution(&g, &f);
        assert_eq!(actual1.len(), expect.len());
        for (a, b) in expect.iter().zip(actual1.iter()) {
            assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
            assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
        }
        assert_eq!(actual2.len(), expect.len());
        for (a, b) in expect.iter().zip(actual2.iter()) {
            assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
            assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
        }
    }
}
