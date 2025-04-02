use super::complex::Complex;
use super::fft::fft;

pub fn fft_range(x: &Vec<Complex>, hz_start: f64, hz_end: f64, len: usize) -> Vec<Complex> {
    let hz_step = (hz_end - hz_start) / (len - 1) as f64;
    let len_ceil = (x.len() + len - 1).next_power_of_two();
    let a: Vec<Complex> = (0..len_ceil)
        .map(|n| match n {
            _ if n < x.len() => {
                let nf = n as f64;
                let hz = hz_start + hz_step * nf / 2.0;
                let theta = -std::f64::consts::TAU * hz * nf / x.len() as f64;
                let w = Complex::expi(theta);
                x[n] * w
            }
            _ => Complex::new(0.0, 0.0),
        })
        .collect();
    let b: Vec<Complex> = (0..len_ceil)
        .map(|n| match n {
            _ if n < len => {
                let nf = n as f64;
                let hz = hz_step * nf / 2.0;
                let theta = -std::f64::consts::TAU * hz * nf / x.len() as f64;
                let w = Complex::expi(theta);
                w.conj()
            }
            _ if len_ceil - n < x.len() => {
                let nf = (len_ceil - n) as f64;
                let hz = hz_step * nf / 2.0;
                let theta = -std::f64::consts::TAU * hz * nf / x.len() as f64;
                let w = Complex::expi(theta);
                w.conj()
            }
            _ => Complex::new(0.0, 0.0),
        })
        .collect();
    let a_fft = fft(&a, false);
    let b_fft = fft(&b, false);
    let y_fft: Vec<Complex> = (a_fft.iter().zip(b_fft.iter()))
        .map(|(a, b)| *a * *b)
        .collect();
    let mut y = fft(&y_fft, true);
    y.truncate(len);
    for k in 0..y.len() {
        y[k] *= b[k].conj();
    }
    return y;
}

#[cfg(test)]
mod tests {
    use super::super::rand::*;
    use super::*;

    // テスト用の DFT 実装
    fn dft_range(x: &Vec<Complex>, hz_start: f64, hz_end: f64, len: usize) -> Vec<Complex> {
        let hz_step = (hz_end - hz_start) / (len - 1) as f64;

        // DFT
        let mut y = vec![Complex::new(0.0, 0.0); len];
        for k in 0..y.len() {
            let hz = hz_start + hz_step * k as f64;
            for n in 0..x.len() {
                let theta = -std::f64::consts::TAU * hz * n as f64 / x.len() as f64;
                let w = Complex::expi(theta);
                y[k] += x[n] * w;
            }
        }

        return y;
    }

    #[test]
    fn test_dft() {
        // 2 Hz の正弦波を用意
        let x: Vec<Complex> = (0..8)
            .map(|n| Complex::expi(2.0 * std::f64::consts::TAU * n as f64 / 8.0))
            .collect();
        let expect_dft: Vec<Complex> = vec![0.0, 8.0, 0.0, 0.0]
            .into_iter()
            .map(|x| Complex::new(x, 0.0))
            .collect();

        // DFT
        let actual_dft = dft_range(&x, 1.0, 4.0, 4);
        assert_eq!(actual_dft.len(), expect_dft.len());
        for (a, b) in expect_dft.iter().zip(actual_dft.iter()) {
            assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
            assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
        }
    }

    #[test]
    fn test_fft_range() {
        // 大小さまざまな素因数を持つ長さの配列に対してテスト
        let len = 2 * 2 * 3 * 23;

        // 乱数列を生成
        let mut rand = XorShift32::default();
        let x: Vec<Complex> = (0..len)
            .map(|_| Complex::new(rand.f64(), rand.f64()))
            .collect();

        // DFT と計算が一致することを確認
        let expect_fft = dft_range(&x, 3.0, 4.0, 4);

        // FFT
        let actual_fft = fft_range(&x, 3.0, 4.0, 4);
        assert_eq!(actual_fft.len(), expect_fft.len());
        for (a, b) in expect_fft.iter().zip(actual_fft.iter()) {
            assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
            assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
        }
    }
}
