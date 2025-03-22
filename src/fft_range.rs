use super::complex::Complex;

pub fn fft_range(x: &Vec<Complex>, hz_start: f64, hz_end: f64, len: usize) -> Vec<Complex> {
    _ = (x, hz_start, hz_end, len);
    todo!();
}

#[cfg(test)]
mod tests {
    use super::super::rand::*;
    use super::*;

    // テスト用の DFT 実装
    fn dft_range(x: &Vec<Complex>, hz_start: f64, hz_end: f64, len: usize) -> Vec<Complex> {
        // 事前に w を計算しておく
        let hz_range = (hz_end - hz_start) * len as f64 / (len - 1) as f64;
        let w: Vec<Complex> = (0..x.len())
            .map(|n| {
                let theta = -std::f64::consts::TAU * (hz_start + hz_range * n as f64) / len as f64;
                Complex::expi(theta)
            })
            .collect();

        // DFT
        let mut y = vec![Complex::new(0.0, 0.0); len];
        for k in 0..y.len() {
            for n in 0..x.len() {
                y[k] += x[n] * w[(n * k) % w.len()];
            }
        }

        return y;
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
        _ = expect_fft;
        /*
        let actual_fft = fft_range(&x, 3.0, 4.0, 4);
        assert_eq!(actual_fft.len(), expect_fft.len());
        for (a, b) in expect_fft.iter().zip(actual_fft.iter()) {
            assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
            assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
        }
        */
    }
}
