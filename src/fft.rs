use super::complex::Complex;

pub fn fft(x: &Vec<Complex>, inverse: bool) -> Vec<Complex> {
    // 事前に w を計算しておく
    let sign = if inverse { 1.0 } else { -1.0 };
    let w: Vec<Complex> = (0..x.len())
        .map(|n| {
            let theta = sign * std::f64::consts::TAU * n as f64 / x.len() as f64;
            Complex::expi(theta)
        })
        .collect();

    // 出力用 (y1 と y2 を交互に使い回す)
    let mut y1 = x.clone();
    let mut y2 = vec![Complex::new(0.0, 0.0); x.len()];

    // FFT
    let mut len = y1.len();
    let mut stride = 1;
    while len >= 2 {
        // len を 2 つの整数の積に分解
        const FACTORS: [usize; 4] = [2, 3, 5, 7];
        let len1 = FACTORS.into_iter().find(|&x| len % x == 0).unwrap_or(len);
        let len2 = len / len1;

        // バタフライ演算
        match len1 {
            2 => {
                // Cooley-Tukey FFT
                // 高速化のため、基数が小さいときは式を展開
                for n2 in 0..len2 {
                    let w = w[(stride * n2) % w.len()];
                    for offset in 0..stride {
                        let a = y1[stride * n2 + offset];
                        let b = y1[stride * (n2 + len2) + offset];
                        y2[stride * (2 * n2) + offset] = a + b;
                        y2[stride * (2 * n2 + 1) + offset] = (a - b) * w;
                    }
                }
            }
            _ if FACTORS.contains(&len1) => {
                // Cooley-Tukey FFT
                y2.iter_mut().for_each(|y| *y = Complex::new(0.0, 0.0));
                for k1 in 0..len1 {
                    for n1 in 0..len1 {
                        for n2 in 0..len2 {
                            let k = len1 * n2 + k1;
                            let n = len2 * n1 + n2;
                            let w = w[(stride * n * k1) % w.len()];
                            for offset in 0..stride {
                                let k = stride * k + offset;
                                let n = stride * n + offset;
                                y2[k] += y1[n] * w;
                            }
                        }
                    }
                }
            }
            _ => {
                // len が小さい素因数を持たない場合は Bluestein's FFT にフォールバック
                assert_eq!(len2, 1);
                let w2: Vec<Complex> = (0..len)
                    .map(|n| {
                        let theta = sign * std::f64::consts::PI * (n * n) as f64 / len as f64;
                        Complex::expi(theta)
                    })
                    .collect();
                let len_ceil = (len * 2 - 1).next_power_of_two();
                for offset in 0..stride {
                    let a: Vec<Complex> = (0..len_ceil)
                        .map(|n| match n {
                            _ if n < len => y1[stride * n + offset] * w2[n],
                            _ => Complex::new(0.0, 0.0),
                        })
                        .collect();
                    let b: Vec<Complex> = (0..len_ceil)
                        .map(|n| match n {
                            _ if n < len => w2[n].conj(),
                            _ if len_ceil - n < len => w2[len_ceil - n].conj(),
                            _ => Complex::new(0.0, 0.0),
                        })
                        .collect();
                    let a_fft = fft(&a, false);
                    let b_fft = fft(&b, false);
                    let y_fft: Vec<Complex> = (a_fft.iter().zip(b_fft.iter()))
                        .map(|(a, b)| *a * *b)
                        .collect();
                    let y = fft(&y_fft, true);
                    for k in 0..len {
                        y2[stride * k + offset] = y[k] * b[k].conj();
                    }
                }
            }
        }

        // y1 と y2 を入れ替えて次のステップへ
        std::mem::swap(&mut y1, &mut y2);
        len = len2;
        stride *= len1;
    }

    // 逆変換の場合は正規化
    if inverse {
        for k in 0..y1.len() {
            y1[k] /= x.len() as f64;
        }
    }

    return y1;
}

pub fn bench() {
    let len_patterns: Vec<usize> = (0..=16).map(|x| 1 << x).collect();
    let count = 30;

    let mut rand = super::rand::XorShift32::default();

    println!("len\tfft");
    for len in len_patterns.into_iter() {
        let x: Vec<Complex> = (0..len)
            .map(|_| Complex::new(rand.f64(), rand.f64()))
            .collect();

        let start = std::time::Instant::now();
        for _ in 0..count {
            let y = fft(&x, false);
            std::hint::black_box(y);
        }
        let fft_dur = start.elapsed().as_secs_f64() / count as f64;

        println!("{}\t{:.9}", len, fft_dur);
    }
}

#[cfg(test)]
mod tests {
    use super::super::rand::*;
    use super::*;

    // テスト用の DFT 実装
    fn dft(x: &Vec<Complex>, inverse: bool) -> Vec<Complex> {
        // 事前に w を計算しておく
        let sign = if inverse { 1.0 } else { -1.0 };
        let w: Vec<Complex> = (0..x.len())
            .map(|n| {
                let theta = sign * std::f64::consts::TAU * n as f64 / x.len() as f64;
                Complex::expi(theta)
            })
            .collect();

        // DFT
        let mut y = vec![Complex::new(0.0, 0.0); x.len()];
        for k in 0..y.len() {
            for n in 0..x.len() {
                y[k] += x[n] * w[(n * k) % w.len()];
            }
            if inverse {
                y[k] /= x.len() as f64;
            }
        }

        return y;
    }

    #[test]
    fn test_dft() {
        // 1 Hz の正弦波を用意
        let x: Vec<Complex> = (0..8)
            .map(|n| Complex::expi(std::f64::consts::TAU * n as f64 / 8.0))
            .collect();
        let expect_dft: Vec<Complex> = vec![0.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            .into_iter()
            .map(|x| Complex::new(x, 0.0))
            .collect();

        // DFT
        let actual_dft = dft(&x, false);
        assert_eq!(actual_dft.len(), expect_dft.len());
        for (a, b) in expect_dft.iter().zip(actual_dft.iter()) {
            assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
            assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
        }

        // iDFT
        let actual_idft = dft(&actual_dft, true);
        assert_eq!(actual_idft.len(), x.len());
        for (a, b) in x.iter().zip(actual_idft.iter()) {
            assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
            assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
        }
    }

    #[test]
    fn test_fft() {
        const LEN_PATTERNS: [usize; 3] = [0, 1, 2 * 2 * 3 * 23];
        for len in LEN_PATTERNS.into_iter() {
            // 乱数列を生成
            let mut rand = XorShift32::default();
            let x: Vec<Complex> = (0..len)
                .map(|_| Complex::new(rand.f64(), rand.f64()))
                .collect();

            // DFT と計算が一致することを確認
            let expect_fft = dft(&x, false);

            // FFT
            let actual_fft = fft(&x, false);
            assert_eq!(actual_fft.len(), expect_fft.len());
            for (a, b) in expect_fft.iter().zip(actual_fft.iter()) {
                assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
                assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
            }

            // iFFT
            let actual_ifft = fft(&actual_fft, true);
            assert_eq!(actual_ifft.len(), x.len());
            for (a, b) in x.iter().zip(actual_ifft.iter()) {
                assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
                assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
            }
        }
    }
}
