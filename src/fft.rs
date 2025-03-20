use super::complex::Complex;

pub fn dft(x: &Vec<Complex>, inverse: bool) -> Vec<Complex> {
    let sign = if inverse { 1.0 } else { -1.0 };
    let mut t = vec![Complex::new(0.0, 0.0); x.len()];

    for ti in 0..t.len() {
        for xi in 0..x.len() {
            let theta = sign * std::f64::consts::TAU * ti as f64 * xi as f64 / x.len() as f64;
            let wn = Complex::expi(theta);
            t[ti] += x[xi] * wn;
        }
    }

    if inverse {
        let t_len = t.len();
        for ti in 0..t_len {
            t[ti] /= t_len as f64;
        }
    }

    t
}

pub fn fft(x: &Vec<Complex>, inverse: bool) -> Vec<Complex> {
    // 事前に w を計算しておく
    let sign = if inverse { 1.0 } else { -1.0 };
    let w: Vec<Complex> = (0..x.len())
        .map(|i| {
            let theta = sign * std::f64::consts::TAU * i as f64 / x.len() as f64;
            Complex::expi(theta)
        })
        .collect();

    // 出力用 (x1 と x2 を交互に使い回す)
    let mut x1 = x.clone();
    let mut x2 = vec![Complex::new(0.0, 0.0); x1.len()];

    // FFT
    let mut stride = 1;
    let mut len = x.len();
    while len > 1 {
        // len を 2 つの整数の積に分解
        let len1 = {
            const FACTORS: [usize; 6] = [2, 3, 5, 7, 11, 13];
            FACTORS.into_iter().find(|f| len % f == 0).unwrap_or(len)
        };
        let len2 = len / len1;

        // バタフライ演算
        match len1 {
            // 高速化のため、基数が小さいときは式を展開
            2 => {
                for n2 in 0..len2 {
                    let w = w[(stride * n2) % w.len()];
                    for offset in 0..stride {
                        let a = x1[stride * n2 + offset];
                        let b = x1[stride * (n2 + len2) + offset];
                        x2[stride * (2 * n2) + offset] = a + b;
                        x2[stride * (2 * n2 + 1) + offset] = (a - b) * w;
                    }
                }
            }
            _ => {
                x2.iter_mut().for_each(|x| *x = Complex::new(0.0, 0.0));
                for k1 in 0..len1 {
                    for n1 in 0..len1 {
                        for n2 in 0..len2 {
                            let k = len1 * n2 + k1;
                            let n = len2 * n1 + n2;
                            let w = w[(stride * n * k1) % w.len()];
                            for offset in 0..stride {
                                x2[stride * k + offset] += x1[stride * n + offset] * w;
                            }
                        }
                    }
                }
            }
        }

        // x1 と x2 を入れ替えて次のステップへ
        std::mem::swap(&mut x1, &mut x2);
        stride *= len1;
        len = len2;
    }

    // 逆変換の場合は正規化
    if inverse {
        for n in 0..x.len() {
            x1[n] /= x.len() as f64;
        }
    }

    x1
}

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

pub fn fft_bluestein(x: &Vec<Complex>, inverse: bool) -> Vec<Complex> {
    let sign = if inverse { 1.0 } else { -1.0 };
    let w: Vec<Complex> = (0..x.len())
        .map(|i| {
            let theta = sign * std::f64::consts::PI * (i * i) as f64 / x.len() as f64;
            Complex::expi(theta)
        })
        .collect();
    let len = (x.len() * 2 - 1).next_power_of_two();
    let a: Vec<Complex> = (0..len)
        .map(|i| match i {
            _ if i < x.len() => x[i] * w[i],
            _ => Complex::new(0.0, 0.0),
        })
        .collect();
    let b: Vec<Complex> = (0..len)
        .map(|i| match i {
            _ if i < x.len() => w[i].conj(),
            _ if len - i < x.len() => w[len - i].conj(),
            _ => Complex::new(0.0, 0.0),
        })
        .collect();
    let a_fft = fft(&a, false);
    let b_fft = fft(&b, false);
    let t_fft: Vec<Complex> = (0..len).map(|i| a_fft[i] * b_fft[i]).collect();
    let mut t = fft(&t_fft, true);
    t.truncate(x.len());
    for ti in 0..t.len() {
        t[ti] *= b[ti].conj();
    }
    if inverse {
        let t_len = t.len();
        for ti in 0..t_len {
            t[ti] /= t_len as f64;
        }
    }
    t
}

pub fn bench() {
    let len_patterns: Vec<usize> = (0..=16).map(|x| 1 << x).collect();
    let mut rand = super::rand::XorShift32::default();

    println!("len,\tfft,\tbluestein");
    for len in len_patterns.into_iter() {
        let x: Vec<Complex> = (0..len)
            .map(|_| Complex::new(rand.f64(), rand.f64()))
            .collect();

        let start = std::time::Instant::now();
        let _ = fft(&x, false);
        let fft_dur = start.elapsed().as_secs_f64();

        let start = std::time::Instant::now();
        let _ = fft_bluestein(&x, false);
        let bluestein_dur = start.elapsed().as_secs_f64();

        println!("{},\t{:.9},\t{:.9}", len, fft_dur, bluestein_dur);
    }
}

#[cfg(test)]
mod tests {
    use super::super::complex::*;
    use super::super::rand::*;
    use super::*;

    #[test]
    fn test_dft() {
        let x: Vec<Complex> = (0..8)
            .map(|x| Complex::new((std::f64::consts::TAU * x as f64 / 8.0).cos(), 0.0))
            .collect();
        let expect_fft = vec![0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0];
        let expect_fft = to_complex(&expect_fft);
        let actual_fft = dft(&x, false);
        let actual_ifft = dft(&actual_fft, true);
        assert_eq!(actual_fft.len(), expect_fft.len());
        for (a, b) in expect_fft.iter().zip(actual_fft.iter()) {
            assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
            assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
        }
        assert_eq!(actual_ifft.len(), x.len());
        for (a, b) in x.iter().zip(actual_ifft.iter()) {
            assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
            assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
        }
    }

    #[test]
    fn test_fft() {
        let mut rand = XorShift32::default();
        const LEN_PATTERNS: [usize; 3] = [128, 2 * 2 * 2 * 3 * 3 * 7, 127];
        for len in LEN_PATTERNS.into_iter() {
            let x: Vec<Complex> = (0..len)
                .map(|_| Complex::new(rand.f64(), rand.f64()))
                .collect();
            let expect_fft = dft(&x, false);
            let actual_fft = fft(&x, false);
            let actual_ifft = fft(&actual_fft, true);
            assert_eq!(actual_fft.len(), expect_fft.len());
            for (a, b) in expect_fft.iter().zip(actual_fft.iter()) {
                assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
                assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
            }
            assert_eq!(actual_ifft.len(), x.len());
            for (a, b) in x.iter().zip(actual_ifft.iter()) {
                assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
                assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
            }
        }
    }

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

    #[test]
    fn test_fft_bluestein() {
        let mut rand = XorShift32::default();
        const LEN_PATTERNS: [usize; 3] = [128, 2 * 2 * 2 * 3 * 3 * 7, 127];
        for len in LEN_PATTERNS.into_iter() {
            let x: Vec<Complex> = (0..len)
                .map(|_| Complex::new(rand.f64(), rand.f64()))
                .collect();
            let expect_fft = dft(&x, false);
            let actual_fft = fft_bluestein(&x, false);
            let actual_ifft = fft_bluestein(&actual_fft, true);
            assert_eq!(actual_fft.len(), expect_fft.len());
            for (a, b) in expect_fft.iter().zip(actual_fft.iter()) {
                assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
                assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
            }
            assert_eq!(actual_ifft.len(), x.len());
            for (a, b) in x.iter().zip(actual_ifft.iter()) {
                assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
                assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
            }
        }
    }
}
