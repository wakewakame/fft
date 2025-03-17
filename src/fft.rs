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

pub fn fft_inplace(x: &Vec<Complex>, inverse: bool) -> Vec<Complex> {
    fn butterfly2(w: &[Complex], x: &mut [Complex], len: usize, len2: usize) {
        for offset in (0..x.len()).step_by(len) {
            for n2 in offset..(offset + len2) {
                let a = x[n2];
                let b = x[n2 + len2];
                let w = w[(n2 * w.len() / len) % w.len()];
                x[n2] = a + b;
                x[n2 + len2] = (a - b) * w;
            }
        }
    }

    fn butterfly_n(w: &[Complex], x: &mut [Complex], len: usize, len1: usize, len2: usize) {
        for offset in (0..x.len()).step_by(len) {
            for n2 in 0..len2 {
                let mut new_x2 = vec![Complex::new(0.0, 0.0); len1];
                for k1 in 0..len1 {
                    for n1 in 0..len1 {
                        let n = len2 * n1 + n2;
                        let w = w[(n * k1 * w.len() / len) % w.len()];
                        new_x2[k1] += x[offset + n] * w;
                    }
                }
                for k1 in 0..len1 {
                    let k = len2 * k1 + n2;
                    x[offset + k] = new_x2[k1];
                }
            }
        }
    }

    // 事前に w を計算しておく
    let sign = if inverse { 1.0 } else { -1.0 };
    let w: Vec<Complex> = (0..x.len())
        .map(|i| {
            let theta = sign * std::f64::consts::TAU * i as f64 / x.len() as f64;
            Complex::expi(theta)
        })
        .collect();

    let mut x = x.clone();

    let mut len = x.len();
    while len > 1 {
        // len を 2 つの整数の積に分解
        let len1 = {
            const FACTORS: [usize; 6] = [13, 11, 7, 5, 3, 2];
            FACTORS.into_iter().find(|f| len % f == 0).unwrap_or(len)
        };
        let len2 = len / len1;

        match len1 {
            2 => butterfly2(&w, &mut x, len, len2),
            3..=7 => butterfly_n(&w, &mut x, len, len1, len2),
            _ => todo!(), // _ => bluestein(&w, &x, len),
        }

        len = len2;
    }

    // 逆変換の場合は正規化
    if inverse {
        let len = x.len();
        for n in 0..len {
            x[n] /= len as f64;
        }
    }

    // TODO: 2 の冪乗以外にも対応する
    bit_reverse(x.len(), &mut x);

    x
}

fn bit_reverse(n: usize, x: &mut [Complex]) {
    let mut i = 0;
    for j in 1..n - 1 {
        let mut k = n >> 1;
        i ^= k;
        while k > i {
            k >>= 1;
            i ^= k;
        }
        if i < j {
            x.swap(i, j);
        }
    }
}

pub fn fft_cooley_tukey(x: &Vec<Complex>, inverse: bool) -> Vec<Complex> {
    // 事前に w を計算しておく
    let sign = if inverse { 1.0 } else { -1.0 };
    let w = (0..x.len())
        .map(|i| {
            let theta = sign * std::f64::consts::TAU * i as f64 / x.len() as f64;
            Complex::expi(theta)
        })
        .collect::<Vec<_>>();

    let mut x = x.clone();

    fn fft_inner(
        w: &Vec<Complex>,
        x: &mut [Complex],
        len: usize,
        offset: usize,
        origin: usize,
        inverse: bool,
    ) {
        if len > 1 {
            let mid = len >> 1;
            for i in offset..mid + offset {
                let wn = w[i * w.len() / len % w.len()];
                let a = x[i];
                let b = x[i + mid];
                x[i] = a + b;
                x[i + mid] = (a - b) * wn;
            }
            let stride = x.len() / len;
            fft_inner(w, x, mid, offset, origin, inverse);
            fft_inner(w, x, mid, offset + mid, origin + stride, inverse);
        } else if len == 1 && offset > origin {
            x.swap(offset, origin);
        }
    }

    let len = x.len();
    fft_inner(&w, &mut x, len, 0, 0, inverse);

    if inverse {
        for ti in 0..len {
            x[ti] /= len as f64;
        }
    }

    x
}

pub fn fft_stockham(x: &Vec<Complex>, inverse: bool) -> Vec<Complex> {
    /*

      T1[0] --+-------------> T2[0] | T2[0] --+-------------> T1[0] | T1[0] --+-------------> T2[0]
             ^|                     |        ^|                     |        ^|
      T1[1] -||--+----------> T2[2] | T2[1] -||--+----------> T1[1] | T1[1] -||--+----------> T2[1]
             || ^|                  |        || ^|                  |        || ^|
      T1[2] -||-||--+-------> T2[4] | T2[2] -||-||--+-------> T1[4] | T1[2] -||-||--+-------> T2[2]
             || || ^|               |        || || ^|               |        || || ^|
      T1[3] -||-||-||--+----> T2[6] | T2[3] -||-||-||--+----> T1[5] | T1[3] -||-||-||--+----> T2[3]
             |v || || ^|            |        |v || || ^|            |        |v || || ^|
      T1[4] -+--||-||-||-W0-> T2[1] | T2[4] -+--||-||-||-W0-> T1[2] | T1[4] -+--||-||-||-W0-> T2[4]
                |v || ||            |           |v || ||            |           |v || ||
      T1[5] ----+--||-||-W1-> T2[3] | T2[5] ----+--||-||-W0-> T1[3] | T1[5] ----+--||-||-W0-> T2[5]
                   |v ||            |              |v ||            |              |v ||
      T1[6] -------+--||-W2-> T2[5] | T2[6] -------+--||-W2-> T1[6] | T1[6] -------+--||-W0-> T2[6]
                      |v            |                 |v            |                 |v
      T1[7] ----------+--W3-> T2[7] | T2[7] ----------+--W2-> T1[7] | T1[7] ----------+--W0-> T2[7]

    */

    let mut t1 = x.clone();
    let mut t2 = vec![Complex::new(0.0, 0.0); t1.len()];
    let sign = if inverse { 1.0 } else { -1.0 };
    let w = (0..t1.len())
        .map(|i| {
            let theta = sign * std::f64::consts::TAU * i as f64 / t1.len() as f64;
            Complex::expi(theta)
        })
        .collect::<Vec<_>>();
    let t_mid = t1.len() >> 1;
    let mut stride = 1;
    while stride < t1.len() {
        for w_i in (0..t_mid).step_by(stride) {
            let wn = w[w_i];
            for s_i in 0..stride {
                let t1_i = w_i + s_i;
                let t2_i = w_i * 2 + s_i;
                let a = t1[t1_i];
                let b = t1[t1_i + t_mid];
                t2[t2_i] = a + b;
                t2[t2_i + stride] = (a - b) * wn;
            }
        }
        stride = stride << 1;
        std::mem::swap(&mut t1, &mut t2);
    }

    if inverse {
        let t_len = t1.len();
        for ti in 0..t_len {
            t1[ti] /= t_len as f64;
        }
    }

    t1
}

pub fn fft_n(x: &Vec<Complex>, inverse: bool) -> Vec<Complex> {
    // 事前に w を計算しておく
    let sign = if inverse { 1.0 } else { -1.0 };
    let w = (0..x.len())
        .map(|i| {
            let theta = sign * std::f64::consts::TAU * i as f64 / x.len() as f64;
            Complex::expi(theta)
        })
        .collect::<Vec<_>>();

    // 出力用 (x1 と x2 を交互に使い回す)
    let mut x1 = x.clone();
    let mut x2 = vec![Complex::new(0.0, 0.0); x1.len()];

    // バタフライ演算
    let mut stride = 1;
    let mut len = x.len();
    while len > 1 {
        // x2 を 0 で初期化
        for k in 0..x2.len() {
            x2[k] = Complex::new(0.0, 0.0);
        }

        // len を 2 つの整数の積に分解
        let len1 = {
            const FACTORS: [usize; 6] = [13, 11, 7, 5, 3, 2];
            FACTORS.into_iter().find(|f| len % f == 0).unwrap_or(len)
        };
        let len2 = len / len1;

        // バタフライ演算
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
    let f_fft = fft_cooley_tukey(&f, false);
    let g_fft = fft_cooley_tukey(&g, false);
    let m_fft = f_fft
        .iter()
        .zip(g_fft.iter())
        .map(|(a, b)| *a * *b)
        .collect::<Vec<Complex>>();
    let mut m = fft_cooley_tukey(&m_fft, true);
    m.truncate(len);
    m
}

pub fn fft_bluestein(x: &Vec<Complex>, inverse: bool) -> Vec<Complex> {
    let sign = if inverse { 1.0 } else { -1.0 };
    let w = (0..x.len())
        .map(|i| {
            let theta = sign * std::f64::consts::PI * (i * i) as f64 / x.len() as f64;
            Complex::expi(theta)
        })
        .collect::<Vec<_>>();
    let len = (x.len() * 2 - 1).next_power_of_two();
    let a = (0..len)
        .map(|i| {
            if i < x.len() {
                x[i] * w[i]
            } else {
                Complex::new(0.0, 0.0)
            }
        })
        .collect::<Vec<_>>();
    let b = (0..len)
        .map(|i| {
            if i < x.len() {
                w[i].conj()
            } else if (len - i) < x.len() {
                w[len - i].conj()
            } else {
                Complex::new(0.0, 0.0)
            }
        })
        .collect::<Vec<_>>();
    let a_fft = fft_cooley_tukey(&a, false);
    let b_fft = fft_cooley_tukey(&b, false);
    let t_fft = a_fft
        .iter()
        .zip(b_fft.iter())
        .map(|(a, b)| *a * *b)
        .collect::<Vec<Complex>>();
    let mut t = fft_cooley_tukey(&t_fft, true);
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
    let mut mt = super::rand::MT19937::default();

    println!("len,\tinplace,\tcooley_tukey,\tstockham,\tn,\tbluestein");
    for len in len_patterns.into_iter() {
        let x: Vec<Complex> = (0..len).map(|_| Complex::new(mt.f64(), mt.f64())).collect();

        let start = std::time::Instant::now();
        let _ = fft_inplace(&x, false);
        let inplace_dur = start.elapsed().as_secs_f64();

        let start = std::time::Instant::now();
        let _ = fft_cooley_tukey(&x, false);
        let cooley_tukey_dur = start.elapsed().as_secs_f64();

        let start = std::time::Instant::now();
        let _ = fft_stockham(&x, false);
        let stockham_dur = start.elapsed().as_secs_f64();

        let start = std::time::Instant::now();
        let _ = fft_n(&x, false);
        let n_dur = start.elapsed().as_secs_f64();

        let start = std::time::Instant::now();
        let _ = fft_bluestein(&x, false);
        let bluestein_dur = start.elapsed().as_secs_f64();

        println!(
            "{},\t{:.9},\t{:.9},\t{:.9},\t{:.9},\t{:.9}",
            len, inplace_dur, cooley_tukey_dur, stockham_dur, n_dur, bluestein_dur
        );
    }
}

#[cfg(test)]
mod tests {
    use super::super::complex::*;
    use super::super::rand::*;
    use super::*;

    #[test]
    fn test_dft() {
        let x = (0..8)
            .map(|x| Complex::new((std::f64::consts::TAU * x as f64 / 8.0).cos(), 0.0))
            .collect::<Vec<Complex>>();
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
        let mut mt = MT19937::default();
        //const LEN_PATTERNS: [usize; 3] = [128, 2 * 2 * 2 * 3 * 3 * 7, 127];
        const LEN_PATTERNS: [usize; 1] = [128];
        for len in LEN_PATTERNS.into_iter() {
            let x: Vec<Complex> = (0..len).map(|_| Complex::new(mt.f64(), mt.f64())).collect();
            let expect_fft = dft(&x, false);
            let actual_fft = fft_inplace(&x, false);
            let actual_ifft = fft_inplace(&actual_fft, true);
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
    fn test_fft_cooley_tukey() {
        let mut mt = MT19937::default();
        let x = (0..128)
            .map(|_| Complex::new(mt.f64(), mt.f64()))
            .collect::<Vec<Complex>>();
        let expect_fft = dft(&x, false);
        let actual_fft = fft_cooley_tukey(&x, false);
        let actual_ifft = fft_cooley_tukey(&actual_fft, true);
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
    fn test_fft_stockham() {
        let mut mt = MT19937::default();
        let x = (0..128)
            .map(|_| Complex::new(mt.f64(), mt.f64()))
            .collect::<Vec<Complex>>();
        let expect_fft = dft(&x, false);
        let actual_fft = fft_stockham(&x, false);
        let actual_ifft = fft_stockham(&actual_fft, true);
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
    fn test_fft_n() {
        let mut mt = MT19937::default();
        const LEN_PATTERNS: [usize; 3] = [128, 2 * 2 * 2 * 3 * 3 * 7, 127];
        for len in LEN_PATTERNS.into_iter() {
            let x: Vec<Complex> = (0..len).map(|_| Complex::new(mt.f64(), mt.f64())).collect();
            let expect_fft = dft(&x, false);
            let actual_fft = fft_n(&x, false);
            let actual_ifft = fft_n(&actual_fft, true);
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
        let mut mt = MT19937::default();
        const LEN_PATTERNS: [usize; 3] = [128, 2 * 2 * 2 * 3 * 3 * 7, 127];
        for len in LEN_PATTERNS.into_iter() {
            let x: Vec<Complex> = (0..len).map(|_| Complex::new(mt.f64(), mt.f64())).collect();
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

    #[test]
    fn test_bit_reverse() {
        let mut x: Vec<Complex> = vec![0, 1, 2, 3, 4, 5, 6, 7]
            .into_iter()
            .map(|x| Complex::new(x as f64, 0.0))
            .collect();
        let expect: Vec<Complex> = vec![0, 4, 2, 6, 1, 5, 3, 7]
            .into_iter()
            .map(|x| Complex::new(x as f64, 0.0))
            .collect();
        bit_reverse(x.len(), &mut x);
        for (a, b) in expect.iter().zip(x.iter()) {
            assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
            assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
        }
    }
}
