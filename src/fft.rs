use super::complex::Complex;

pub fn dft(x: &Vec<Complex>, invererse: bool) -> Vec<Complex> {
    let sign = if invererse { 1.0 } else { -1.0 };
    let mut t = vec![Complex::new(0.0, 0.0); x.len()];

    for ti in 0..t.len() {
        for xi in 0..x.len() {
            let theta = sign * std::f64::consts::TAU * ti as f64 * xi as f64 / x.len() as f64;
            let wn = Complex::expi(theta);
            t[ti] += x[xi] * wn;
        }
    }

    if invererse {
        let t_len = t.len();
        for ti in 0..t_len {
            t[ti] /= t_len as f64;
        }
    }

    t
}

pub fn fft_cooley_tukey(x: &Vec<Complex>, invererse: bool) -> Vec<Complex> {
    /*

      T[0] ---+---------------+---------+-----> T[0]
             ^|              ^|        ^v
      T[1] --||--+-----------||--+-----+--W0--> T[4]
             || ^|           |v ^|
      T[2] --||-||--+--------+--||-W0---+-----> T[2]
             || || ^|           |v     ^v
      T[3] --||-||-||--+--------+--W2--+--W0--> T[6]
             |v || || ^|
      T[4] --+--||-||-||-W0---+---------+-----> T[1]
                |v || ||     ^|        ^v
      T[5] -----+--||-||-W1--||--+-----+--W0--> T[5]
                   |v ||     |v ^|
      T[6] --------+--||-W2--+--||-W0---+-----> T[3]
                      |v        |v     ^v
      T[7] -----------+--W3-----+--W2--+--W0--> T[7]

    */

    let mut t = x.clone();

    fn fft_inner(
        t: &mut [Complex],
        t_len: usize,
        t_offset: usize,
        t_origin: usize,
        invererse: bool,
    ) {
        if t_len > 1 {
            let sign = if invererse { 1.0 } else { -1.0 };
            let t_mid = t_len >> 1;
            for t_i in t_offset..t_mid + t_offset {
                let theta = sign * std::f64::consts::TAU * t_i as f64 / t_len as f64;
                let wn = Complex::expi(theta);
                let a = t[t_i];
                let b = t[t_i + t_mid];
                t[t_i] = a + b;
                t[t_i + t_mid] = (a - b) * wn;
            }
            let stride = t.len() / t_len;
            fft_inner(t, t_mid, t_offset, t_origin, invererse);
            fft_inner(t, t_mid, t_offset + t_mid, t_origin + stride, invererse);
        } else if t_len == 1 && t_offset > t_origin {
            t.swap(t_offset, t_origin);
        }
    }

    let t_len = t.len();
    fft_inner(&mut t, t_len, 0, 0, invererse);

    if invererse {
        for ti in 0..t_len {
            t[ti] /= t_len as f64;
        }
    }

    t
}

pub fn fft_stockham(x: &Vec<Complex>, invererse: bool) -> Vec<Complex> {
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
    let sign = if invererse { 1.0 } else { -1.0 };
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

    if invererse {
        let t_len = t1.len();
        for ti in 0..t_len {
            t1[ti] /= t_len as f64;
        }
    }

    t1
}

pub fn fft_n(x: &Vec<Complex>, invererse: bool) -> Vec<Complex> {
    let sign = if invererse { 1.0 } else { -1.0 };
    let w = (0..x.len())
        .map(|i| {
            let theta = sign * std::f64::consts::TAU * i as f64 / x.len() as f64;
            Complex::expi(theta)
        })
        .collect::<Vec<_>>();
    let mut x1 = x.clone();
    let mut x2 = vec![Complex::new(0.0, 0.0); x1.len()];
    let mut stride = 1;
    let mut len = x1.len();
    while len > 1 {
        let len1 = {
            //const FACTORS: [usize; 6] = [7, 6, 5, 4, 3, 2];
            const FACTORS: [usize; 1] = [2];
            FACTORS.into_iter().find(|f| len % f == 0).unwrap_or(len)
        };
        let len2 = len / len1;
        for offset in 0..stride {
            for n2 in 0..len2 {
                for k1 in 0..len1 {
                    for n1 in 0..len1 {
                        x2[offset + stride * (len1 * n2 + k1)] += x1
                            [offset + stride * (len2 * n1 + n2)]
                            * w[(stride * (n2 + len2 * n1) * k1) % w.len()];
                    }
                }
            }
        }
        stride *= len1;
        len = len2;
        std::mem::swap(&mut x1, &mut x2);
    }
    if invererse {
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

pub fn fft_bluestein(x: &Vec<Complex>, invererse: bool) -> Vec<Complex> {
    let sign = if invererse { 1.0 } else { -1.0 };
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
    if invererse {
        let t_len = t.len();
        for ti in 0..t_len {
            t[ti] /= t_len as f64;
        }
    }
    t
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
        let x = (0..128)
            .map(|_| Complex::new(mt.f64(), mt.f64()))
            .collect::<Vec<Complex>>();
        let expect_fft = dft(&x, false);
        let actual_fft = fft_n(&x, false);
        //let actual_ifft = fft_n(&actual_fft, true);
        assert_eq!(actual_fft.len(), expect_fft.len());
        for (a, b) in expect_fft.iter().zip(actual_fft.iter()) {
            assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
            assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
        }
        //assert_eq!(actual_ifft.len(), x.len());
        //for (a, b) in x.iter().zip(actual_ifft.iter()) {
        //    assert!((a.re - b.re).abs() < 1e-10, "{} != {}", a.re, b.re);
        //    assert!((a.im - b.im).abs() < 1e-10, "{} != {}", a.im, b.im);
        //}
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
        let x = (0..127)
            .map(|_| Complex::new(mt.f64(), mt.f64()))
            .collect::<Vec<Complex>>();
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
