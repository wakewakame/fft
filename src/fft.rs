#[derive(Debug, Clone, Copy)]
pub struct Complex(pub f64, pub f64);

pub fn dft(x: &Vec<Complex>, invererse: bool) -> Vec<Complex> {
    let sign = if invererse { 1.0 } else { -1.0 };
    let t_len = x.len();
    let mut t = vec![Complex(0.0, 0.0); t_len];

    for ti in 0..t.len() {
        for xi in 0..x.len() {
            let theta = sign * std::f64::consts::TAU * ti as f64 * xi as f64 / x.len() as f64;
            let wn = Complex(theta.cos(), theta.sin());
            t[ti].0 += x[xi].0 * wn.0 - x[xi].1 * wn.1;
            t[ti].1 += x[xi].0 * wn.1 + x[xi].1 * wn.0;
        }
    }

    if invererse {
        t.iter_mut().for_each(|x| {
            x.0 /= t_len as f64;
            x.1 /= t_len as f64;
        });
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
    let t_len = t.len();

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
                let wn = Complex(theta.cos(), theta.sin());
                let a = t[t_i];
                let b = t[t_i + t_mid];
                let ba = Complex(a.0 - b.0, a.1 - b.1);
                t[t_i] = Complex(a.0 + b.0, a.1 + b.1);
                t[t_i + t_mid] = Complex(ba.0 * wn.0 - ba.1 * wn.1, ba.0 * wn.1 + ba.1 * wn.0);
            }
            let stride = t.len() / t_len;
            fft_inner(t, t_mid, t_offset, t_origin, invererse);
            fft_inner(t, t_mid, t_offset + t_mid, t_origin + stride, invererse);
        } else if t_len == 1 && t_offset > t_origin {
            t.swap(t_offset, t_origin);
        }
    }
    fft_inner(&mut t, t_len, 0, 0, invererse);

    if invererse {
        t.iter_mut().for_each(|x| {
            x.0 /= t_len as f64;
            x.1 /= t_len as f64;
        });
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

    let t_len = x.len();
    let mut t1 = x.clone();
    let mut t2 = vec![Complex(0.0, 0.0); t_len];
    let sign = if invererse { 1.0 } else { -1.0 };
    let w = (0..t_len)
        .map(|i| {
            let theta = sign * std::f64::consts::TAU * i as f64 / t_len as f64;
            Complex(theta.cos(), theta.sin())
        })
        .collect::<Vec<_>>();
    let t_mid = t_len >> 1;
    let mut stride = 1;
    while stride < t_len {
        for w_i in (0..t_mid).step_by(stride) {
            let wn = w[w_i];
            for s_i in 0..stride {
                let t1_i = w_i + s_i;
                let t2_i = w_i * 2 + s_i;
                let a = t1[t1_i];
                let b = t1[t1_i + t_mid];
                let ba = Complex(a.0 - b.0, a.1 - b.1);
                t2[t2_i] = Complex(a.0 + b.0, a.1 + b.1);
                t2[t2_i + stride] = Complex(ba.0 * wn.0 - ba.1 * wn.1, ba.0 * wn.1 + ba.1 * wn.0);
            }
        }
        stride = stride << 1;
        std::mem::swap(&mut t1, &mut t2);
    }

    if invererse {
        t1.iter_mut().for_each(|x| {
            x.0 /= t_len as f64;
            x.1 /= t_len as f64;
        });
    }

    t1
}

pub fn fft_blestein(_: &Vec<f64>, _: &Vec<f64>) -> Vec<f64> {
    todo!();
}

pub fn convolution(_: &Vec<f64>, _: &Vec<f64>) -> Vec<f64> {
    todo!();
}

// 便利関数
pub fn add_sin(vec: &mut Vec<f64>, freq: f64, phase: f64, amp: f64) {
    let len = vec.len() as f64;
    vec.iter_mut().enumerate().for_each(|(i, x)| {
        *x += amp * (std::f64::consts::TAU * freq * i as f64 / len + phase).sin();
    });
}

pub fn to_complex(vec: &Vec<f64>) -> Vec<Complex> {
    vec.iter().map(|&x| Complex(x, 0.0)).collect()
}

pub fn to_re(vec: &Vec<Complex>) -> Vec<f64> {
    vec.iter().map(|x| x.0).collect()
}

pub fn to_im(vec: &Vec<Complex>) -> Vec<f64> {
    vec.iter().map(|x| x.1).collect()
}

pub fn to_abs(vec: &Vec<Complex>) -> Vec<f64> {
    vec.iter()
        .map(|x| (x.0.powi(2) + x.1.powi(2)).sqrt())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dft() {
        let x = (0..8)
            .map(|x| Complex((std::f64::consts::TAU * x as f64 / 8.0).cos(), 0.0))
            .collect::<Vec<Complex>>();
        let y = dft(&x, false);
        let z = dft(&y, true);
        assert!((y[0].0 - 0.0).abs() < 1e-10);
        assert!((y[1].0 - 4.0).abs() < 1e-10);
        assert!((y[2].0 - 0.0).abs() < 1e-10);
        assert!((y[3].0 - 0.0).abs() < 1e-10);
        assert!((y[4].0 - 0.0).abs() < 1e-10);
        assert!((y[5].0 - 0.0).abs() < 1e-10);
        assert!((y[6].0 - 0.0).abs() < 1e-10);
        assert!((y[7].0 - 4.0).abs() < 1e-10);
        assert!((y[0].1 - 0.0).abs() < 1e-10);
        assert!((y[1].1 - 0.0).abs() < 1e-10);
        assert!((y[2].1 - 0.0).abs() < 1e-10);
        assert!((y[3].1 - 0.0).abs() < 1e-10);
        assert!((y[4].1 - 0.0).abs() < 1e-10);
        assert!((y[5].1 - 0.0).abs() < 1e-10);
        assert!((y[6].1 - 0.0).abs() < 1e-10);
        assert!((y[7].1 - 0.0).abs() < 1e-10);
        for (a, b) in x.iter().zip(z.iter()) {
            assert!((a.0 - b.0).abs() < 1e-10);
            assert!((a.1 - b.1).abs() < 1e-10);
        }
    }

    #[test]
    fn test_fft_cooley_tukey() {
        let mut mt = super::super::rand::MT19937::default();
        let x = (0..128)
            .map(|_| Complex(mt.f64(), mt.f64()))
            .collect::<Vec<Complex>>();
        let expect_fft = dft(&x, false);
        let actual_fft = fft_cooley_tukey(&x, false);
        for (a, b) in expect_fft.iter().zip(actual_fft.iter()) {
            assert!((a.0 - b.0).abs() < 1e-10);
            assert!((a.1 - b.1).abs() < 1e-10);
        }
        let actual_ifft = fft_cooley_tukey(&actual_fft, true);
        for (a, b) in x.iter().zip(actual_ifft.iter()) {
            assert!((a.0 - b.0).abs() < 1e-10);
            assert!((a.1 - b.1).abs() < 1e-10);
        }
    }

    #[test]
    fn test_fft_stockham() {
        let mut mt = super::super::rand::MT19937::default();
        let x = (0..128)
            .map(|_| Complex(mt.f64(), mt.f64()))
            .collect::<Vec<Complex>>();
        let expect_fft = dft(&x, false);
        let actual_fft = fft_stockham(&x, false);
        for (a, b) in expect_fft.iter().zip(actual_fft.iter()) {
            assert!((a.0 - b.0).abs() < 1e-10);
            assert!((a.1 - b.1).abs() < 1e-10);
        }
        let actual_ifft = fft_cooley_tukey(&actual_fft, true);
        for (a, b) in x.iter().zip(actual_ifft.iter()) {
            assert!((a.0 - b.0).abs() < 1e-10);
            assert!((a.1 - b.1).abs() < 1e-10);
        }
    }
}
