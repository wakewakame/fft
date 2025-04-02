use super::complex::*;

pub fn plot(v: &Vec<Complex>, width: usize, height: usize, autofit: bool) {
    let len = (v.len() - 1).max(1);
    let (w2, h3) = (width * 2, height * 3);
    let y: Vec<f64> = (0..w2)
        .into_iter()
        .map(|x| {
            let i = len * x / w2;
            let p = (len * x) as f64 / w2 as f64 - i as f64;
            let y = v[i] * (1.0 - p) as f64 + v[i + 1] * p as f64;
            return y.abs();
        })
        .collect();
    let (y_min, y_range) = if autofit {
        let y_min = y.iter().fold(0.0, |a: f64, b: &f64| a.min(*b));
        let y_max = y.iter().fold(0.0, |a: f64, b: &f64| a.max(*b));
        let y_range = y_max - y_min;
        (y_min, y_range)
    } else {
        (0.0, 1.0)
    };
    let y = y
        .iter()
        .map(|y| (h3 as f64) * (1.0 - (y - y_min) / y_range))
        .collect::<Vec<f64>>();
    let mut result = "".to_string();
    for h in 0..height {
        for w in 0..width {
            let mut dots = 0x2800;
            for w2 in 0..2 {
                for h2 in 0..3 {
                    if y[w * 2 + w2] <= (3 * h + h2) as f64 {
                        dots |= 1 << (w2 * 3 + h2);
                    }
                }
            }
            result.push(std::char::from_u32(dots).unwrap());
        }
        result.push('\n');
    }
    println!("{}", result);
}
