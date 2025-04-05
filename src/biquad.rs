pub struct BiQuad {
    sample_rate: f64,
    a0: f64,
    a1: f64,
    a2: f64,
    b0: f64,
    b1: f64,
    b2: f64,
    i1: f64,
    i2: f64,
    o1: f64,
    o2: f64,
}

// TODO: 複数のフィルタを直列に繋げられるようにしたい
// TODO: フィルタの種類を増やしたい
// TODO: カットオフ周波数の正弦波を入力したときのテストを追加したい
impl BiQuad {
    pub fn new(sample_rate: f64) -> Self {
        Self {
            sample_rate,
            a0: 1.0,
            a1: 0.0,
            a2: 0.0,
            b0: 1.0,
            b1: 0.0,
            b2: 0.0,
            i1: 0.0,
            i2: 0.0,
            o1: 0.0,
            o2: 0.0,
        }
    }

    pub fn low_pass(&mut self, cutoff: f64, q: f64) {
        let omega = std::f64::consts::TAU * cutoff / self.sample_rate;
        let alpha = omega.sin() / (2.0 * q);
        self.a0 = 1.0 + alpha;
        self.a1 = -2.0 * omega.cos();
        self.a2 = 1.0 - alpha;
        self.b0 = (1.0 - omega.cos()) / 2.0;
        self.b1 = 1.0 - omega.cos();
        self.b2 = (1.0 - omega.cos()) / 2.0;
    }

    pub fn high_pass(&mut self, cutoff: f64, q: f64) {
        let omega = std::f64::consts::TAU * cutoff / self.sample_rate;
        let alpha = omega.sin() / (2.0 * q);

        self.a0 = 1.0 + alpha;
        self.a1 = -2.0 * omega.cos();
        self.a2 = 1.0 - alpha;
        self.b0 = (1.0 + omega.cos()) / 2.0;
        self.b1 = -(1.0 + omega.cos());
        self.b2 = (1.0 + omega.cos()) / 2.0;
    }

    pub fn band_pass(&mut self, cutoff: f64, bw: f64) {
        let omega = std::f64::consts::TAU * cutoff / self.sample_rate;
        let alpha = omega.sin() * 2.0f64.ln().sinh() / 2.0 * bw * omega / omega.sin();

        self.a0 = 1.0 + alpha;
        self.a1 = -2.0 * omega.cos();
        self.a2 = 1.0 - alpha;
        self.b0 = alpha;
        self.b1 = 0.0;
        self.b2 = -alpha;
    }

    pub fn process(&mut self, input: f64) -> f64 {
        let output = (self.b0 * input + self.b1 * self.i1 + self.b2 * self.i2
            - self.a1 * self.o1
            - self.a2 * self.o2)
            / self.a0;

        self.i2 = self.i1;
        self.i1 = input;
        self.o2 = self.o1;
        self.o1 = output;

        output
    }

    // 周波数特性の配列を返す
    pub fn freq_response(&self, start_hz: f64, end_hz: f64, num_points: usize) -> Vec<f64> {
        use super::complex::Complex;
        let mut result = Vec::with_capacity(num_points);
        for i in 0..num_points {
            let hz = start_hz + (end_hz - start_hz) * i as f64 / num_points as f64;
            let omega = std::f64::consts::TAU * hz / self.sample_rate;
            let z = Complex::expi(omega);
            let h = (Complex::new(self.b0, 0.0) + z * self.b1 + z * z * self.b2)
                / (Complex::new(self.a0, 0.0) + z * self.a1 + z * z * self.a2);
            result.push(h.abs());
        }
        result
    }
}
