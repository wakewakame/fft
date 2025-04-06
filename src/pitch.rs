pub struct Pitch {}

impl Pitch {
    pub fn new() -> Self {
        Pitch {}
    }

    pub fn pitch_shift(&self, input: &Vec<f64>, shift: f64) -> Vec<f64> {
        std::hint::black_box(shift);
        input.clone()
    }
}
