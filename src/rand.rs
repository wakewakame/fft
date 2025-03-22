pub struct XorShift32(u32);

impl XorShift32 {
    pub fn new(seed: u32) -> XorShift32 {
        assert_ne!(seed, 0);
        XorShift32(seed)
    }
    pub fn u32(&mut self) -> u32 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 17;
        self.0 ^= self.0 << 5;
        self.0
    }
    pub fn f64(&mut self) -> f64 {
        self.u32() as f64 / std::u32::MAX as f64
    }
}

impl Default for XorShift32 {
    fn default() -> Self {
        XorShift32(1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_xor_shift32() {
        let mut rand = XorShift32::default();
        for _ in 0..9999 {
            rand.u32();
        }
        assert_eq!(rand.u32(), 1799336688);
    }

    #[test]
    #[should_panic]
    fn test_xor_shift32_zero() {
        XorShift32::new(0);
    }
}
