pub struct MT19937 {
    state: [u32; 624],
    index: usize,
}

impl MT19937 {
    pub fn new(seed: u32) -> MT19937 {
        let mut mt = MT19937 {
            state: [0; 624],
            index: 624,
        };
        mt.state[0] = seed;
        for i in 1..mt.state.len() {
            mt.state[i] = 1812433253u32
                .wrapping_mul(mt.state[i - 1] ^ (mt.state[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        mt
    }
    pub fn u32(&mut self) -> u32 {
        const M: usize = 397;
        const UPPER_MASK: u32 = 0x80000000;
        const LOWER_MASK: u32 = !UPPER_MASK;
        const MAG01: [u32; 2] = [0x00000000, 0x9908b0df];
        let n = self.state.len();
        if self.index >= n {
            self.index = 0;
            for i in 0..n - M {
                let y = (self.state[i] & UPPER_MASK) | (self.state[i + 1] & LOWER_MASK);
                self.state[i] = self.state[i + M] ^ (y >> 1) ^ MAG01[y as usize & 1];
            }
            for i in n - M..n - 1 {
                let y = (self.state[i] & UPPER_MASK) | (self.state[i + 1] & LOWER_MASK);
                self.state[i] = self.state[i + M - n] ^ (y >> 1) ^ MAG01[y as usize & 1];
            }
            let y = (self.state[n - 1] & UPPER_MASK) | (self.state[0] & LOWER_MASK);
            self.state[n - 1] = self.state[M - 1] ^ (y >> 1) ^ MAG01[y as usize & 1];
        }
        let mut y = self.state[self.index];
        y ^= y.wrapping_shr(11);
        y ^= y.wrapping_shl(7) & 0x9d2c5680;
        y ^= y.wrapping_shl(15) & 0xefc60000;
        y ^= y.wrapping_shr(18);
        self.index += 1;
        return y;
    }
    pub fn f64(&mut self) -> f64 {
        self.u32() as f64 / std::u32::MAX as f64
    }
}

impl Default for MT19937 {
    fn default() -> Self {
        MT19937::new(5489)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mt19937() {
        let mut mt = MT19937::default();
        for _ in 0..9999 {
            mt.u32();
        }
        assert_eq!(mt.u32(), 4123659995);
    }
}
