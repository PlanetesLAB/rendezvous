use super::SimdF64;

#[derive(Clone, Copy)]
pub struct ScalarF64x8(pub [f64; 8]);

impl SimdF64 for ScalarF64x8 {
    #[inline]
    fn splat(x: f64) -> Self {
        Self([x; 8])
    }

    #[inline]
    fn load(values: &[f64; 8]) -> Self {
        Self(*values)
    }

    #[inline]
    fn store(self, values: &mut [f64; 8]) {
        *values = self.0;
    }

    #[inline]
    fn add(self, rhs: Self) -> Self {
        let mut result = [0.0; 8];

        for i in 0..8 {
            result[i] = self.0[i] + rhs.0[i];
        }

        Self(result)
    }

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        let mut result = [0.0; 8];

        for i in 0..8 {
            result[i] = self.0[i] - rhs.0[i];
        }

        Self(result)
    }

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        let mut result = [0.0; 8];

        for i in 0..8 {
            result[i] = self.0[i] * rhs.0[i];
        }

        Self(result)
    }

    #[inline]
    fn div(self, rhs: Self) -> Self {
        let mut result = [0.0; 8];

        for i in 0..8 {
            result[i] = self.0[i] / rhs.0[i];
        }

        Self(result)
    }

    #[inline]
    fn sqrt(self) -> Self {
        let mut result = [0.0; 8];

        for i in 0..8 {
            result[i] = self.0[i].sqrt();
        }

        Self(result)
    }
}
