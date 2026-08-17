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
        Self(std::array::from_fn(|i| self.0[i] + rhs.0[i]))
    }

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self(std::array::from_fn(|i| self.0[i] - rhs.0[i]))
    }

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self(std::array::from_fn(|i| self.0[i] * rhs.0[i]))
    }

    #[inline]
    fn div(self, rhs: Self) -> Self {
        Self(std::array::from_fn(|i| self.0[i] / rhs.0[i]))
    }

    #[inline]
    fn sqrt(self) -> Self {
        Self(std::array::from_fn(|i| self.0[i].sqrt()))
    }
}
