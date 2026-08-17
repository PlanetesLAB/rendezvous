use std::arch::x86_64::*;

use super::SimdF64;

#[derive(Clone, Copy)]
#[repr(transparent)]
pub struct Avx512F64x8(pub __m512d);

impl SimdF64 for Avx512F64x8 {
    #[inline]
    fn splat(x: f64) -> Self {
        unsafe { Self(_mm512_set1_pd(x)) }
    }

    #[inline]
    fn load(values: &[f64; 8]) -> Self {
        unsafe { Self(_mm512_loadu_pd(values.as_ptr())) }
    }

    #[inline]
    fn store(self, values: &mut [f64; 8]) {
        unsafe { _mm512_storeu_pd(values.as_mut_ptr(), self.0) }
    }

    #[inline]
    fn add(self, rhs: Self) -> Self {
        unsafe { Self(_mm512_add_pd(self.0, rhs.0)) }
    }

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        unsafe { Self(_mm512_sub_pd(self.0, rhs.0)) }
    }

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        unsafe { Self(_mm512_mul_pd(self.0, rhs.0)) }
    }

    #[inline]
    fn div(self, rhs: Self) -> Self {
        unsafe { Self(_mm512_div_pd(self.0, rhs.0)) }
    }

    #[inline]
    fn sqrt(self) -> Self {
        unsafe { Self(_mm512_sqrt_pd(self.0)) }
    }
}
