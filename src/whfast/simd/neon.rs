#![allow(clippy::wildcard_imports)]

use std::arch::aarch64::*;

use super::SimdF64;

#[derive(Clone, Copy)]
#[repr(C)]
pub struct NeonF64x8(float64x2_t, float64x2_t, float64x2_t, float64x2_t);

impl SimdF64 for NeonF64x8 {
    #[inline]
    fn splat(x: f64) -> Self {
        unsafe {
            Self(
                vdupq_n_f64(x),
                vdupq_n_f64(x),
                vdupq_n_f64(x),
                vdupq_n_f64(x),
            )
        }
    }

    #[inline]
    fn load(values: &[f64; 8]) -> Self {
        unsafe {
            let ptr = values.as_ptr();
            Self(
                vld1q_f64(ptr),
                vld1q_f64(ptr.add(2)),
                vld1q_f64(ptr.add(4)),
                vld1q_f64(ptr.add(6)),
            )
        }
    }

    #[inline]
    fn store(self, values: &mut [f64; 8]) {
        unsafe {
            let ptr = values.as_mut_ptr();
            vst1q_f64(ptr, self.0);
            vst1q_f64(ptr.add(2), self.1);
            vst1q_f64(ptr.add(4), self.2);
            vst1q_f64(ptr.add(6), self.3);
        }
    }

    #[inline]
    fn add(self, rhs: Self) -> Self {
        unsafe {
            Self(
                vaddq_f64(self.0, rhs.0),
                vaddq_f64(self.1, rhs.1),
                vaddq_f64(self.2, rhs.2),
                vaddq_f64(self.3, rhs.3),
            )
        }
    }

    #[inline]
    fn sub(self, rhs: Self) -> Self {
        unsafe {
            Self(
                vsubq_f64(self.0, rhs.0),
                vsubq_f64(self.1, rhs.1),
                vsubq_f64(self.2, rhs.2),
                vsubq_f64(self.3, rhs.3),
            )
        }
    }

    #[inline]
    fn mul(self, rhs: Self) -> Self {
        unsafe {
            Self(
                vmulq_f64(self.0, rhs.0),
                vmulq_f64(self.1, rhs.1),
                vmulq_f64(self.2, rhs.2),
                vmulq_f64(self.3, rhs.3),
            )
        }
    }

    #[inline]
    fn div(self, rhs: Self) -> Self {
        unsafe {
            Self(
                vdivq_f64(self.0, rhs.0),
                vdivq_f64(self.1, rhs.1),
                vdivq_f64(self.2, rhs.2),
                vdivq_f64(self.3, rhs.3),
            )
        }
    }

    #[inline]
    fn sqrt(self) -> Self {
        unsafe {
            Self(
                vsqrtq_f64(self.0),
                vsqrtq_f64(self.1),
                vsqrtq_f64(self.2),
                vsqrtq_f64(self.3),
            )
        }
    }

    // load/store omitted here
}
