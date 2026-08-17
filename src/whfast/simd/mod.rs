mod scalar;

#[cfg(target_arch = "x86_64")]
mod avx512;

#[cfg(target_arch = "aarch64")]
mod neon;

pub use scalar::ScalarF64x8;

#[cfg(target_arch = "x86_64")]
pub use avx512::Avx512F64x8;

#[cfg(target_arch = "aarch64")]
pub use neon::NeonF64x8;

pub trait SimdF64: Copy {
    fn splat(x: f64) -> Self;

    fn load(values: &[f64; 8]) -> Self;
    fn store(self, values: &mut [f64; 8]);

    #[must_use]
    fn add(self, rhs: Self) -> Self;

    #[must_use]
    fn sub(self, rhs: Self) -> Self;

    #[must_use]
    fn mul(self, rhs: Self) -> Self;

    #[must_use]
    fn div(self, rhs: Self) -> Self;

    #[must_use]
    fn sqrt(self) -> Self;
}
