use super::simd::SimdF64;

pub struct RebParticle<V: SimdF64> {
    pub m: V,
    pub x: V,
    pub y: V,
    pub z: V,
    pub vx: V,
    pub vy: V,
    pub vz: V,
}
