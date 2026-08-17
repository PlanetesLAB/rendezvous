use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Coordinates {
    /// Jacobi coordinates (default)
    Jacobi,
    /// Democratic Heliocentric coordinates
    DemocraticHeliocentric,
    /// WHDS coordinates (Hernandez and Dehnen 2017)
    Whds,
    /// Barycentric coordinates
    Barycentric,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kernel {
    Default,
    ModifiedKick,
    Composition,
    Lazy,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SafeMode {
    DriftKickDrift,
    Combine,
}

#[derive(Debug, Error)]
pub enum WHFastError {
    #[error("WHFast/MEGNO only supports first-order variational equations")]
    UnsupportedVariationalOrder,

    #[error("Test particle variational particles are not supported in WHFast")]
    TestParticleVariationalNotSupported,

    #[error("Variational particles require Jacobi coordinates in WHFast")]
    VariationalJacobiOnly,

    #[error("Non-standard kernels require Jacobi coordinates in WHFast")]
    NonStandardKernelJacobiOnly,

    #[error("Variational particles are only supported with the standard kernel in WHFast")]
    VariationalStandardKernelOnly,

    #[error("Symplectic correctors require Jacobi or Barycentric coordinates in WHFast")]
    SymplecticCorrectorJacobiOrBarycentricOnly,

    #[error("Cannot keep unsynchronized particles when using SafeMode::Combine in WHFast")]
    InvalidSafeModeCombination,
}
