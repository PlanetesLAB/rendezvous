pub mod corrector;
pub mod particle;
pub mod simd;
use thiserror::Error;

use rayon::iter::{IndexedParallelIterator, ParallelIterator};

use crate::gravity::{Gravity, IgnoreGravityTerms};
use crate::integrator::{ForceSplit, StepContext, SyncContext, Synchronize};
use crate::particle::{Particles, TestParticleKind, Transformations};
use crate::rendezvous::VariationalOrder;

pub use corrector::Corrector;
pub use particle::RebParticle;

pub struct WHFast {
    recalculate_coordinates_this_time_step: bool,
    coordinates: Coordinates,
    kernel: Kernel,
    corrector: Corrector,
    use_secondary_corrector: bool,
    keep_unsynchoronized: bool,
    is_synchronized: bool,
    safe_mode: SafeMode,
    particles: Particles,
    _time_step_warning: usize,
}

impl WHFast {
    fn init(&mut self, ctx: &mut SyncContext<'_>) -> Result<(), WHFastError> {
        if let Some(var_cfg) = ctx.var_cfg {
            for v in var_cfg {
                if !matches!(v.order, VariationalOrder::First) {
                    return Err(WHFastError::UnsupportedVariationalOrder);
                }
                if v.is_test_particle {
                    return Err(WHFastError::TestParticleVariationalNotSupported);
                }
            }
        }

        if ctx.var_cfg.is_some() && !matches!(self.coordinates, Coordinates::Jacobi) {
            return Err(WHFastError::VariationalJacobiOnly);
        }

        if !matches!(self.kernel, Kernel::Default)
            && !matches!(self.coordinates, Coordinates::Jacobi)
        {
            return Err(WHFastError::NonStandardKernelJacobiOnly);
        }

        if ctx.var_cfg.is_some() && !matches!(self.kernel, Kernel::Default) {
            return Err(WHFastError::VariationalStandardKernelOnly);
        }

        if !matches!(self.corrector, Corrector::None)
            && !matches!(
                self.coordinates,
                Coordinates::Jacobi | Coordinates::Barycentric
            )
        {
            return Err(WHFastError::SymplecticCorrectorJacobiOrBarycentricOnly);
        }

        if self.keep_unsynchoronized && matches!(self.safe_mode, SafeMode::Combine) {
            return Err(WHFastError::InvalidSafeModeCombination);
        }

        if matches!(self.kernel, Kernel::ModifiedKick | Kernel::Lazy) {
            *ctx.gravity = Gravity::Jacobi;
        } else {
            match self.coordinates {
                Coordinates::Jacobi => {
                    *ctx.ignore_gravity_terms = IgnoreGravityTerms::IgnoreWHFastwithJacobi;
                }
                Coordinates::Barycentric => {
                    *ctx.ignore_gravity_terms = IgnoreGravityTerms::IgnoreAll;
                }
                _ => {
                    *ctx.ignore_gravity_terms = IgnoreGravityTerms::IgnoreWHFastwithDHC;
                }
            }
        }

        self.recalculate_coordinates();

        Ok(())
    }

    pub fn recalculate_coordinates(&mut self) {
        self.recalculate_coordinates_this_time_step = true;
    }

    fn kepler_step(&mut self, ctx: &SyncContext<'_>, _dt: f64) {
        let _m0 = ctx.particles[0].m;
        let _n_active = if ctx.particles.are_all_active()
            || matches!(ctx.test_particle_kind, TestParticleKind::Massive)
        {
            ctx.particles.n_real()
        } else {
            ctx.particles.active.len()
        };

        let mut _eta = _m0;

        match self.coordinates {
            Coordinates::Jacobi => {
                self.particles
                    .par_iter_mut()
                    .enumerate()
                    .for_each(|(_i, _p)| {});
            }

            Coordinates::DemocraticHeliocentric => {}
            Coordinates::Whds => {}
            Coordinates::Barycentric => {}
        }
    }

    fn com_step(&mut self, dt: f64) {
        self.particles[0].x += self.particles[0].vx * dt;
        self.particles[0].y += self.particles[0].vy * dt;
        self.particles[0].z += self.particles[0].vz * dt;
    }

    fn apply_corrector(&mut self) {}

    fn apply_secondary_corrector(&mut self) {}
}

impl Synchronize for WHFast {
    fn synchronize(&mut self, ctx: &mut SyncContext<'_>) {
        if self.init(ctx).is_err() {
            return;
        }

        if !self.is_synchronized {
            let n_real = ctx.particles.n_real();
            let n_active = if ctx.particles.are_all_active()
                || matches!(ctx.test_particle_kind, TestParticleKind::Massive)
            {
                n_real
            } else {
                ctx.particles.active.len()
            };

            let sync_particles = if self.keep_unsynchoronized {
                Some(self.particles.clone())
            } else {
                None
            };

            match self.kernel {
                Kernel::Default | Kernel::ModifiedKick => {}
                Kernel::Lazy => {
                    self.kepler_step(ctx, ctx.dt / 2.0);
                    self.com_step(ctx.dt / 2.0);
                }
                Kernel::Composition => {
                    self.kepler_step(ctx, 3.0 * ctx.dt / 8.0);
                    self.com_step(3.0 * ctx.dt / 8.0);
                }
            }

            if self.use_secondary_corrector {
                self.apply_secondary_corrector();
            }

            if !matches!(self.corrector, Corrector::None) {
                self.apply_corrector();
            }

            match self.coordinates {
                Coordinates::Jacobi => {
                    let masses = ctx
                        .particles
                        .active
                        .iter()
                        .map(|p| p.m)
                        .collect::<Vec<f64>>();
                    ctx.particles.active.transform_jacobi_to_inertial_posvel(
                        &self.particles.active,
                        &masses,
                        n_real,
                        n_active,
                    );
                }
                Coordinates::DemocraticHeliocentric => {
                    ctx.particles.active.transform_dhc_to_inertial_posvel(
                        &self.particles.active,
                        n_real,
                        n_active,
                    );
                }
                Coordinates::Whds => {
                    ctx.particles.active.transform_whds_to_inertial_posvel(
                        &self.particles.active,
                        n_real,
                        n_active,
                    );
                }
                Coordinates::Barycentric => {
                    ctx.particles
                        .active
                        .transform_barycentric_to_inertial_posvel(
                            &self.particles.active,
                            n_real,
                            n_active,
                        );
                }
            }

            if let Some(_vcs) = ctx.var_cfg {
                todo!()
            }

            if let Some(sp) = sync_particles {
                self.particles = sp;
            }
        }
    }
}

impl ForceSplit for WHFast {
    fn pre_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }

    fn post_force(&mut self, _ctx: &mut StepContext<'_>) {
        todo!()
    }
}

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
