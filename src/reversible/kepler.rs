use std::f64::consts::PI;
use std::mem;

use super::stiefel;
use crate::particle::Particle;

pub const NMAX_QUARTIC: usize = 64;
pub const NMAX_NEWTON: usize = 32;

pub struct KeplerSystem<'a> {
    pub mass: f64,
    pub r0: f64,
    pub r0i: f64,
    pub v2: f64,
    pub beta: f64,
    pub eta0: f64,
    pub zeta0: f64,
    pub gs: [f64; 6],
    pub p: &'a mut Particle,
}

impl<'a> KeplerSystem<'a> {
    pub fn new(p: &'a mut Particle, mass: f64) -> Self {
        let r0 = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
        let r0i = 1.0 / r0;
        let v2 = p.vx * p.vx + p.vy * p.vy + p.vz * p.vz;
        let beta = 2.0 * mass * r0i - v2;
        let eta0 = p.x * p.vx + p.y * p.vy + p.z * p.vz;
        let zeta0 = mass - beta * r0;

        Self {
            mass,
            r0,
            r0i,
            v2,
            beta,
            eta0,
            zeta0,
            gs: [0.0; 6],
            p,
        }
    }

    #[allow(clippy::too_many_lines)]
    pub fn solve(&mut self, dt: f64) {
        struct BetaVals {
            x_per_period: f64,
            inv_period: f64,
        }

        let betvals = if self.beta > 0.0 {
            Some(BetaVals {
                x_per_period: 2.0 * PI * self.beta.sqrt(),
                inv_period: self.beta.sqrt() * self.beta / (2.0 * PI * self.mass),
            })
        } else {
            None
        };

        let mut nwarns = 0;

        let mut x = if let Some(bv) = &betvals {
            if dt.abs() * bv.inv_period > 1.0 && nwarns == 0 {
                nwarns += 1;
                eprintln!(
                    "WARNING: WHFast is having convergence issues because the time step is \
                     comparable to or larger than the orbital period (dt * inv_period = {}). \
                     Consider reducing the time step.",
                    dt.abs() * bv.inv_period
                );
            }

            let dtr0i = dt * self.r0i;
            dtr0i * (1.0 - dtr0i * self.eta0 * 0.5 * self.r0i)
        } else {
            0.0
        };

        let mut converged = false;
        let mut old_x = x;

        stiefel::gs3(&mut self.gs, self.beta, x);

        let val = self.eta0 * self.gs[1] + self.zeta0 * self.gs[2];

        let mut ri = 1.0 / (self.r0 + val);
        x = ri * (x * val - self.eta0 * self.gs[2] - self.zeta0 * self.gs[3] + dt);

        if let Some(bv) = &betvals {
            if (x - old_x).abs() > 0.01 * bv.x_per_period {
                // Quartic solver
                // Linear initial guess
                x = self.beta * dt / self.mass;
                let mut prev_x = [0.0; NMAX_QUARTIC + 1];
                for n_lag in 1..NMAX_QUARTIC {
                    stiefel::gs3(&mut self.gs, self.beta, x);
                    let f = self.r0 * x + self.eta0 * self.gs[2] + self.zeta0 * self.gs[3] - dt;
                    let fp = self.r0 + self.eta0 * self.gs[1] + self.zeta0 * self.gs[2];
                    let fpp = self.eta0 * self.gs[0] + self.zeta0 * self.gs[1];
                    let denom = fp + (16.0 * fp * fp - 20.0 * f * fpp).abs().sqrt();
                    let x_new = (x * denom - 5.0 * f) / denom;

                    for &pxi in prev_x.iter().take(n_lag).skip(1) {
                        if (x_new - pxi).abs() < 1e-12 {
                            converged = true;
                            break;
                        }
                    }

                    if converged {
                        break;
                    }

                    prev_x[n_lag] = x_new;
                    x = x_new;
                }
                let val = self.eta0 * self.gs[1] + self.zeta0 * self.gs[2];
                ri = 1.0 / (self.r0 + val);
            } else {
                // Newton's method

                for _ in 1..NMAX_NEWTON {
                    let old_x2 = old_x;
                    old_x = x;
                    stiefel::gs3(&mut self.gs, self.beta, x);
                    let val = self.eta0 * self.gs[1] + self.zeta0 * self.gs[2];
                    let ri = 1.0 / (self.r0 + val);
                    x = ri * (x * val - self.eta0 * self.gs[2] - self.zeta0 * self.gs[3] + dt);
                    if (x - old_x).abs() < 1e-12 || (x - old_x2).abs() < 1e-12 {
                        converged = true;
                        break;
                    }
                }
            }
        }

        if !converged {
            let (mut xmin, mut xmax) = if let Some(bv) = &betvals {
                let xmin = bv.x_per_period * (dt * bv.inv_period).floor();
                let xmax = xmin + bv.x_per_period;
                (xmin, xmax)
            } else {
                // Hyperbolic
                let h2 = self.r0 * self.r0 * self.v2 - self.eta0 * self.eta0;
                let q = h2
                    / self.mass
                    / (1.0 + (1.0 - h2 * self.beta / (self.mass * self.mass)).sqrt());
                let vq = (h2.sqrt() / q).copysign(dt);
                let mut xmin = dt / ((vq * dt).abs() + self.r0);
                let mut xmax = dt / q;
                if dt < 0.0 {
                    mem::swap(&mut xmin, &mut xmax);
                }
                (xmin, xmax)
            };
            x = 0.5 * (xmin + xmax);

            loop {
                stiefel::gs3(&mut self.gs, self.beta, x);
                let s = self.r0 * x + self.eta0 * self.gs[2] + self.zeta0 * self.gs[3] - dt;
                if s >= 0.0 {
                    xmax = x;
                } else {
                    xmin = x;
                }
                x = 0.5 * (xmin + xmax);

                if (xmax - xmin).abs() > ((xmax + xmin) * 1e-15).abs() {
                    break;
                }
            }
            let val = self.eta0 * self.gs[1] + self.zeta0 * self.gs[2];
            ri = 1.0 / (self.r0 + val);
        }

        if ri.is_nan() {
            ri = 0.0;
            self.gs[1] = 0.0;
            self.gs[2] = 0.0;
            self.gs[3] = 0.0;
        }

        let _f = -self.mass * self.gs[2] * self.r0i;
        let _g = dt - self.mass * self.gs[3];
        let _fd = -self.mass * self.gs[1] * self.r0i * ri;
        let _gd = -self.mass * self.gs[2] * ri;
    }
}

pub struct KeplerResult {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
    pub nwarns: usize,
    pub val: f64,
}

// fn kepler_solver_variations(
//     ctx: &SyncContext<'_>,
//     dp1: &mut Particle,
//     kpl: &mut KeplerSystem,
//     i: usize,
//     x: f64,
//     mass: f64,
//     dt: f64,
// ) {
//     // Variations
//     if let Some(vcs) = ctx.var_cfg {
//         for vc in vcs.iter() {
//             stiefel::gs(&mut kpl.gs, kpl.beta, x);
//             let dr0 = (dp1.x * x0 + dp1.y * y0 + dp1.z * z0) * r0i;
//             let dbeta =
//                 -2.0 * mass * dr0 * r0i * r0i - 2.0 * (dp1.vx * vx0 + dp1.vy * vy0 + dp1.vz * vz0);
//             let deta0 =
//                 dp1.x * vx0 + x0 * dp1.vx + dp1.y * vy0 + y0 * dp1.vy + dp1.z * vz0 + z0 * dp1.vz;
//             let dzeta0 = -beta * dr0 - r0 * dbeta;
//             let g3_beta = 0.5 * (3.0 * gs[5] - x * gs[4]);
//             let g2_beta = 0.5 * (2.0 * gs[4] - x * gs[3]);
//             let g1_beta = 0.5 * (gs[3] - x * gs[2]);
//             let tbeta = eta0 * g2_beta + zeta0 * g3_beta;
//             let dx = -ri * (x * dr0 + gs[2] * deta0 + gs[3] * dzeta0 + tbeta * dbeta);
//             let dg1 = gs[0] * dx + g1_beta * dbeta;
//             let dg2 = gs[1] * dx + g2_beta * dbeta;
//             let dg3 = gs[2] * dx + g3_beta * dbeta;
//             let dr = dr0 + gs[1] * deta0 + gs[2] * dzeta0 + eta0 * dg1 + zeta0 * dg2;
//             let df = mass * gs[2] * dr0 * r0i * r0i - mass * dg2 * r0i;
//             let dg = -mass * dg3;
//             let dfd = -mass * dg1 * r0i * ri + mass * gs[1] * (dr0 * r0i + dr * ri) * r0i * ri;
//             let dgd = -mass * dg2 * ri + mass * gs[2] * dr * ri * ri;
//
//             dp1.x += f * dp1.x + g * dp1.vx + df * x0 + dg * vx0;
//             dp1.y += f * dp1.y + g * dp1.vy + df * y0 + dg * vy0;
//             dp1.z += f * dp1.z + g * dp1.vz + df * z0 + dg * vz0;
//
//             dp1.vx += fd * dp1.x + gd * dp1.vx + dfd * x0 + dgd * vx0;
//             dp1.vy += fd * dp1.y + gd * dp1.vy + dfd * y0 + dgd * vy0;
//             dp1.vz += fd * dp1.z + gd * dp1.vz + dfd * z0 + dgd * vz0;
//         }
//     }
// }
