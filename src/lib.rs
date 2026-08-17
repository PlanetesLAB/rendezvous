#![doc(
    html_logo_url = "https://raw.githubusercontent.com/PlanetesLAB/documentation/refs/heads/main/logo.jpg"
)]
#![allow(clippy::missing_panics_doc)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::struct_excessive_bools)]
#![allow(clippy::similar_names)]
#![allow(clippy::missing_errors_doc)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::cast_lossless)]

pub mod boundary;
pub mod collision;
pub mod eos;
pub mod gbs;
pub mod gravity;
pub mod ias15;
pub mod integrator;
pub mod janus;
pub mod leapfrog;
pub mod mercurius;
pub mod ode;
pub mod particle;
pub mod rendezvous;
pub mod reversible;
pub mod saba;
pub mod sei;
pub mod trace;
pub mod tree;
pub mod utils;
pub mod whfast;
