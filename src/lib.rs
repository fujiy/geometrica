// #![feature(generic_const_exprs)]
// #![allow(incomplete_features)]

pub mod astrodynamics;
pub mod euclidean;
pub mod kinematics;
pub mod linear_space;
pub mod manifold;

pub mod lie;
pub mod linear;

// pub use linear_space::{AffineFrame, AffineSpace, Basis, LinearSpace};
pub use manifold::Scalar;
