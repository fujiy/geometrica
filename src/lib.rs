// #![feature(generic_const_exprs)]
// #![allow(incomplete_features)]

#[macro_use]
pub mod linear_space;
pub mod affine_space;
pub mod astrodynamics;
pub mod euclidean;
pub mod kinematics;
pub mod lie_group;

// pub use linear_space::{AffineFrame, AffineSpace, Basis, LinearSpace};
pub use lie_group::Scalar;
