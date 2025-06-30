// #![feature(generic_const_exprs)]
// #![allow(incomplete_features)]

pub mod manifold;
pub mod linear_space;
pub mod euclidean;
pub mod kinematics;
pub mod astrodynamics;

pub use manifold::{Manifold, LieGroup};

