use nalgebra as na;

use super::super::space::LinearSpace;
use crate::lie::LieGroup;
use crate::manifold::Manifold;

pub type GL<const N: usize, V> = GeneralLinearGroup<N, V>;

pub struct GeneralLinearGroup<const N: usize, V: LinearSpace<N>> {
    matrix: na::SMatrix<V::Field, N, N>,
}

macro_rules! impl_general_linear_group {
    ($N:literal, $D:literal) => {
        impl<V: LinearSpace<$N>> Manifold<$D> for GeneralLinearGroup<$N, V> {
            type Field = V::Field;
        }

        impl<V: LinearSpace<$N>> LieGroup<$D> for GeneralLinearGroup<$N, V> {
            fn identity() -> Self {
                Self {
                    matrix: na::SMatrix::<V::Field, $N, $N>::identity(),
                }
            }

            fn multiply(&self, other: &Self) -> Self {
                Self {
                    matrix: self.matrix * other.matrix,
                }
            }

            fn inverse(&self) -> Self {
                Self {
                    matrix: self.matrix.try_inverse().unwrap(),
                }
            }
        }
    };
}

impl_general_linear_group!(1, 1);
impl_general_linear_group!(2, 4);
impl_general_linear_group!(3, 9);
impl_general_linear_group!(4, 16);
