use nalgebra as na;

use crate::lie::Torsor;
use crate::manifold::Manifold;

pub trait LinearSpace<const N: usize>: Manifold<N> {
    type DualSpace: LinearSpace<N, DualSpace = Self, Field = <Self as Manifold<N>>::Field>;

    fn zero() -> Self {
        Self::_from_raw(na::SVector::zeros())
    }
    fn scale(&self, scalar: <Self as Manifold<N>>::Field) -> Self {
        Self::_from_raw(self._get_raw() * scalar)
    }

    fn pair_with(&self, dual: &Self::DualSpace) -> <Self as Manifold<N>>::Field {
        self._get_raw().dot(&dual._get_raw())
    }
    fn dual_basis(basis: &[Self; N]) -> [Self::DualSpace; N] {
        let raws: [na::SVector<Self::Field, N>; N] = std::array::from_fn(|i| basis[i].raw.clone());
        let matrix = na::SMatrix::<Self::Field, N, N>::from_columns(&raws);
        let inverse = matrix.try_inverse().expect("Basis is singular");

        let inverse_t = inverse.transpose();
        std::array::from_fn(|i| Self::_from_raw(inverse_t.column(i).into_owned()))
    }

    // fn _induced_basis(basis: &Basis<N, Self>) -> [na::SVector<<Self as Manifold<N>>::Field, N>; N];

    fn _get_raw(&self) -> &na::SVector<<Self as Manifold<N>>::Field, N>;
    fn _from_raw(raw: na::SVector<<Self as Manifold<N>>::Field, N>) -> Self;
}

#[macro_export]
macro_rules! impl_vector_ops {
    ($name:ident) => {
        impl<'a, const N: usize, K: Scalar> std::ops::Add<&'a $name<N, K>> for &'a $name<N, K> {
            type Output = $name<N, K>;

            fn add(self, rhs: Self) -> Self::Output {
                Self::_from_raw(self._get_raw() + rhs._get_raw())
            }
        }
    };
}

pub trait InnerProductSpace<const N: usize>: LinearSpace<N> {
    fn dot(&self, other: &Self) -> <Self as Manifold<N>>::Field;
}

pub type AffineSpace<const N: usize, V: LinearSpace<N>> = Torsor<N, V>;
