use nalgebra::DimName;
use nalgebra::{DimMul, DimProd, OMatrix};

use crate::lie::{GroupAction, LieGroup};
use crate::linear::space::{Allocator, DefaultAllocator, LinearSpace};
use crate::manifold::Manifold;

pub type GL<V> = GeneralLinearGroup<V>;

pub struct GeneralLinearGroup<V: LinearSpace>
where
    DefaultAllocator: Allocator<V::Dim> + Allocator<V::Dim, V::Dim>,
{
    matrix: OMatrix<V::Field, V::Dim, V::Dim>,
}

pub trait DimOfGL {
    type Dim: DimName;
}

impl<N: DimName> DimOfGL for N
where
    N: DimMul<N>,
    DimProd<N, N>: DimName,
{
    type Dim = DimProd<N, N>;
}

impl<V: LinearSpace> Manifold for GeneralLinearGroup<V>
where
    V::Dim: DimOfGL,
    DefaultAllocator: Allocator<V::Dim> + Allocator<V::Dim, V::Dim>,
{
    type Field = V::Field;
    type Dim = <V::Dim as DimOfGL>::Dim;
}

impl<V: LinearSpace> LieGroup for GeneralLinearGroup<V>
where
    V::Dim: DimOfGL,
    DefaultAllocator: Allocator<V::Dim> + Allocator<V::Dim, V::Dim>,
{
    fn identity() -> Self {
        Self {
            matrix: OMatrix::<V::Field, V::Dim, V::Dim>::identity(),
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        Self {
            matrix: &self.matrix * &other.matrix,
        }
    }

    fn inverse(&self) -> Self {
        Self {
            matrix: self
                .matrix
                .clone()
                .try_inverse()
                .expect("Matrix is not invertible."),
        }
    }
}

impl<V: LinearSpace> GroupAction<V> for GeneralLinearGroup<V>
where
    V::Dim: DimOfGL,
    DefaultAllocator: Allocator<V::Dim> + Allocator<V::Dim, V::Dim>,
{
    fn act_on(&self, point: &V) -> V {
        V::_from_raw(&self.matrix * point._get_raw())
    }
}
