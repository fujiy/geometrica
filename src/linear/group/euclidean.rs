use nalgebra::{
    DimAdd, DimDiff, DimDiv, DimMul, DimName, DimProd, DimQuot, DimSub, DimSum, Rotation2, U1, U2,
    U3, UnitDualQuaternion,
};

use crate::lie::LieGroup;
use crate::linear::space::{Allocator, DefaultAllocator, InnerProductSpace};
use crate::manifold::*;

pub type SE<V> = SpecialEuclideanGroup<V>;

pub trait SERepr<V: InnerProductSpace>
where
    DefaultAllocator: Allocator<V::Dim>,
{
    type Repr;
    fn identity() -> Self::Repr;
    fn multiply(a: &Self::Repr, b: &Self::Repr) -> Self::Repr;
    fn inverse(a: &Self::Repr) -> Self::Repr;
}

impl<V: InnerProductSpace<Dim = U2>> SERepr<V> for U2
where
    DefaultAllocator: Allocator<V::Dim>,
{
    type Repr = (V, V::Field);
    fn identity() -> Self::Repr {
        (V::zero(), V::Field::zero())
    }
    fn multiply(a: &Self::Repr, b: &Self::Repr) -> Self::Repr {
        let (t1, theta1) = a;
        let (t2, theta2) = b;
        let rotation_matrix = Rotation2::new(*theta1);
        (
            V::_from_raw(rotation_matrix * t2._get_raw() + t1._get_raw()),
            *theta1 + *theta2,
        )
    }
    fn inverse(a: &Self::Repr) -> Self::Repr {
        let (t, theta) = a;
        let rotation_matrix = Rotation2::new(-*theta);
        (V::_from_raw(rotation_matrix * -t._get_raw()), -*theta)
    }
}

impl<V: InnerProductSpace> SERepr<V> for U3
where
    DefaultAllocator: Allocator<V::Dim>,
{
    type Repr = UnitDualQuaternion<V::Field>;
    fn identity() -> Self::Repr {
        UnitDualQuaternion::identity()
    }
    fn multiply(a: &Self::Repr, b: &Self::Repr) -> Self::Repr {
        a * b
    }
    fn inverse(a: &Self::Repr) -> Self::Repr {
        a.inverse()
    }
}

pub trait DimOfSE {
    type Dim: DimName;
}

impl<N: DimName> DimOfSE for N
where
    N: DimSub<U1> + DimMul<DimDiff<N, U1>>,
    DimProd<N, DimDiff<N, U1>>: DimDiv<U2>,
    DimQuot<DimProd<N, DimDiff<N, U1>>, U2>: DimAdd<N>,
    DimSum<DimQuot<DimProd<N, DimDiff<N, U1>>, U2>, N>: DimName,
{
    type Dim = DimSum<DimQuot<DimProd<N, DimDiff<N, U1>>, U2>, N>;
}

pub struct SpecialEuclideanGroup<V: InnerProductSpace>
where
    V::Dim: SERepr<V>,
    DefaultAllocator: Allocator<V::Dim>,
{
    repr: <V::Dim as SERepr<V>>::Repr,
}

impl<V: InnerProductSpace> Manifold for SpecialEuclideanGroup<V>
where
    V::Dim: SERepr<V> + DimOfSE,
    DefaultAllocator: Allocator<V::Dim>,
{
    type Field = V::Field;
    type Dim = <V::Dim as DimOfSE>::Dim;
}

impl<V: InnerProductSpace> LieGroup for SpecialEuclideanGroup<V>
where
    V::Dim: SERepr<V> + DimOfSE,
    DefaultAllocator: Allocator<V::Dim>,
{
    fn identity() -> Self {
        Self {
            repr: <V::Dim as SERepr<V>>::identity(),
        }
    }
    fn multiply(&self, other: &Self) -> Self {
        Self {
            repr: <V::Dim as SERepr<V>>::multiply(&self.repr, &other.repr),
        }
    }
    fn inverse(&self) -> Self {
        Self {
            repr: <V::Dim as SERepr<V>>::inverse(&self.repr),
        }
    }
}
