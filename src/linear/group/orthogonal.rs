use crate::lie::{GroupAction, LieGroup};
use crate::linear::space::{Allocator, DefaultAllocator, InnerProductSpace};
use crate::manifold::*;
use nalgebra::{
    DimDiff, DimDiv, DimMul, DimName, DimProd, DimQuot, DimSub, Rotation2, U1, U2, U3, Unit,
    UnitQuaternion,
};

pub type SO<V> = SpecialOrthogonalGroup<V>;

pub trait SORepr<V: InnerProductSpace>
where
    DefaultAllocator: Allocator<<V as Manifold>::Dim>,
{
    type Repr;

    fn identity() -> Self::Repr;
    fn multiply(a: &Self::Repr, b: &Self::Repr) -> Self::Repr;
    fn inverse(a: &Self::Repr) -> Self::Repr;

    fn act(a: &Self::Repr, point: &V) -> V;
}

impl<V: InnerProductSpace<Dim = U2>> SORepr<V> for U2
where
    DefaultAllocator: Allocator<<V as Manifold>::Dim>,
{
    type Repr = V::Field;
    fn identity() -> Self::Repr {
        V::Field::zero()
    }
    fn multiply(a: &Self::Repr, b: &Self::Repr) -> Self::Repr {
        *a + *b
    }
    fn inverse(a: &Self::Repr) -> Self::Repr {
        a.neg()
    }
    fn act(a: &Self::Repr, point: &V) -> V {
        V::_from_raw(Rotation2::new(*a) * point._get_raw())
    }
}

impl<V: InnerProductSpace<Dim = U3>> SORepr<V> for U3
where
    DefaultAllocator: Allocator<<V as Manifold>::Dim>,
{
    type Repr = UnitQuaternion<V::Field>;

    fn identity() -> Self::Repr {
        UnitQuaternion::identity()
    }
    fn multiply(a: &Self::Repr, b: &Self::Repr) -> Self::Repr {
        a * b
    }
    fn inverse(a: &Self::Repr) -> Self::Repr {
        a.inverse()
    }
    fn act(a: &Self::Repr, point: &V) -> V {
        V::_from_raw(a * point._get_raw())
    }
}

pub struct SpecialOrthogonalGroup<V: InnerProductSpace>
where
    V::Dim: SORepr<V>,
    DefaultAllocator: Allocator<V::Dim>,
{
    repr: <V::Dim as SORepr<V>>::Repr,
}

// pub type DimOfSO<N: DimName> = DimQuot<DimProd<N, DimDiff<N, U1>>, U2>;

pub trait DimOfSO {
    type Dim: DimName;
}

impl<N: DimName> DimOfSO for N
where
    N: DimSub<U1> + DimMul<DimDiff<N, U1>>,
    DimProd<N, DimDiff<N, U1>>: DimDiv<U2>,
    DimQuot<DimProd<N, DimDiff<N, U1>>, U2>: DimName,
{
    type Dim = DimQuot<DimProd<N, DimDiff<N, U1>>, U2>;
}

impl<V: InnerProductSpace> Manifold for SpecialOrthogonalGroup<V>
where
    V::Dim: SORepr<V> + DimOfSO,
    DefaultAllocator: Allocator<V::Dim>,
{
    type Field = V::Field;
    type Dim = <V::Dim as DimOfSO>::Dim;
}

impl<V: InnerProductSpace> LieGroup for SpecialOrthogonalGroup<V>
where
    V::Dim: SORepr<V> + DimOfSO,
    DefaultAllocator: Allocator<V::Dim>,
{
    fn identity() -> Self {
        Self {
            repr: <V::Dim as SORepr<V>>::identity(),
        }
    }
    fn multiply(&self, other: &Self) -> Self {
        Self {
            repr: <V::Dim as SORepr<V>>::multiply(&self.repr, &other.repr),
        }
    }
    fn inverse(&self) -> Self {
        Self {
            repr: <V::Dim as SORepr<V>>::inverse(&self.repr),
        }
    }
}

impl<V: InnerProductSpace> GroupAction<V> for SpecialOrthogonalGroup<V>
where
    V::Dim: SORepr<V> + DimOfSO,
    DefaultAllocator: Allocator<V::Dim>,
{
    fn act_on(&self, point: &V) -> V {
        <V::Dim as SORepr<V>>::act(&self.repr, point)
    }
}

pub trait SpecialOrthogonalGroup2D<V: InnerProductSpace>
where
    DefaultAllocator: Allocator<V::Dim>,
{
    fn angle(&self) -> V::Field;
}

impl<V: InnerProductSpace<Dim = U2>> SpecialOrthogonalGroup<V> {
    pub fn from_angle(angle: V::Field) -> Self {
        Self { repr: angle }
    }
}

impl<V: InnerProductSpace<Dim = U2>> SpecialOrthogonalGroup2D<V> for SpecialOrthogonalGroup<V> {
    fn angle(&self) -> V::Field {
        self.repr
    }
}

pub trait SpecialOrthogonalGroup3D<V: InnerProductSpace>
where
    DefaultAllocator: Allocator<V::Dim>,
{
    fn angle(&self) -> V::Field;
}

impl<V: InnerProductSpace<Dim = U3>> SpecialOrthogonalGroup<V> {
    pub fn from_axis_angle(axis: &V, angle: V::Field) -> Self {
        Self {
            repr: UnitQuaternion::from_axis_angle(
                &Unit::new_normalize(axis._get_raw().clone()),
                angle,
            ),
        }
    }

    pub fn axis(&self) -> Option<V> {
        let raw = self.repr.axis()?;
        Some(V::_from_raw(*raw))
    }

    pub fn scaled_axis(&self) -> V {
        V::_from_raw(self.repr.scaled_axis())
    }
}

impl<V: InnerProductSpace<Dim = U3>> SpecialOrthogonalGroup3D<V> for SpecialOrthogonalGroup<V> {
    fn angle(&self) -> V::Field {
        self.repr.angle()
    }
}

// pub struct SpecialOrthogonalGroup2D<V: LinearSpace> {
//     angle: V::Field,
// }

// impl<V: LinearSpace<2>> Manifold<1> for SpecialOrthogonalGroup2D<V> {
//     type Field = V::Field;
// }

// impl<V: LinearSpace<2>> LieGroup<1> for SpecialOrthogonalGroup2D<V> {
//     fn identity() -> Self {
//         Self {
//             angle: V::Field::zero(),
//         }
//     }

//     fn multiply(&self, other: &Self) -> Self {
//         Self {
//             angle: self.angle + other.angle,
//         }
//     }

//     fn inverse(&self) -> Self {
//         Self { angle: -self.angle }
//     }
// }

// impl<V: LinearSpace<2>> SpecialOrthogonalGroup2D<V> {
//     pub fn from_angle(angle: V::Field) -> Self {
//         Self { angle }
//     }

//     pub fn angle(&self) -> V::Field {
//         self.angle
//     }
// }

// pub struct SpecialOrthogonalGroup3D<V: LinearSpace<3>> {
//     quaternion: na::UnitQuaternion<V::Field>,
// }

// impl<V: LinearSpace<3>> Manifold<3> for SpecialOrthogonalGroup3D<V> {
//     type Field = V::Field;
// }

// impl<V: LinearSpace<3>> LieGroup<3> for SpecialOrthogonalGroup3D<V> {
//     fn identity() -> Self {
//         SpecialOrthogonalGroup3D {
//             quaternion: na::UnitQuaternion::identity(),
//         }
//     }

//     fn multiply(&self, other: &Self) -> Self {
//         SpecialOrthogonalGroup3D {
//             quaternion: self.quaternion * other.quaternion,
//         }
//     }

//     fn inverse(&self) -> Self {
//         SpecialOrthogonalGroup3D {
//             quaternion: self.quaternion.inverse(),
//         }
//     }
// }

// impl<V: LinearSpace<3>> SpecialOrthogonalGroup3D<V> {
//     pub fn from_axis_angle(axis: &V, angle: V::Field) -> Self {
//         SpecialOrthogonalGroup3D {
//             quaternion: na::UnitQuaternion::from_axis_angle(
//                 &na::Unit::new_normalize(axis.raw),
//                 angle,
//             ),
//         }
//     }

//     // pub fn from_quaternion(basis: &Basis<3, Vector<3, K>>, w: K, x: K, y: K, z: K) -> Self {
//     //     let vector = basis.from_local(&[x, y, z]);
//     //     Rotation3D {
//     //         raw: na::UnitQuaternion::from_quaternion(na::Quaternion::new(
//     //             w,
//     //             vector.raw[0],
//     //             vector.raw[1],
//     //             vector.raw[2],
//     //         )),
//     //     }
//     // }

//     pub fn angle(&self) -> V::Field {
//         self.quaternion.angle()
//     }

//     pub fn axis(&self) -> Option<V> {
//         let raw = self.quaternion.axis()?;
//         Some(V::_from_raw(*raw))
//     }

//     pub fn scaled_axis(&self) -> V {
//         V::_from_raw(self.raw.scaled_axis())
//     }

//     // // w, x, y, z
//     // pub fn quaternion(&self, basis: &Basis<3, Vector<3, K>>) -> [K; 4] {
//     //     let w = self.raw.scalar();
//     //     let vector = Vector {
//     //         raw: self.raw.vector().into(),
//     //     };
//     //     let xyz = basis.to_local(&vector);
//     //     [w, xyz[0], xyz[1], xyz[2]]
//     // }

//     pub fn _quaternion_raw(&self) -> &na::UnitQuaternion<V::Field> {
//         &self.quaternion
//     }
// }
