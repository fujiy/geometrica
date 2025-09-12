use nalgebra as na;

use crate::lie_group::{Chart, LieAlgebra, LieAlgebraDual, LieGroup, Manifold, Torsor, Zero};
use crate::linear_space::{GeneralLinearGroup, InnerProductSpace, LinearSpace, OrthonormalBasis};

// Affine Space

pub trait AffineSpace<const N: usize> {
    type V: LinearSpace<N>;
}

impl<const N: usize, V: InnerProductSpace<N> + LieGroup<N>, T: AffineSpace<N, V = V>> Torsor<N, V>
    for T
{
}

// Affine Frame for affine space

pub trait Frame<const N: usize> {
    type VectorSpace: LinearSpace<N>;
}

type Aff<const N: usize, V> = GeneralAffineGroup<N, V>;

pub struct GeneralAffineGroup<const N: usize, V: LinearSpace<N>> {
    // translation: <LinearSpace<N> as LieGroup<N>>::Field,
    translation: V,
    linear_transformation: GeneralLinearGroup<N, V>,
}

macro_rules! impl_general_affine_group {
    ($N:literal, $D:literal) => {
        impl<V: LinearSpace<$N>> Manifold for GeneralAffineGroup<$N, V> {
            type Field = V::Field;
        }

        impl<V: LinearSpace<$N> + LieGroup<$N>> LieGroup<$D> for GeneralAffineGroup<$N, V> {
            fn identity() -> Self {
                Self {
                    translation: V::zero(),
                    linear_transformation: GeneralLinearGroup::identity(),
                }
            }

            fn multiply(&self, other: &Self) -> Self {
                Self {
                    translation: self
                        .linear_transformation
                        .apply(&other.translation)
                        .add(&self.translation),
                    linear_transformation: self
                        .linear_transformation
                        .multiply(&other.linear_transformation),
                }
            }

            fn inverse(&self) -> Self {
                let lin_inv = self.linear_transformation.inverse();
                Self {
                    translation: lin_inv.apply(&self.translation.inverse()),
                    linear_transformation: lin_inv,
                }
            }
        }

        impl<V: LinearSpace<$N> + LieGroup<$N>, T: Frame<$N>> Torsor<$D, GeneralAffineGroup<$N, V>>
            for T
        {
        }

        impl<V: LinearSpace<$N> + LieGroup<$N>> Chart<$N, V> for GeneralAffineGroup<$N, V> {
            fn to_coords(&self, of: &V) -> [V::Field; $N] {
                let translated = of.difference(&self.translation);
                self.linear_transformation.to_coords(&translated)
            }

            fn element_from(&self, coords: &[V::Field; $N]) -> V {
                let lin_part = self.linear_transformation.element_from(coords);
                self.translation.add(&lin_part)
            }
        }
    };
}

impl_general_affine_group!(1, 2);
impl_general_affine_group!(2, 6);
impl_general_affine_group!(3, 12);

// Special Euclidean Group SE(N)

pub trait OrthonormalFrame<const N: usize> {
    type VectorSpace: LinearSpace<N>;
}

type SE2<V> = SpecialEuclideanGroup2D<V>;
type SE3<V> = SpecialEuclideanGroup3D<V>;

pub struct SpecialEuclideanGroup2D<V: LinearSpace<2>> {
    angle: V::Field,
    translation: na::SVector<V::Field, 2>,
}

impl<V: LinearSpace<2>> Manifold for SpecialEuclideanGroup2D<V> {
    type Field = V::Field;
}

impl<V: LinearSpace<2>> LieGroup<3> for SpecialEuclideanGroup2D<V> {
    fn identity() -> Self {
        Self {
            angle: V::Field::zero(),
            translation: na::SVector::<V::Field, 2>::zeros(),
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        let rotation_matrix = na::Rotation2::new(self.angle);
        Self {
            angle: self.angle + other.angle,
            translation: rotation_matrix * other.translation + self.translation,
        }
    }

    fn inverse(&self) -> Self {
        let rotation_matrix = na::Rotation2::new(-self.angle);
        Self {
            angle: -self.angle,
            translation: rotation_matrix * -self.translation,
        }
    }
}

impl<V: LinearSpace<2>, T: OrthonormalFrame<2>> Torsor<3, SpecialEuclideanGroup2D<V>> for T {}

impl<V: LinearSpace<2> + LieGroup<2>> Chart<2, V> for SpecialEuclideanGroup2D<V> {
    fn element_from(&self, coords: &[V::Field; 2]) -> V {
        let rotation_matrix = na::Rotation2::new(self.angle);
        V::_from_raw(
            rotation_matrix * na::SVector::<V::Field, 2>::from_column_slice(coords)
                + self.translation,
        )
    }
    fn to_coords(&self, of: &V) -> [V::Field; 2] {
        let rotation_matrix = na::Rotation2::new(-self.angle);
        (rotation_matrix * (of._get_raw() - self.translation)).into()
    }
}

pub struct SpecialEuclideanGroup3D<V: LinearSpace<3>> {
    quaternion: na::UnitQuaternion<V::Field>,
    translation: na::SVector<V::Field, 3>,
}

impl<V: LinearSpace<3>> Manifold for SpecialEuclideanGroup3D<V> {
    type Field = V::Field;
}

impl<V: LinearSpace<3>> LieGroup<6> for SpecialEuclideanGroup3D<V> {
    fn identity() -> Self {
        Self {
            quaternion: na::UnitQuaternion::identity(),
            translation: na::SVector::<V::Field, 3>::zeros(),
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        Self {
            quaternion: self.quaternion * other.quaternion,
            translation: self.quaternion * other.translation + self.translation,
        }
    }

    fn inverse(&self) -> Self {
        let quat_inv = self.quaternion.inverse();
        Self {
            quaternion: quat_inv,
            translation: quat_inv * -self.translation,
        }
    }
}

impl<V: LinearSpace<3>, T: OrthonormalFrame<3>> Torsor<6, SpecialEuclideanGroup3D<V>> for T {}

impl<V: LinearSpace<3> + LieGroup<3>> Chart<3, V> for SpecialEuclideanGroup3D<V> {
    fn element_from(&self, coords: &[V::Field; 3]) -> V {
        V::_from_raw(
            self.quaternion * na::SVector::<V::Field, 3>::from_column_slice(coords)
                + self.translation,
        )
    }
    fn to_coords(&self, of: &V) -> [V::Field; 3] {
        (self.quaternion.inverse() * (of._get_raw() - self.translation)).into()
    }
}
