use nalgebra as na;

use super::super::space::LinearSpace;
use super::orthogonal::SpecialOrthogonalGroup3D;
use crate::lie::LieGroup;
use crate::manifold::Manifold;

pub type SO2<V> = SpecialEuclideanGroup2D<V>;
pub type SO3<V> = SpecialEuclideanGroup3D<V>;

pub struct SpecialEuclideanGroup2D<V: LinearSpace<2>> {
    translation: V,
    rotation_angle: V::Field,
}

impl<V: LinearSpace<2>> Manifold<2> for SpecialEuclideanGroup2D<V> {
    type Field = V::Field;
}

impl<V: LinearSpace<2>> LieGroup<3> for SpecialEuclideanGroup2D<V> {
    fn identity() -> Self {
        Self {
            translation: V::zero(),
            rotation_angle: V::Field::zero(),
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        let rotation_matrix = na::Rotation2::new(self.angle);
        Self {
            translation: V::_from_raw(
                rotation_matrix * other.translation._get_raw() + self.translation._get_raw(),
            ),
            rotation_angle: self.rotation_angle + other.rotation_angle,
        }
    }

    fn inverse(&self) -> Self {
        let rotation_matrix = na::Rotation2::new(-self.angle);
        Self {
            translation: V::_from_raw(rotation_matrix * -self.translation._get_raw()),
            rotation_angle: -self.rotation_angle,
        }
    }
}

impl<V: LinearSpace<2>> SpecialEuclideanGroup2D<V> {
    pub fn new(translation: V, rotation_angle: V::Field) -> Self {
        Self {
            translation,
            rotation_angle,
        }
    }

    pub fn translation(&self) -> V {
        self.translation
    }

    pub fn rotation_angle(&self) -> V::Field {
        self.rotation_angle
    }
}

pub struct SpecialEuclideanGroup3D<V: LinearSpace<3>> {
    translation: V,
    rotation_quaternion: na::UnitQuaternion<V::Field>,
}

impl<V: LinearSpace<3>> Manifold<3> for SpecialEuclideanGroup3D<V> {
    type Field = V::Field;
}

impl<V: LinearSpace<3>> LieGroup<6> for SpecialEuclideanGroup3D<V> {
    fn identity() -> Self {
        Self {
            translation: V::zero(),
            rotation_quaternion: na::UnitQuaternion::identity(),
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        let rotation_matrix = self.rotation_quaternion.to_rotation_matrix();
        Self {
            translation: V::_from_raw(
                rotation_matrix.matrix() * other.translation._get_raw()
                    + self.translation._get_raw(),
            ),
            rotation_quaternion: self.rotation_quaternion * other.rotation_quaternion,
        }
    }

    fn inverse(&self) -> Self {
        let rotation_matrix = self.rotation_quaternion.to_rotation_matrix().inverse();
        Self {
            translation: V::_from_raw(rotation_matrix.matrix() * -self.translation._get_raw()),
            rotation_quaternion: self.rotation_quaternion.inverse(),
        }
    }
}

impl<V: LinearSpace<3>> SpecialEuclideanGroup3D<V> {
    pub fn new(translation: V, rotation: SpecialOrthogonalGroup3D<V>) -> Self {
        Self {
            translation,
            rotation_quaternion: rotation._quaternion_raw(),
        }
    }

    pub fn translation(&self) -> V {
        self.translation
    }

    pub fn rotation(&self) -> na::UnitQuaternion<V::Field> {
        self.rotation_quaternion
    }
}
