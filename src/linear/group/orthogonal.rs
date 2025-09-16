use nalgebra as na;

use crate::lie::LieGroup;
use crate::linear::space::LinearSpace;
use crate::manifold::Manifold;

pub type SO2<V> = SpecialOrthogonalGroup2D<V>;
pub type SO3<V> = SpecialOrthogonalGroup3D<V>;

pub struct SpecialOrthogonalGroup2D<V: LinearSpace<2>> {
    angle: V::Field,
}

impl<V: LinearSpace<2>> Manifold<1> for SpecialOrthogonalGroup2D<V> {
    type Field = V::Field;
}

impl<V: LinearSpace<2>> LieGroup<1> for SpecialOrthogonalGroup2D<V> {
    fn identity() -> Self {
        Self {
            angle: V::Field::zero(),
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        Self {
            angle: self.angle + other.angle,
        }
    }

    fn inverse(&self) -> Self {
        Self { angle: -self.angle }
    }
}

impl<V: LinearSpace<2>> SpecialOrthogonalGroup2D<V> {
    pub fn from_angle(angle: V::Field) -> Self {
        Self { angle }
    }

    pub fn angle(&self) -> V::Field {
        self.angle
    }
}

pub struct SpecialOrthogonalGroup3D<V: LinearSpace<3>> {
    quaternion: na::UnitQuaternion<V::Field>,
}

impl<V: LinearSpace<3>> Manifold<3> for SpecialOrthogonalGroup3D<V> {
    type Field = V::Field;
}

impl<V: LinearSpace<3>> LieGroup<3> for SpecialOrthogonalGroup3D<V> {
    fn identity() -> Self {
        SpecialOrthogonalGroup3D {
            quaternion: na::UnitQuaternion::identity(),
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        SpecialOrthogonalGroup3D {
            quaternion: self.quaternion * other.quaternion,
        }
    }

    fn inverse(&self) -> Self {
        SpecialOrthogonalGroup3D {
            quaternion: self.quaternion.inverse(),
        }
    }
}

impl<V: LinearSpace<3>> SpecialOrthogonalGroup3D<V> {
    pub fn from_axis_angle(axis: &V, angle: V::Field) -> Self {
        SpecialOrthogonalGroup3D {
            quaternion: na::UnitQuaternion::from_axis_angle(
                &na::Unit::new_normalize(axis.raw),
                angle,
            ),
        }
    }

    // pub fn from_quaternion(basis: &Basis<3, Vector<3, K>>, w: K, x: K, y: K, z: K) -> Self {
    //     let vector = basis.from_local(&[x, y, z]);
    //     Rotation3D {
    //         raw: na::UnitQuaternion::from_quaternion(na::Quaternion::new(
    //             w,
    //             vector.raw[0],
    //             vector.raw[1],
    //             vector.raw[2],
    //         )),
    //     }
    // }

    pub fn angle(&self) -> V::Field {
        self.quaternion.angle()
    }

    pub fn axis(&self) -> Option<V> {
        let raw = self.quaternion.axis()?;
        Some(V::_from_raw(*raw))
    }

    pub fn scaled_axis(&self) -> V {
        V::_from_raw(self.raw.scaled_axis())
    }

    // // w, x, y, z
    // pub fn quaternion(&self, basis: &Basis<3, Vector<3, K>>) -> [K; 4] {
    //     let w = self.raw.scalar();
    //     let vector = Vector {
    //         raw: self.raw.vector().into(),
    //     };
    //     let xyz = basis.to_local(&vector);
    //     [w, xyz[0], xyz[1], xyz[2]]
    // }

    pub fn _quaternion_raw(&self) -> &na::UnitQuaternion<V::Field> {
        &self.quaternion
    }
}
