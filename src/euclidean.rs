use nalgebra as na;

pub use crate::lie_group::{LieGroup, Scalar};

use crate::affine_space::{SpecialEuclideanGroup2D, SpecialEuclideanGroup3D};
use crate::lie_group::Manifold;
use crate::linear_space::{
    InnerProductSpace, LinearSpace, SpecialOrthogonalGroup2D, SpecialOrthogonalGroup3D,
};
use crate::{
    affine_space::AffineSpace,
    lie_group::{LieAlgebra, Torsor},
};

pub struct Vector<const N: usize, K: Scalar = f64> {
    raw: na::SVector<K, N>,
}

impl<const N: usize, K: Scalar> Manifold for Vector<N, K> {
    type Field = K;
}

impl<const N: usize, K: Scalar> LieGroup<N> for Vector<N, K> {
    fn identity() -> Self {
        Self {
            raw: na::SVector::from_element(K::zero()),
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        Self {
            raw: self.raw.clone() + other.raw.clone(),
        }
    }

    fn inverse(&self) -> Self {
        Self {
            raw: -self.raw.clone(),
        }
    }
}

impl<const N: usize, K: Scalar> LinearSpace<N> for Vector<N, K> {
    type DualSpace = Covector<N, K>;

    fn _get_raw(&self) -> &na::SVector<K, N> {
        &self.raw
    }

    fn _from_raw(raw: na::SVector<K, N>) -> Self {
        Self { raw }
    }
}

impl<const N: usize, K: Scalar> InnerProductSpace<N> for Vector<N, K> {
    fn flat(&self) -> Self::DualSpace {
        Covector {
            raw: self.raw.clone(),
        }
    }
    fn sharp(covector: &Self::DualSpace) -> Self {
        Vector {
            raw: covector.raw.clone(),
        }
    }
}

pub struct Covector<const N: usize, K: Scalar = f64> {
    raw: na::SVector<K, N>,
}

impl<const N: usize, K: Scalar> Manifold for Covector<N, K> {
    type Field = K;
}

impl<const N: usize, K: Scalar> LinearSpace<N> for Covector<N, K> {
    type DualSpace = Vector<N, K>;

    fn _get_raw(&self) -> &na::SVector<K, N> {
        &self.raw
    }

    fn _from_raw(raw: na::SVector<K, N>) -> Self {
        Self { raw }
    }
}

impl<const N: usize, K: Scalar> InnerProductSpace<N> for Covector<N, K> {
    fn flat(&self) -> Self::DualSpace {
        Vector {
            raw: self.raw.clone(),
        }
    }
    fn sharp(covector: &Self::DualSpace) -> Self {
        Covector {
            raw: covector.raw.clone(),
        }
    }
}

pub trait Basis<const N: usize> {}
impl<const N: usize, T: Basis<N>> crate::linear_space::Basis<N> for T {
    type Space = Vector<N>;
}

pub trait OrthonormalBasis<const N: usize> {}
impl<const N: usize, T: OrthonormalBasis<N>> crate::linear_space::OrthonormalBasis<N> for T {
    type Space = Vector<N>;
}

// // Euclidean space

pub trait EuclideanSpace<const N: usize> {}

pub trait Frame<const N: usize> {}
impl<const N: usize, T: Frame<N>> crate::affine_space::Frame<N> for T {
    type VectorSpace = Vector<N>;
}

pub trait OrthonormalFrame<const N: usize> {}
impl<const N: usize, T: OrthonormalFrame<N>> crate::affine_space::OrthonormalFrame<N> for T {
    type VectorSpace = Vector<N>;
}

impl<const N: usize, T: EuclideanSpace<N>> AffineSpace<N> for T {
    type V = Vector<N>;
}

// // Translation group T(N)

pub type Translation<const N: usize, K = f64> = Vector<N, K>;

// Rotation group or special orthogonal group SO(N)

pub type Rotation2D<K = f64> = SpecialOrthogonalGroup2D<Vector<2, K>>;
pub type Rotation3D<K = f64> = SpecialOrthogonalGroup3D<Vector<3, K>>;

// impl<K: Scalar> Rotation3D<K> {
// pub fn from_axis_angle(axis: &Vector<3, K>, angle: K) -> Self {
//     Rotation3D {
//         raw: na::UnitQuaternion::from_axis_angle(&na::Unit::new_normalize(axis.raw), angle),
//     }
// }

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

// pub fn angle(&self) -> K {
//     self.raw.angle()
// }

// pub fn axis(&self) -> Option<Vector<3, K>> {
//     let raw = self.raw.axis()?;
//     Some(Vector { raw: *raw })
// }

// pub fn scaled_axis(&self) -> Vector<3, K> {
//     Vector {
//         raw: self.raw.scaled_axis(),
//     }
// }

// // w, x, y, z
// pub fn quaternion(&self, basis: &Basis<3, Vector<3, K>>) -> [K; 4] {
//     let w = self.raw.scalar();
//     let vector = Vector {
//         raw: self.raw.vector().into(),
//     };
//     let xyz = basis.to_local(&vector);
//     [w, xyz[0], xyz[1], xyz[2]]
// }
// }

// // Motion group or special Euclidean group SE(N)

pub type Motion2D<K = f64> = SpecialEuclideanGroup2D<Vector<2, K>>;
pub type Motion3D<K = f64> = SpecialEuclideanGroup3D<Vector<3, K>>;

// pub trait Motion<const N: usize, K: Scalar = f64>: From<Translation<N, K>> {
//     type RotationType: Rotation<N, K>;

//     fn translation(&self) -> Translation<N, K>;
//     fn rotation(&self) -> Self::RotationType;

//     fn new(
//         translation: &Translation<N, K>,
//         rotation: &Self::RotationType,
//         rotational_origin: &EuclideanSpace<N, K>,
//     ) -> Self;

//     fn from_rotation(rotation: &Self::RotationType, origin: &EuclideanSpace<N, K>) -> Self {
//         Self::new(&Translation::zero(), rotation, origin)
//     }
// }

// impl<K: Scalar> From<Translation<2, K>> for Motion2D<K> {
//     fn from(translation: Translation<2, K>) -> Self {
//         Motion2D {
//             translation_raw: translation.raw,
//             rotation_raw: K::zero(),
//         }
//     }
// }

// impl<K: Scalar> Motion<2, K> for Motion2D<K> {
//     type RotationType = Rotation2D<K>;

//     fn translation(&self) -> Translation<2, K> {
//         Translation {
//             raw: self.translation_raw,
//         }
//     }

//     fn rotation(&self) -> Self::RotationType {
//         Rotation2D {
//             raw: self.rotation_raw,
//         }
//     }

//     fn new(
//         translation: &Translation<2, K>,
//         rotation: &Rotation2D<K>,
//         rotational_origin: &EuclideanSpace<2, K>,
//     ) -> Self {
//         let rotation_matrix = na::Rotation2::new(rotation.raw);
//         Motion2D {
//             translation_raw: translation.raw - rotation_matrix * rotational_origin.raw
//                 + rotational_origin.raw,
//             rotation_raw: rotation.raw,
//         }
//     }
// }

// // 3D Motion group SE(3)

// impl<K: Scalar> From<Translation<3, K>> for Motion3D<K> {
//     fn from(translation: Translation<3, K>) -> Self {
//         Motion3D {
//             translation_raw: translation.raw,
//             rotation_raw: na::UnitQuaternion::identity(),
//         }
//     }
// }

// impl<K: Scalar> Motion<3, K> for Motion3D<K> {
//     type RotationType = Rotation3D<K>;

//     fn translation(&self) -> Translation<3, K> {
//         Translation {
//             raw: self.translation_raw,
//         }
//     }

//     fn rotation(&self) -> Self::RotationType {
//         Rotation3D {
//             raw: self.rotation_raw,
//         }
//     }

//     fn new(
//         translation: &Translation<3, K>,
//         rotation: &Rotation3D<K>,
//         rotational_origin: &EuclideanSpace<3, K>,
//     ) -> Self {
//         Motion3D {
//             translation_raw: translation.raw - rotation.raw * rotational_origin.raw
//                 + rotational_origin.raw,
//             rotation_raw: rotation.raw,
//         }
//     }
// }
