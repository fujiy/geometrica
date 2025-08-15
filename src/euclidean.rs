use nalgebra as na;

pub use crate::linear_space::{
    AffineFrame, AffineSpace, Basis, InnerProductSpace, LinearSpace, OrthogonalAffineFrame,
    OrthogonalBasis, OrthonormalAffineFrame, OrthonormalBasis,
};
pub use crate::manifold::{Chart, LieGroup, Manifold, Scalar, Torsor};

// Vector space associated with Euclidean space

define_vector_and_covector!(Vector, Covector);

impl<const N: usize, K: Scalar> InnerProductSpace<N> for Vector<N, K> {
    fn dot(&self, other: &Self) -> <Self as Manifold<N>>::Field {
        self.raw.dot(&other.raw)
    }
}

// Euclidean space

#[derive(Clone, Debug)]
pub struct EuclideanSpace<const N: usize, K = f64> {
    raw: na::SVector<K, N>,
}

impl<const N: usize, K: Scalar> EuclideanSpace<N, K> {
    pub fn reference_frame() -> OrthonormalAffineFrame<N, Self> {
        OrthonormalAffineFrame::new(
            EuclideanSpace {
                raw: na::SVector::from_element(K::zero()),
            },
            OrthonormalBasis::new(std::array::from_fn(|i| {
                let mut repr = na::SVector::from_element(K::zero());
                repr[i] = K::one();
                Vector { raw: repr }
            })),
        )
    }
}

impl<const N: usize, K: Scalar> Manifold<N> for EuclideanSpace<N, K> {
    type Field = K;
}

impl<const N: usize, K: Scalar> Torsor<N, Vector<N, K>> for EuclideanSpace<N, K> {
    fn act(&self, group_element: &Vector<N, K>) -> Self {
        EuclideanSpace {
            raw: self.raw + group_element.raw,
        }
    }

    fn difference(&self, other: &Self) -> Vector<N, K> {
        Vector {
            raw: self.raw - other.raw,
        }
    }
}

impl<const N: usize, K: Scalar> AffineSpace<N> for EuclideanSpace<N, K> {
    type VectorSpace = Vector<N, K>;
}

pub struct PolarChart<const N: usize, K: Scalar = f64> {
    base: AffineFrame<N, EuclideanSpace<N, K>>,
}

// Translation group T(N)

pub type Translation<const N: usize, K = f64> = Vector<N, K>;

// Rotation group or special orthogonal group SO(N)

pub trait Rotation<const N: usize, K: Scalar = f64>: Manifold<N> + LieGroup<N> {}

pub struct Rotation2D<K: Scalar = f64> {
    raw: K,
}

impl<K: Scalar> Rotation2D<K> {
    pub fn from_angle(angle: K) -> Self {
        Rotation2D { raw: angle }
    }

    pub fn angle(&self) -> K {
        self.raw
    }
}

impl<K: Scalar> Rotation<2, K> for Rotation2D<K> {}

impl<K: Scalar> Manifold<2> for Rotation2D<K> {
    type Field = K;
}

impl<K: Scalar> LieGroup<2> for Rotation2D<K> {
    fn identity() -> Self {
        Rotation2D { raw: K::zero() }
    }

    fn multiply(&self, other: &Self) -> Self {
        Rotation2D {
            raw: self.raw + other.raw,
        }
    }

    fn inverse(&self) -> Self {
        Rotation2D { raw: -self.raw }
    }
}

impl<K: Scalar> Torsor<2, Rotation2D<K>> for OrthonormalBasis<2, Vector<2, K>> {
    fn act(&self, group_element: &Rotation2D<K>) -> Self {
        let rotation_matrix = na::Rotation2::new(group_element.raw);
        OrthonormalBasis::new(std::array::from_fn(|i| Vector {
            raw: rotation_matrix * self[i].raw,
        }))
    }

    fn difference(&self, other: &Self) -> Rotation2D<K> {
        Rotation2D::from_angle(na::Rotation2::rotation_between(&self[0].raw, &other[0].raw).angle())
    }
}
pub struct Rotation3D<K: Scalar = f64> {
    raw: na::UnitQuaternion<K>,
}

impl<K: Scalar> Rotation3D<K> {
    pub fn from_axis_angle(axis: &Vector<3, K>, angle: K) -> Self {
        Rotation3D {
            raw: na::UnitQuaternion::from_axis_angle(&na::Unit::new_normalize(axis.raw), angle),
        }
    }

    pub fn from_quaternion(basis: &Basis<3, Vector<3, K>>, w: K, x: K, y: K, z: K) -> Self {
        let vector = basis.from_local(&[x, y, z]);
        Rotation3D {
            raw: na::UnitQuaternion::from_quaternion(na::Quaternion::new(
                w,
                vector.raw[0],
                vector.raw[1],
                vector.raw[2],
            )),
        }
    }

    pub fn angle(&self) -> K {
        self.raw.angle()
    }

    pub fn axis(&self) -> Option<Vector<3, K>> {
        let raw = self.raw.axis()?;
        Some(Vector { raw: *raw })
    }

    pub fn scaled_axis(&self) -> Vector<3, K> {
        Vector {
            raw: self.raw.scaled_axis(),
        }
    }

    // w, x, y, z
    pub fn quaternion(&self, basis: &Basis<3, Vector<3, K>>) -> [K; 4] {
        let w = self.raw.scalar();
        let vector = Vector {
            raw: self.raw.vector().into(),
        };
        let xyz = basis.to_local(&vector);
        [w, xyz[0], xyz[1], xyz[2]]
    }
}

impl<K: Scalar> Manifold<3> for Rotation3D<K> {
    type Field = K;
}

impl<K: Scalar> LieGroup<3> for Rotation3D<K> {
    fn identity() -> Self {
        Rotation3D {
            raw: na::UnitQuaternion::identity(),
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        Rotation3D {
            raw: self.raw * other.raw,
        }
    }

    fn inverse(&self) -> Self {
        Rotation3D {
            raw: self.raw.inverse(),
        }
    }
}

impl<K: Scalar> Torsor<3, Rotation3D<K>> for OrthonormalBasis<3, Vector<3, K>> {
    fn act(&self, group_element: &Rotation3D<K>) -> Self {
        let rotation_matrix = group_element.raw.to_rotation_matrix();
        OrthonormalBasis::new(std::array::from_fn(|i| Vector {
            raw: rotation_matrix * self[i].raw,
        }))
    }

    fn difference(&self, other: &Self) -> Rotation3D<K> {
        Rotation3D {
            raw: na::UnitQuaternion::rotation_between(&self[0].raw, &other[0].raw)
                .expect("Rotation3D::difference failed"),
        }
    }
}

impl<K: Scalar> Rotation<3, K> for Rotation3D<K> {}

// Motion group or special Euclidean group SE(N)

pub trait Motion<const N: usize, K: Scalar = f64>: From<Translation<N, K>> {
    type RotationType: Rotation<N, K>;

    fn translation(&self) -> Translation<N, K>;
    fn rotation(&self) -> Self::RotationType;

    fn new(
        translation: &Translation<N, K>,
        rotation: &Self::RotationType,
        rotational_origin: &EuclideanSpace<N, K>,
    ) -> Self;

    fn from_rotation(rotation: &Self::RotationType, origin: &EuclideanSpace<N, K>) -> Self {
        Self::new(&Translation::zero(), rotation, origin)
    }
}

pub struct Motion2D<K: Scalar = f64> {
    translation_raw: na::SVector<K, 2>,
    rotation_raw: K,
}

impl<K: Scalar> Manifold<3> for Motion2D<K> {
    type Field = K;
}

impl<K: Scalar> LieGroup<3> for Motion2D<K> {
    fn identity() -> Self {
        Motion2D {
            translation_raw: na::SVector::from_element(K::zero()),
            rotation_raw: K::zero(),
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        let self_rotation_matrix = na::Rotation2::new(self.rotation_raw);
        Motion2D {
            translation_raw: self.translation_raw + self_rotation_matrix * other.translation_raw,
            rotation_raw: self.rotation_raw + other.rotation_raw,
        }
    }

    fn inverse(&self) -> Self {
        let rotation_matrix = na::Rotation2::new(-self.rotation_raw);
        Motion2D {
            translation_raw: -(rotation_matrix * self.translation_raw),
            rotation_raw: -self.rotation_raw,
        }
    }
}

impl<K: Scalar> From<Translation<2, K>> for Motion2D<K> {
    fn from(translation: Translation<2, K>) -> Self {
        Motion2D {
            translation_raw: translation.raw,
            rotation_raw: K::zero(),
        }
    }
}

impl<K: Scalar> Motion<2, K> for Motion2D<K> {
    type RotationType = Rotation2D<K>;

    fn translation(&self) -> Translation<2, K> {
        Translation {
            raw: self.translation_raw,
        }
    }

    fn rotation(&self) -> Self::RotationType {
        Rotation2D {
            raw: self.rotation_raw,
        }
    }

    fn new(
        translation: &Translation<2, K>,
        rotation: &Rotation2D<K>,
        rotational_origin: &EuclideanSpace<2, K>,
    ) -> Self {
        let rotation_matrix = na::Rotation2::new(rotation.raw);
        Motion2D {
            translation_raw: translation.raw - rotation_matrix * rotational_origin.raw
                + rotational_origin.raw,
            rotation_raw: rotation.raw,
        }
    }
}

impl<K: Scalar> Torsor<3, Motion2D<K>> for OrthonormalAffineFrame<2, EuclideanSpace<2, K>> {
    fn act(&self, group_element: &Motion2D<K>) -> Self {
        OrthonormalAffineFrame::new(
            self.origin.act(&group_element.translation()),
            self.basis.act(&group_element.rotation()),
        )
    }

    fn difference(&self, other: &Self) -> Motion2D<K> {
        let translation = self.origin.difference(&other.origin);
        let rotation = Rotation2D::from_angle(
            na::Rotation2::rotation_between(&self.basis[0].raw, &other.basis[0].raw).angle(),
        );
        Motion2D {
            translation_raw: translation.raw,
            rotation_raw: rotation.raw,
        }
    }
}

// 3D Motion group SE(3)

pub struct Motion3D<K: Scalar = f64> {
    translation_raw: na::SVector<K, 3>,
    rotation_raw: na::UnitQuaternion<K>,
}

impl<K: Scalar> Manifold<6> for Motion3D<K> {
    type Field = K;
}

impl<K: Scalar> LieGroup<6> for Motion3D<K> {
    fn identity() -> Self {
        Motion3D {
            translation_raw: na::SVector::from_element(K::zero()),
            rotation_raw: na::UnitQuaternion::identity(),
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        Motion3D {
            translation_raw: self.translation_raw + self.rotation_raw * other.translation_raw,
            rotation_raw: self.rotation_raw * other.rotation_raw,
        }
    }

    fn inverse(&self) -> Self {
        let rotation_raw_inverse = self.rotation_raw.inverse();
        Motion3D {
            translation_raw: -(rotation_raw_inverse * self.translation_raw),
            rotation_raw: rotation_raw_inverse,
        }
    }
}

impl<K: Scalar> From<Translation<3, K>> for Motion3D<K> {
    fn from(translation: Translation<3, K>) -> Self {
        Motion3D {
            translation_raw: translation.raw,
            rotation_raw: na::UnitQuaternion::identity(),
        }
    }
}

impl<K: Scalar> Motion<3, K> for Motion3D<K> {
    type RotationType = Rotation3D<K>;

    fn translation(&self) -> Translation<3, K> {
        Translation {
            raw: self.translation_raw,
        }
    }

    fn rotation(&self) -> Self::RotationType {
        Rotation3D {
            raw: self.rotation_raw,
        }
    }

    fn new(
        translation: &Translation<3, K>,
        rotation: &Rotation3D<K>,
        rotational_origin: &EuclideanSpace<3, K>,
    ) -> Self {
        Motion3D {
            translation_raw: translation.raw - rotation.raw * rotational_origin.raw
                + rotational_origin.raw,
            rotation_raw: rotation.raw,
        }
    }
}

impl<K: Scalar> Torsor<6, Motion3D<K>> for OrthonormalAffineFrame<3, EuclideanSpace<3, K>> {
    fn act(&self, group_element: &Motion3D<K>) -> Self {
        OrthonormalAffineFrame::new(
            self.origin.act(&group_element.translation()),
            self.basis.act(&group_element.rotation()),
        )
    }
    fn difference(&self, other: &Self) -> Motion3D<K> {
        let translation = self.origin.difference(&other.origin);
        let rotation = Rotation3D {
            raw: na::UnitQuaternion::rotation_between(&self.basis[0].raw, &other.basis[0].raw)
                .expect("Motion3D::difference failed"),
        };
        Motion3D {
            translation_raw: translation.raw,
            rotation_raw: rotation.raw,
        }
    }
}
