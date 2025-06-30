use nalgebra as na;

pub use crate::manifold::{Scalar, Chart, LieGroup, GroupAction, Manifold};

pub use crate::linear_space::{AffineSpace, AffineFrame, LinearSpace, Basis};

// Vector space associated with Euclidean space

#[derive(Clone)]
pub struct Vector<const N: usize, K: Scalar = f64> {
    raw: na::SVector<K, N>,
}

impl<const N: usize, K: Scalar> Vector<N, K> {

}

impl<const N: usize, K: Scalar> Manifold<N> for Vector<N, K> {
    type Field = K;
}

impl<const N: usize, K: Scalar> AffineSpace<N> for Vector<N, K> {
    type VectorSpace = Self;

    fn translate(&self, other: &Self) -> Self {
        Vector { raw: self.raw + other.raw }
    }

    fn difference(&self, other: &Self) -> Self::VectorSpace {
        Vector { raw: self.raw - other.raw }
    }

    fn reference_frame() -> AffineFrame<N, Self> {
        AffineFrame {
            origin: Vector { raw: na::SVector::from_element(K::zero()) },
            basis: Self::VectorSpace::reference_basis(),
        }
    }
}

impl<const N: usize, K: Scalar> LinearSpace<N> for Vector<N, K> {
    type DualSpace = Covector<N, K>;

    fn zero() -> Self {
        Vector { raw: na::SVector::from_element(K::zero()) }
    }

    fn scale(&self, scalar: K) -> Self {
        Vector { raw: self.raw * scalar }
    }

    fn pair_with(&self, dual: &Self::DualSpace) -> <Self as Manifold<N>>::Field {
        self.raw.dot(&dual.raw)
    }

    fn dual_basis(basis: &[Self; N]) -> [Self::DualSpace; N] {
        let raws: [na::SVector<K, N>; N] = std::array::from_fn(|i| basis[i].raw.clone());
        let matrix = na::SMatrix::<K, N, N>::from_columns(&raws);
        let inverse = matrix.try_inverse().expect("Basis is singular");

        let inverse_t = inverse.transpose();
        std::array::from_fn(|i| Covector { raw: inverse_t.column(i).into_owned() })
    }

    fn reference_basis() -> Basis<N, Self> {
        let vectors = std::array::from_fn(|i| {
            let mut repr = na::SVector::from_element(K::zero());
            repr[i] = K::one();
            Vector { raw: repr }
        });
        Basis {
            dual_basis: Self::dual_basis(&vectors),
            basis: vectors,
        }
    }
}

impl<'a, const N: usize, K: Scalar> std::ops::Add<&'a Vector<N, K>> for &'a Vector<N, K> {
    type Output = Vector<N, K>;

    fn add(self, rhs: Self) -> Self::Output {
        Vector { raw: self.raw + rhs.raw }
    }
}
impl<'a, const N: usize, K: Scalar> std::ops::Sub<&'a Vector<N, K>> for &'a Vector<N, K> {
    type Output = Vector<N, K>;

    fn sub(self, rhs: Self) -> Self::Output {
        Vector { raw: self.raw - rhs.raw }
    }
}
impl<'a, const N: usize, K: Scalar> std::ops::Neg for &'a Vector<N, K> {
    type Output = Vector<N, K>;

    fn neg(self) -> Self::Output {
        Vector { raw: -self.raw }
    }
}
impl<'a, const N: usize, K: Scalar> std::ops::Mul<K> for &'a Vector<N, K> {
    type Output = Vector<N, K>;

    fn mul(self, scalar: K) -> Self::Output {
        Vector { raw: self.raw * scalar }
    }
}

impl<const N: usize, K: Scalar> LieGroup<N> for Vector<N, K> {
    type LieAlgebra = Self;

    fn identity() -> Self {
        Self::zero()
    }

    fn multiply(&self, other: &Self) -> Self {
        self + other
    }

    fn inverse(&self) -> Self {
        - self
    }
}

impl<const N: usize, K: Scalar> GroupAction<N, EuclideanSpace<N, K>> for Vector<N, K> {
    fn act_on(&self, point: &EuclideanSpace<N, K>) -> EuclideanSpace<N, K> {
        point.translate(self)
    }
}


// Covector space associated with Euclidean space

#[derive(Clone)]
pub struct Covector<const N: usize, K: Scalar = f64> {
    raw: na::SVector<K, N>,
}

impl<const N: usize, K: Scalar> Manifold<N> for Covector<N, K> {
    type Field = K;
}

impl<const N: usize, K: Scalar> AffineSpace<N> for Covector<N, K> {
    type VectorSpace = Self;

    fn translate(&self, other: &Self) -> Self {
        Covector { raw: self.raw + other.raw }
    }
    fn difference(&self, other: &Self) -> Self::VectorSpace {
        Covector { raw: self.raw - other.raw }
    }

    fn reference_frame() -> AffineFrame<N, Self> {
        AffineFrame {
            origin: Covector { raw: na::SVector::from_element(K::zero()) },
            basis: Self::VectorSpace::reference_basis(),
        }
    }
}

impl<const N: usize, K: Scalar> LinearSpace<N> for Covector<N, K> {
    type DualSpace = Vector<N, K>;


    fn zero() -> Self {
        Covector { raw: na::SVector::from_element(K::zero()) }
    }

    fn scale(&self, scalar: K) -> Self {
        Covector { raw: self.raw * scalar }
    }

    fn pair_with(&self, dual: &Self::DualSpace) -> <Self as Manifold<N>>::Field {
        self.raw.dot(&dual.raw)
    }

    fn dual_basis(basis: &[Self; N]) -> [Self::DualSpace; N] {
        let raws: [na::SVector<K, N>; N] = std::array::from_fn(|i| basis[i].raw.clone());
        let matrix = na::SMatrix::<K, N, N>::from_columns(&raws);
        let inverse = matrix.try_inverse().expect("Basis is singular");

        let inverse_t = inverse.transpose();
        std::array::from_fn(|i| Vector { raw: inverse_t.column(i).into_owned() })
    }

    fn reference_basis() -> Basis<N, Self> {
        let vectors = std::array::from_fn(|i| {
            let mut repr = na::SVector::from_element(K::zero());
            repr[i] = K::one();
            Covector { raw: repr }
        });
        Basis {
            dual_basis: Self::dual_basis(&vectors),
            basis: vectors,
        }
    }
}

impl<const N: usize, K: Scalar> std::ops::Mul<K> for Covector<N, K> {
    type Output = Self;

    fn mul(self, scalar: K) -> Self::Output {
        Covector { raw: self.raw * scalar }
    }
}


// Euclidean space 


// #[derive(Clone)]
pub struct EuclideanSpace<const N: usize, K = f64> {
    raw: na::SVector<K, N>,
}

impl<const N: usize, K> EuclideanSpace<N, K> {

}

impl<const N: usize, K: Scalar> Manifold<N> for EuclideanSpace<N, K> {
    type Field = K;

}

impl<const N: usize, K: Scalar> AffineSpace<N> for EuclideanSpace<N, K> {
    type VectorSpace = Vector<N, K>;

    fn translate(&self, other: &Self::VectorSpace) -> Self {
        EuclideanSpace { raw: self.raw + other.raw }
    }

    fn difference(&self, other: &Self) -> Self::VectorSpace {
        Vector { raw: self.raw - other.raw }
    }

    fn reference_frame() -> AffineFrame<N, Self> {
        AffineFrame {
            origin: EuclideanSpace { raw: na::SVector::from_element(K::zero()) },
            basis: Self::VectorSpace::reference_basis(),
        }
    }
}


pub struct PolarChart<const N: usize, K: Scalar = f64> {
    base: AffineFrame<N, EuclideanSpace<N, K>>,
}


// Translation group T(N)

type Translation<const N: usize, K: Scalar = f64> = Vector<N, K>;

// Rotation group SO(N)

pub struct Rotation2D<K: Scalar = f64> {
    angle: K,
}

impl<K: Scalar> Rotation2D<K> {
    pub fn from_angle(angle: K) -> Self {
        Rotation2D { angle }
    }

    pub fn to_angle(&self) -> K {
        self.angle
    }
}

impl<K: Scalar> Manifold<2> for Rotation2D<K> {
    type Field = K;
}

impl<K: Scalar> LieGroup<2> for Rotation2D<K> {
    type LieAlgebra = Vector<2, K>;

    fn identity() -> Self {
        Rotation2D { angle: K::zero() }
    }

    fn multiply(&self, other: &Self) -> Self {
        Rotation2D { angle: self.angle + other.angle }
    }

    fn inverse(&self) -> Self {
        Rotation2D { angle: -self.angle }
    }
}

impl<K: Scalar> GroupAction<2, Vector<2, K>> for Rotation2D<K> {
    fn act_on(&self, point: &Vector<2, K>) -> Vector<2, K> {
        let rotation_matrix = na::Rotation2::new(self.angle);
        Vector {
            raw: rotation_matrix * point.raw
        }
    }
}


struct Rotation3D<K: Scalar = f64> {
    quaternion: na::UnitQuaternion<K>,
}

impl<K: Scalar> Manifold<3> for Rotation3D<K> {
    type Field = K;
}

impl<K: Scalar> LieGroup<3> for Rotation3D<K> {
    type LieAlgebra = Vector<3, K>;

    fn identity() -> Self {
        Rotation3D { quaternion: na::UnitQuaternion::identity() }
    }

    fn multiply(&self, other: &Self) -> Self {
        Rotation3D { quaternion: self.quaternion * other.quaternion }
    }

    fn inverse(&self) -> Self {
        Rotation3D { quaternion: self.quaternion.inverse() }
    }
}

impl<K: Scalar> GroupAction<3, Vector<3, K>> for Rotation3D<K> {
    fn act_on(&self, point: &Vector<3, K>) -> Vector<3, K> {
        Vector { raw: self.quaternion * point.raw }
    }
}


// Motion group SE(N)

pub struct Motion2D<K: Scalar = f64> {
    translation: Translation<2, K>,
    rotation: Rotation2D<K>,
}

impl<K: Scalar> Manifold<2> for Motion2D<K> {
    type Field = K;
}

impl<K: Scalar> LieGroup<2> for Motion2D<K> {
    type LieAlgebra = Vector<2, K>;

    fn identity() -> Self {
        Motion2D { 
            translation: Translation::zero(), 
            rotation: Rotation2D::identity() 
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        Motion2D { 
            translation: &self.translation + &other.translation, 
            rotation: self.rotation.multiply(&other.rotation) 
        }
    }

    fn inverse(&self) -> Self {
        Motion2D { 
            translation: - &self.translation, 
            rotation: self.rotation.inverse() 
        }
    }
    
}


pub struct Motion3D<K: Scalar = f64> {
    translation: Translation<3, K>,
    rotation: Rotation3D<K>,
}

impl<K: Scalar> Manifold<3> for Motion3D<K> {
    type Field = K;
}

impl<K: Scalar> LieGroup<3> for Motion3D<K> {
    type LieAlgebra = Vector<3, K>;

    fn identity() -> Self {
        Motion3D { 
            translation: Translation::zero(), 
            rotation: Rotation3D::identity() 
        }
    }

    fn multiply(&self, other: &Self) -> Self {
        Motion3D { 
            translation: &self.translation + &other.translation, 
            rotation: self.rotation.multiply(&other.rotation) 
        }
    }

    fn inverse(&self) -> Self {
        Motion3D { 
            translation: - &self.translation, 
            rotation: self.rotation.inverse() 
        }
    }
}


