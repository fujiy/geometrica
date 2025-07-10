use core::marker::PhantomData;

/// 不変参照 (`&'id mut ()`) でライフタイム `'id` を *不変* に保持
// #[derive(Copy, Clone)]
// struct InvariantLifetime<'id>(PhantomData<fn(&'id mut ())>);

#[derive(Clone, Copy, Default)]
struct InvariantLifetime<'id>(PhantomData<*mut &'id ()>);

struct CoordinateSystem<const N: usize, T = f64> {
    basis: [Vector<N, T>; N],
}

impl<const N: usize, T> CoordinateSystem<N, T> where T: Copy + Default {
    fn vector(&self, repr: [T; N]) -> Vector<N, T> {
        Vector { raw: repr }
    }

    fn reference() -> CoordinateSystem<N, T> {
        CoordinateSystem::<N, T> {
            basis: std::array::from_fn(|i| {
                let mut repr = [T::default(); N];
                repr[i] = T::default();
                Vector { raw: repr }
            }),
        }
    }
}

// #[derive(Copy, Clone)]
struct Vector<const N: usize, T = f64> {
    raw: [T; N],
}

impl<const N: usize, T> Vector<N, T> {
    // fn get_raw(&self) -> &'id i32 {
    //     &self.value
    // }


    fn in_local<F, R>(&self, coord: CoordinateSystem<N, T>, f: F) -> R where
        T: Copy,
        F: for<'id> FnOnce(Local<'id, N, T>) -> R {
        let local = Local { 
            value: self.raw,
            _marker: InvariantLifetime(PhantomData),
        };
        f(local)
    }

    // fn in_local_ref<F, R>(&self, coord: &CoordinateSystem<N, T>) -> Local<'id, N, T> where {
    //     Local { 
    //         value: self.raw,
    //         _marker: InvariantLifetime(PhantomData),
    //     }
    // }

    // fn in_local_2<F, R>(&self, coord: CoordinateSystem<N, T>, f: F) -> R where
    //     T: Copy,
    //     F: for<'id> FnOnce(&'id Local<'id, N, T>, &'id Local<'id, N, T>) -> R {
    //     let local1 = Local { 
    //         value: self.raw,
    //         _marker: InvariantLifetime(PhantomData),
    //     };
    //     let local2 = Local { 
    //         value: self.raw,
    //         _marker: InvariantLifetime(PhantomData),
    //     };
    //     f(&local1, &local2)
    // }
}

impl<const N: usize, T> std::ops::Add for &Vector<N, T> where
    T: std::ops::Add<Output = T> + Copy,{
    type Output = Vector<N, T>;
    fn add(self, rhs: Self) -> Self::Output {
        Vector { 
            raw: std::array::from_fn(|i| self.raw[i] + rhs.raw[i]) 
        }
    }
}

struct Local<'id, const N: usize, T> {
    value: [T; N],
    _marker: InvariantLifetime<'id>,
}

impl<'id, const N: usize, T> std::ops::Add for Local<'id, N, T> where
    T: std::ops::Add<Output = T> + Copy,{
    type Output = Local<'id, N, T>;
    fn add(self, rhs: Self) -> Self::Output {
        Local { 
            value: std::array::from_fn(|i| self.value[i] + rhs.value[i]),
            _marker: self._marker,
        }
    }
}

// impl<'id> std::ops::Add<&B<'id>> for &B<'id> {
//     type Output = B<'id>;
//     fn add(self, rhs: &B<'id>) -> B<'id> {
//         B {
//             value: self.value + rhs.value,
//             _brand: self._brand,
//         }
//     }
// }

/// 呼び出しごとに *全く新しい* `'brand` を生成

#[test]
fn test() {
    let ref_coord = CoordinateSystem::<3>::reference();
    let other_coord = CoordinateSystem::<3>::reference();

    let v1 = ref_coord.vector([1.0, 2.0, 3.0]);
    let v2 = ref_coord.vector([4.0, 5.0, 6.0]);
    let v3 = other_coord.vector([4.0, 5.0, 6.0]);

    // v1.in_local(ref_coord, |l1| {
    //     v2.in_local(ref_coord, |l2| {
    //         let _ = l1 + l2; // ✔ 同じ座標系なので OK
    //     });
    // });

    // with_a(|a1| {
    //     let b1 = a1.make_b(3);
    //     let b2 = a1.make_b(4);
    //     let b3 = a1.make_b(5);
    //     let b4 = a1.make_b(5);
    //     let b5 = a1.make_b(5);
    //     let _ = b1 + b2;           // ✔ 同じブランドなので OK
    //     // let _ = &b1 + &b2;           // ✔ 同じブランドなので OK

    //     with_a(|a2| {
    //         let c1 = a2.make_b(1);
    //         let c2 = a2.make_b(2);
    //         let c3 = a2.make_b(3);
    //         let _ = c1 + c2;    // ✔ 同じブランドなので OK
    //         // let _ = c3 + b3;    // ❌ コンパイルエラー:
    //                                //     型が `B<'brand1>` と `B<'brand2>`
    //     });

    //     // let _ = b4 + b5;           // ✔ 同じブランドなので OK
    // });
}

// #[test]
// fn test2() {
//     let a1 = A { _brand: InvariantLifetime(PhantomData) };
//     let b1 = a1.make_b(3);
//     let b2 = a1.make_b(4);
//     let _ = b1 + b2;           // ✔ 同じブランドなので OK
//     let _ = b1 + b1;           // ✔ 同じブランドなので OK

//     let a2 = A { _brand: InvariantLifetime(PhantomData) };
//     let c1 = a2.make_b(1);
//     // let _ = b1 + c1;    // ❌ コンパイルエラー:
//     //     型が `B<'brand1>` と `B<'brand2>`
// }