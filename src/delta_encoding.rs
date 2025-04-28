#![allow(unused)] // TODO: Drop this

use pa_types::{Cost, I};

pub trait VEncoding<Base> {
    fn zero() -> Self;
    fn one() -> Self;
    fn from(p: Base, m: Base) -> Self;
    fn value(&self) -> Cost;
    fn p(&self) -> Base;
    fn m(&self) -> Base;
    fn pm(&self) -> (Base, Base) {
        (self.p(), self.m())
    }
    fn value_of_prefix(&self, j: I) -> Cost;
    fn value_of_suffix(&self, j: I) -> Cost;
    fn value_to(v: &Vec<Self>, j: I) -> Cost
    where
        Self: Sized;
    fn value_from(v: &Vec<Self>, j: I) -> Cost
    where
        Self: Sized;
}

#[derive(Clone, Default, Copy, PartialEq, Eq, Debug)]
pub struct V<Base>(pub Base, pub Base);

macro_rules! impl_vencoding {
    ($($t:ty),+) => {
        $(
            impl VEncoding<$t> for V<$t> {
                #[inline(always)]
                fn zero() -> Self {
                    V(0, 0)
                }
                #[inline(always)]
                fn one() -> Self {
                    V(<$t>::MAX, 0)
                }
                #[inline(always)]
                fn from(p: $t, m: $t) -> Self {
                    V(p, m)
                }
                #[inline(always)]
                fn value(&self) -> Cost {
                    self.0.count_ones() as Cost - self.1.count_ones() as Cost
                }
                /// Value of the first `j` bits.
                /// NOTE: Requires `0 <= j < $t::BITS`.
                #[inline(always)]
                fn value_of_prefix(&self, j: I) -> Cost {
                    debug_assert!(0 <= j && j < <$t>::BITS as I);
                    let mask = (1 << j) - 1;
                    (self.0 & mask).count_ones() as Cost - (self.1 & mask).count_ones() as Cost
                }
                /// Value of the last `j` bits.
                /// NOTE: Requires `j > 0`.
                #[inline(always)]
                fn value_of_suffix(&self, j: I) -> Cost {
                    debug_assert!(0 < j && j <= <$t>::BITS as I);
                    let mask = !(((1 as $t) << (<$t>::BITS as I - j)).wrapping_sub(1));
                    (self.0 & mask).count_ones() as Cost - (self.1 & mask).count_ones() as Cost
                }
                #[inline(always)]
                fn pm(&self) -> ($t, $t) {
                    (self.0, self.1)
                }
                #[inline(always)]
                fn p(&self) -> $t {
                    self.0
                }
                #[inline(always)]
                fn m(&self) -> $t {
                    self.1
                }
                fn value_to(v: &Vec<Self>, j: I) -> Cost {
                    let mut s = 0;
                    for vj in &v[0..j as usize / 64] {
                        s += vj.value();
                    }
                    if j % 64 != 0 {
                        s += v[j as usize / 64].value_of_prefix(j % 64);
                    }
                    s
                }
                fn value_from(v: &Vec<Self>, j: I) -> Cost {
                    let mut s = 0;
                    if j % 64 != 0 {
                        s += v[j as usize / 64].value_of_suffix(64 - j % 64);
                    }
                    for vj in &v[j.div_ceil(64) as usize..] {
                        s += vj.value();
                    }
                    s
                }
            }
        )+
    }
}
impl_vencoding!(u8, u16, u32, u64);

pub trait HEncoding<Base>: Copy {
    fn zero() -> Self;
    fn one() -> Self;
    fn from(p: Base, m: Base) -> Self;
    fn value(&self) -> Cost;
    fn p(&self) -> Base;
    fn m(&self) -> Base;
    #[inline(always)]
    fn pm(&self) -> (Base, Base) {
        (self.p(), self.m())
    }
}

impl HEncoding<u8> for i8 {
    #[inline(always)]
    fn zero() -> Self {
        0
    }
    #[inline(always)]
    fn one() -> Self {
        1
    }
    #[inline(always)]
    fn from(p: u8, m: u8) -> Self {
        p as i8 - m as i8
    }
    #[inline(always)]
    fn value(&self) -> Cost {
        *self as Cost
    }
    #[inline(always)]
    fn p(&self) -> u8 {
        (*self > 0) as u8
    }
    #[inline(always)]
    fn m(&self) -> u8 {
        (*self < 0) as u8
    }
}

// implement HEncoding for all unsigned types.
macro_rules! impl_unsigned {
    ($($t:ty),+) => {
        $(
            impl HEncoding<$t> for ($t, $t) {
                #[inline(always)]
                fn zero() -> Self {
                    (0, 0)
                }
                #[inline(always)]
                fn one() -> Self {
                    (1, 0)
                }
                #[inline(always)]
                fn from(p: $t, m: $t) -> Self {
                    (p as $t, m as $t)
                }
                #[inline(always)]
                fn value(&self) -> Cost {
                    self.0 as Cost - self.1 as Cost
                }
                #[inline(always)]
                fn p(&self) -> $t {
                    self.0 as $t
                }
                #[inline(always)]
                fn m(&self) -> $t {
                    self.1 as $t
                }
            }
        )+
    }
}
impl_unsigned!(u8, u16, u32, u64);
