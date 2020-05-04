use crate::aligned::{Packet, P4, A32};

use core::ops::{Add, AddAssign, Mul, MulAssign, SubAssign, Sub, Neg};
use num_traits::{Zero, One, FromPrimitive};
use num_traits::Num;
use rand::Rng;
use rand_distr::{StandardNormal, Distribution};
use std::ops::{Div, DivAssign};
use alga::general::{Ring, AbstractGroup, Additive, TwoSidedInverse, AbstractMagma, AbstractMonoid,
                    AbstractSemigroup, AbstractQuasigroup, AbstractLoop, Identity, AbstractRing,
                    Multiplicative, AbstractRingCommutative, AbstractGroupAbelian};
use std::fmt;
use std::fmt::Formatter;

#[derive(Copy, Clone, Debug, PartialEq)]
#[repr(C, align(32))]
pub struct Aligned4xf64{
    pub dat: [f64; 4]
}

impl Aligned4xf64{
    pub fn sum_reduce(&self) -> f64 {
        self.dat.iter().sum()
    }
    pub fn mean_reduce(&self) -> f64 {
        self.sum_reduce() / 4.0
    }
    pub fn sum_sq_reduce(&self) -> f64{
        self.dat.iter().map(|&x| x*x).sum()
    }
    pub fn mean_sq_reduce(&self) -> f64{
        self.sum_sq_reduce() / 4.0
    }

}

impl fmt::Display for Aligned4xf64{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "[ {}, {}, {}, {} ]", self.dat[0], self.dat[1], self.dat[2], self.dat[3] )
    }
}

//impl<'a> Packet<f64, P4, A32> for Aligned4xf64{
//    type ArrayT = [f64; 4];
//    type SliceT = &'a [f64; 4];
//    type SliceMutT = &'a mut [f64; 4];
//
//    fn into_data(self) -> Self::ArrayT {
//        self.dat
//    }
//}

#[derive(Copy, Clone, Debug, PartialEq)]
#[repr(C, align(64))]
pub struct Aligned8xf64{
    pub dat: [f64; 8]
}

#[cfg(target_feature="avx")]
unsafe fn add_assign_4xf64_f64(a: &mut Aligned4xf64, b: f64){
    use core::arch::x86_64;
    let mma = x86_64::_mm256_load_pd(a.dat.as_ptr());
    let mmb = x86_64::_mm256_broadcast_sd(&b);
    let mm_sum = x86_64::_mm256_add_pd(mma, mmb);
    x86_64::_mm256_store_pd(a.dat.as_mut_ptr(), mm_sum);
}

#[cfg(target_feature="avx")]
unsafe fn add_assign_4xf64(a: &mut Aligned4xf64, b: &Aligned4xf64){
    use core::arch::x86_64;
    let mma = x86_64::_mm256_load_pd(a.dat.as_ptr());
    let mmb = x86_64::_mm256_load_pd(b.dat.as_ptr());
    let mm_add = x86_64::_mm256_add_pd(mma, mmb);
    x86_64::_mm256_store_pd(a.dat.as_mut_ptr(), mm_add);
    return;
}

#[cfg(target_feature="avx")]
unsafe fn sub_assign_4xf64(a: &mut Aligned4xf64, b: &Aligned4xf64){
    use core::arch::x86_64;
    let mma = x86_64::_mm256_load_pd(a.dat.as_ptr());
    let mmb = x86_64::_mm256_load_pd(b.dat.as_ptr());
    let mm_sub = x86_64::_mm256_sub_pd(mma, mmb);
    x86_64::_mm256_store_pd(a.dat.as_mut_ptr(), mm_sub);
    return;
}

#[cfg(target_feature="+avx")]
unsafe fn mul_assign_4xf64(a: &mut Aligned4xf64, b: &Aligned4xf64){
    use core::arch::x86_64;
    let mma = x86_64::_mm256_load_pd(a.dat.as_ptr());
    let mmb = x86_64::_mm256_load_pd(b.dat.as_ptr());
    let mm_mul = x86_64::_mm256_mul_pd(mma, mmb);
    x86_64::_mm256_store_pd(a.dat.as_mut_ptr(), mm_mul);
    return;
}

impl Aligned4xf64{
    pub fn map<F>(&self, f: F) -> Aligned4xf64
        where F: Fn(f64) -> f64{
        let mut c = Aligned4xf64::default();
        for (c, &a) in c.dat.iter_mut().zip(self.dat.iter()){
            *c = f(a);
        }
        c
    }
}

impl Distribution<Aligned4xf64> for StandardNormal{
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Aligned4xf64 {
        let mut sample = Aligned4xf64::default();
        for s in sample.dat.iter_mut(){
            *s = StandardNormal.sample(rng);
        }
        sample
    }
}

impl From<[f64; 4]> for Aligned4xf64{
    fn from(dat: [f64; 4]) -> Self {
        Self{dat}
    }
}

impl From<f64> for Aligned4xf64{
    fn from(real: f64) -> Self {
        Self{dat: [real, real, real ,real]}
    }
}

impl Default for Aligned4xf64{
    fn default() -> Self {
        Self::zero()
    }
}
impl Zero for Aligned4xf64{
    fn zero() -> Self {
        Self{dat: [0.0; 4]}
    }

    fn is_zero(&self) -> bool {
        self.dat.iter().all(|a| a.is_zero())
    }
}

impl One for Aligned4xf64{
    fn one() -> Self {
        Self{dat: [1.0; 4]}
    }

    fn is_one(&self ) -> bool{
        self.dat.iter().all(|a| a.is_one())
    }
}

impl Add for Aligned4xf64{
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let mut me = self;
        me += rhs;
        me
    }
}
impl AddAssign<f64> for Aligned4xf64{
    #[allow(unreachable_code)]
    fn add_assign(&mut self, rhs: f64) {
//        #[cfg(target_arch="x86_64")]
//        #[target_feature(enable="avx2")]
//        unsafe {
//            add_assign_4xf64_f64(me, rhs);
//            return;
//        }
        for a in self.dat.iter_mut(){
            *a += rhs;
        }
    }
}
impl Add<f64> for Aligned4xf64{
    type Output = Self;

    fn add(self, rhs: f64) -> Self {
        let mut me = self;
        me += rhs;
        me
    }
}
impl SubAssign for Aligned4xf64{
    #[allow(unreachable_code)]
    fn sub_assign(&mut self, rhs: Self) {
        for (a, &b) in self.dat.iter_mut().zip(rhs.dat.iter()){
            *a -= b;
        }
    }
}
impl Sub for Aligned4xf64{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut me = self;
        me -= rhs;
        me
    }
}
impl AddAssign for Aligned4xf64{
    #[allow(unreachable_code)]
    fn add_assign(&mut self, rhs: Self) {
//        let arr_a = &mut self.dat[0..4];
//        let arr_b = & rhs.dat[0..4];
        for (a, &b) in self.dat.iter_mut().zip(rhs.dat.iter()){
            *a += b
        }
    }
}

impl MulAssign<Aligned4xf64> for Aligned4xf64{
    #[allow(unreachable_code)]
    fn mul_assign(&mut self, rhs: Aligned4xf64) {
        for (a, &b) in self.dat.iter_mut().zip(rhs.dat.iter()){
            *a *= b;
        }
    }
}
impl Mul<Aligned4xf64> for Aligned4xf64{
    type Output = Self;

    fn mul(self, rhs: Aligned4xf64) -> Self::Output {
        let mut me = self;
        me *= rhs;
        me
    }
}
impl MulAssign<f64> for Aligned4xf64{
    #[allow(unreachable_code)]
    fn mul_assign(&mut self, rhs: f64) {
        for a in self.dat.iter_mut(){ *a *= rhs} ;
    }
}
impl Mul<f64> for Aligned4xf64{
    type Output = Self;
    fn mul(self, rhs: f64) -> Self::Output {
        let mut me = self;
        me *= rhs;
        me
    }
}
impl DivAssign<Aligned4xf64> for Aligned4xf64{
    #[allow(unreachable_code)]
    fn div_assign(&mut self, rhs: Aligned4xf64) {
        for (a, &b) in self.dat.iter_mut().zip(rhs.dat.iter()){
            *a /= b;
        }
    }
}
impl Div<Aligned4xf64> for Aligned4xf64{
    type Output = Self;

    fn div(self, rhs: Aligned4xf64) -> Self::Output {
        let mut me = self;
        me /= rhs;
        me
    }
}
impl DivAssign<f64> for Aligned4xf64{
    #[allow(unreachable_code)]
    fn div_assign(&mut self, rhs: f64) {
        for a in self.dat.iter_mut(){ *a /= rhs} ;
    }
}
impl Div<f64> for Aligned4xf64{
    type Output = Self;
    fn div(self, rhs: f64) -> Self::Output {
        let mut me = self;
        me /= rhs;
        me
    }
}

impl Neg for Aligned4xf64{
    type Output = Self;

    fn neg(self) -> Self::Output {
        let mut me = self;
//        #[cfg(target_arch="x86_64")]
//        #[target_feature(enable="avx2")]
//        unsafe {
//            use std::arch::x86_64;
//            let mma = x86_64::_mm256_load_pd(me.dat.as_ptr());
//            let mmb =x86_64::_mm256_broadcast_sd(&(-1.0));
//            let mm_neg = x86_64::_mm256_mul_pd(mma, mmb);
//            x86_64::_mm256_store_pd(me.dat.as_mut_ptr(), mm_neg);
//            return me;
//        }

        for a in me.dat.iter_mut(){
            *a *= -1.0;
        }
        me
    }
}

impl FromPrimitive for Aligned4xf64{
    fn from_i64(n: i64) -> Option<Self> {
        let n  = f64::from_i64(n);
        n.map(|f| Aligned4xf64::from(f))
    }

    fn from_u64(n: u64) -> Option<Self> {
        let n  = f64::from_u64(n);
        n.map(|f| Aligned4xf64::from(f))
    }

    fn from_f64(n: f64) -> Option<Self> {
        Some(Aligned4xf64::from(n))
    }
}


impl Identity<Additive> for Aligned4xf64{
    fn identity() -> Self {
        Self::zero()
    }
}
impl Identity<Multiplicative> for Aligned4xf64{
    fn identity() -> Self {
        Self::one()
    }
}
impl AbstractMagma<Additive> for Aligned4xf64{
    fn operate(&self, right: &Self) -> Self {
        *self + *right
    }
}
impl AbstractMagma<Multiplicative> for Aligned4xf64{
    fn operate(&self, right: &Self) -> Self {
        *self * *right
    }
}
impl TwoSidedInverse<Additive> for Aligned4xf64{
    fn two_sided_inverse(&self) -> Self {
        -*self
    }
}
impl AbstractQuasigroup<Additive> for Aligned4xf64{ }
impl AbstractSemigroup<Additive> for Aligned4xf64{ } impl AbstractSemigroup<Multiplicative> for Aligned4xf64{ }
impl AbstractLoop<Additive> for Aligned4xf64{ }
impl AbstractMonoid<Additive> for Aligned4xf64{ } impl AbstractMonoid<Multiplicative> for Aligned4xf64{ }
impl AbstractGroup<Additive> for Aligned4xf64{}
impl AbstractGroupAbelian<Additive> for Aligned4xf64{}

impl AbstractRing for Aligned4xf64{ }
impl AbstractRingCommutative for Aligned4xf64{ }


#[cfg(test)]
pub mod tests{
    use crate::vf64::Aligned4xf64;

    #[test]
    pub fn test_aligned_4xf64(){
        let mut a = Aligned4xf64::from([1.0, 2.0, 3.0, 4.0]);
        let b = Aligned4xf64::from([8.0, 7.0, 6.0, 5.0]);
        //unsafe{let c1 = add_assign_4xf64(&mut a, &b);}
        let c = a + b;
    }
}