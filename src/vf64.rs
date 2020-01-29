use core::ops::{Add, AddAssign, Mul, MulAssign, SubAssign, Sub, Neg};
use num_traits::{Zero};
use rand::Rng;
use rand_distr::{StandardNormal, Distribution};

#[derive(Copy, Clone, Debug, PartialEq)]
#[repr(C, align(32))]
pub struct Aligned4xf64{
    pub dat: [f64; 4]
}

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
impl MulAssign<f64> for Aligned4xf64{
    #[allow(unreachable_code)]
    fn mul_assign(&mut self, rhs: f64) {
        for a in self.dat.iter_mut(){ *a *= rhs} ;
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
impl Mul<f64> for Aligned4xf64{
    type Output = Self;
    fn mul(self, rhs: f64) -> Self::Output {
        let mut me = self;
        me *= rhs;
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



//#[cfg(test)]
pub mod tests{
    use crate::vf64::Aligned4xf64;
    use crate::vf64::add_assign_4xf64;
    #[allow(dead_code)]
    //#[test]
    pub fn test_aligned_4xf64(){
        let mut a = Aligned4xf64::from([1.0, 2.0, 3.0, 4.0]);
        let b = Aligned4xf64::from([8.0, 7.0, 6.0, 5.0]);
        unsafe{let c1 = add_assign_4xf64(&mut a, &b);}
        let c = a + b;
    }
}