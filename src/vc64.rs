use core::ops::{Add, AddAssign, Mul, MulAssign, SubAssign, Sub, Neg};
use num_traits::{Zero};
use num_complex::Complex64 as c64;
use crate::vf64::Aligned4xf64;
use std::fmt;

#[inline]
#[cfg(target_feature="avx")]
unsafe fn mul_2xc64(a0: c64, a1: c64, b0: c64, b1: c64) -> (c64, c64){
    use core::arch::x86_64::*;
    let mut c_arr = Aligned4xf64::zero();

    let a = _mm256_set_pd(a0.re, a0.im, a1.re, a1.im);
    let b = _mm256_set_pd(b0.re, b0.im, b1.re, b1.im);
    // [a1 a0 a3 a2]
    let ap = _mm256_permute_pd(a, 0x05);

    // [a0b0  a1b1  a2b2  a3b3]
    let d0 = _mm256_mul_pd(a, b);
    // [a1b0  a0b1  a3b2  a2b3]
    let d1 = _mm256_mul_pd(ap, b);

    // [a0b0  a1b0  a2b2  a3b2]
    let el = _mm256_unpacklo_pd(d0, d1);
    // [a1b1  a0b3  a3b3  a2b3]
    let eh = _mm256_unpackhi_pd(d0, d1);
    //let e0 = _mm256_shuffle_pd(d0, d1, 0x00);
    //let e1 = _mm256_shuffle_pd(d0, d1, 0x0F);

    let c = _mm256_addsub_pd(el, eh);
    _mm256_store_pd(c_arr.dat.as_mut_ptr(), c);


    return (c64{re: c_arr.dat[0], im: c_arr.dat[1]}, c64{re: c_arr.dat[2], im: c_arr.dat[3]});
}

#[derive(Copy, Clone, Debug, PartialEq)]
#[repr(C, align(32))]
pub struct Aligned4xc64{
    pub dat: [c64; 4]
}

impl fmt::Display for Aligned4xc64{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[ {}, {}, {}, {} ]", self.dat[0], self.dat[1], self.dat[2], self.dat[3] )
    }
}

impl MulAssign<f64> for Aligned4xc64{
    #[allow(unreachable_code)]
    fn mul_assign(&mut self, rhs: f64) {
        for a in self.dat.iter_mut(){ *a *= rhs} ;
    }
}

impl MulAssign<c64> for Aligned4xc64{
    #[allow(unreachable_code)]
    fn mul_assign(&mut self, rhs: c64) {
        for a in self.dat.iter_mut(){ *a *= rhs} ;
    }
}

impl MulAssign for Aligned4xc64{
    #[allow(unreachable_code)]
    fn mul_assign(&mut self, rhs: Self) {
        #[cfg(target_feature="avx")]
        {
            for (a, b) in self.dat.chunks_exact_mut(2)
                .zip(rhs.dat.chunks_exact(2)) {

                let c = unsafe{mul_2xc64(a[0], a[1], b[0], b[1])};
                a[0] = c.0;
                a[1] = c.1;

            } ;
            return;
        }
        //fallback
        for (a, &b) in self.dat.iter_mut()
            .zip(rhs.dat.iter())
        {
            *a *= b;
        }
    }
}

impl Mul for Aligned4xc64{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        let mut me = self;
        me *= rhs;
        return me;
    }
}

pub fn multiply(a_arr: &[Aligned4xc64], b_arr: &[Aligned4xc64], c_arr: &mut [Aligned4xc64]){
    for ((&a, &b), c) in a_arr.iter().zip(b_arr.iter()).zip(c_arr.iter_mut()){
        *c = a * b;
    }
}

pub fn multiply_single(a_arr: &[c64], b_arr: &[c64], c_arr: &mut [c64]){
    let a_ch = a_arr.chunks_exact(4);
    let b_ch = b_arr.chunks_exact(4);
    let c_ch = c_arr.chunks_exact_mut(4);

    for ((a, b), c) in a_ch.clone().zip(b_ch.clone()).zip(c_ch){
        for ((&ai, &bi), ci) in a.iter().zip(b.iter()).zip(c.iter_mut()){
            *ci = ai * bi;
        }
    }
    for((&ai, &bi), ci) in a_ch.remainder().iter().zip(
        b_ch.remainder().iter()).zip(c_arr.chunks_exact_mut(4).into_remainder()){
        *ci = ai * bi;
    }
}