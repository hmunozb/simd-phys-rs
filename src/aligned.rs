pub use num_traits::Num;

pub enum EnumAlignment{
    A8, A16, A32, A64
}
pub fn alignment_of(a: EnumAlignment) -> u32{
    match a{
        EnumAlignment::A8 => 8,
        EnumAlignment::A16 => 16,
        EnumAlignment::A32 => 32,
        EnumAlignment::A64 => 64
    }
}

pub trait Alignment {}
pub struct A8; impl Alignment for A8{}
pub struct A16; impl Alignment for A16{}
pub struct A32; impl Alignment for A32{}
pub struct A64; impl Alignment for A64{}

pub trait AlignedPacket<T: Copy + Clone>{
    type ArrayT;
    fn alignment() -> u32;
    fn data_arr(&self) -> &Self::ArrayT;
    fn data_arr_mut(&mut self) -> &mut Self::ArrayT;
}

pub trait StaticOpDispatch<T: Copy + Clone>{

}

pub trait NumericAlignedPacket<T: Copy + Clone>: AlignedPacket<T> + Num { }

#[derive(Copy, Clone, Debug, PartialEq)]
#[repr(C, align(32))]
pub struct Aligned4x64<T> where T: Copy+Clone{
    pub dat: [T ; 4]
}

impl<T: Copy + Clone> AlignedPacket<T> for Aligned4x32<T>{
    type ArrayT = [T; 4];
    fn alignment() -> u32{
        alignment_of(EnumAlignment::A32)
    }
    fn data_arr(&self) -> &Self::ArrayT{
        &self.dat
    }
    fn data_arr_mut(&mut self) -> &mut Self::ArrayT{
        &mut self.dat
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
#[repr(C, align(64))]
pub struct Aligned8x64<T> where T: Copy+Clone{
    pub dat: [T ; 8]
}

#[derive(Copy, Clone, Debug, PartialEq)]
#[repr(C, align(16))]
pub struct Aligned4x32<T> where T: Copy+Clone{
    pub dat: [T ; 4]
}

#[derive(Copy, Clone, Debug, PartialEq)]
#[repr(C, align(32))]
pub struct Aligned8x32<T> where T: Copy+Clone{
    pub dat: [T ; 8]
}

#[derive(Copy, Clone, Debug, PartialEq)]
#[repr(C, align(64))]
pub struct Aligned16x64<T> where T: Copy+Clone{
    pub dat: [T ; 16]
}