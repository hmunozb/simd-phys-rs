pub use num_traits::Num;
use std::cell::Ref;
use std::borrow::Borrow;
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

pub trait Alignment {fn val()->u32;}
pub struct A8; impl Alignment for A8{fn val() -> u32{8}}
pub struct A16; impl Alignment for A16{fn val() -> u32{16}}
pub struct A32; impl Alignment for A32{fn val() -> u32{32}}
pub struct A64; impl Alignment for A64{fn val() -> u32{64}}

pub trait Packing {fn val()->u32;}
pub struct P2; impl Packing for P2{fn val() ->u32{2}}
pub struct P4; impl Packing for P4{fn val() ->u32{4}}
pub struct P8; impl Packing for P8{fn val() ->u32{8}}
pub struct P16; impl Packing for P16{fn val() ->u32{16}}
pub struct P32; impl Packing for P32{fn val() ->u32{32}}
pub struct P64; impl Packing for P64{fn val() ->u32{64}}

pub trait Packet<T, P: Packing, A: Alignment>
where for<'b> &'b<Self as Packet<T, P, A>>::ArrayT : IntoIterator,
      for<'b> &'b mut <Self as Packet<T, P, A>>::ArrayT : IntoIterator,
      for<'b> &'b<<Self as Packet<T, P, A>>::ArrayT as IntoIterator>::IntoIter : ExactSizeIterator,
      for<'b> &'b mut <<Self as Packet<T, P, A>>::ArrayT as IntoIterator>::IntoIter : ExactSizeIterator
{
    type ArrayT : IntoIterator<Item=T>;
    //type SliceT : IntoIterator<Item=T>;
    //type SliceMutT : IntoIterator<Item=T>;

    fn into_data(self) -> Self::ArrayT;
    fn data_iter(&self) -> &(<Self::ArrayT as IntoIterator>::IntoIter);
    fn data_iter_mut(&mut self) -> &mut (<Self::ArrayT as IntoIterator>::IntoIter);
//    fn data_arr(&self) -> &Self::ArrayT;
//    fn data_arr_mut(&self) -> &mut Self::ArrayT;
//
//    fn into_iter(self) -> Self::IterT;
//    fn data_iter(&self) -> Self::IterT;
}


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