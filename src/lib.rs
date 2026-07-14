// Copyright 2017 10x Genomics

//! # debruijn: a De Bruijn graph library for DNA seqeunces in Rust.
//! This library provides tools for efficient construction DeBruijn graphs (dBG)
//! from DNA sequences, tracking arbitrary metadata associated with kmers in the
//! graph, and performing path-compression of unbranched graph paths to improve
//! speed and reduce memory consumption.
//!
//! Most applications of `debruijn` will follow this general workflow:
//! 1. You generate a set of sequences to make a dBG from.
//! 2. You pass those sequences to the `filter_kmers` function, which converts the sequences into kmers, while tracking 'metadata' about each kmer in a very customizable way. The metadata could be read count, a set of colors, a set of read counts split by haplotype, a UMI count, etc.
//! 3. The the library will convert the kmers to a compressed dBG. You can also customize the rules for how to compress the dBG and how to 'combine' the per-kmer metadata.
//!
//! Then you can use the final compressed dBG how you like. There are some methods for simplifying and re-building the  graph, but those could be developed more.
//!
//! ## Examples
//! - [Local phased SV assembly tool in our Long Ranger package](https://github.com/10XGenomics/longranger/blob/master/lib/pvc/src/asm_caller.rs#L205)
//! - [Single-cell VDJ assember](https://github.com/10XGenomics/cellranger/blob/master/lib/rust/vdj_asm/src/asm.rs#L191)
//! - [Build a colored, compressed dBG of a transcriptome reference](https://github.com/10XGenomics/rust-pseudoaligner/blob/master/src/build_index.rs#L40)
//!
//! All the data structures in debruijn-rs are specialized to the 4 base DNA alphabet,
//! and use 2-bit packed encoding of base-pairs into integer types, and efficient methods for
//! reverse complement, enumerating kmers from longer sequences, and transfering data between
//! sequences.
//!
//! ## Encodings
//! Most methods for ingesting sequence data into the library have a form named 'bytes',
//! which expects bases encoded as the integers 0,1,2,3, and a separate form names 'ascii',
//! which expects bases encoded as the ASCII letters A,C,G,T.

use bimap::BiMap;
use clap::ValueEnum;
use rand::Rng;
use rand::rngs::ThreadRng;
use serde_derive::{Deserialize, Serialize};
use summarizer::Marker;
use std::fmt::{self, Debug, Display};
use std::hash::Hash;
use std::marker::PhantomData;
use std::{array, mem};
use std::ops::Range;

use crate::compression::{CheckCompress, compress_kmers_with_hash};
use crate::dna_string::DnaString;
use crate::filter::filter_kmers;
use crate::reads::{ReadData, Reads, ReadsPaired};
use crate::serde::{SerGraph, SerKmers, SerReads};
use crate::summarizer::{ID, SampleInfo, SummaryConfig, SummaryData, Tag, Translator};

pub mod clean_graph;
pub mod compression;
pub mod dna_string;
pub mod reads;
pub mod filter;
pub mod summarizer;
pub mod graph;
pub mod kmer;
pub mod msp;
pub mod neighbors;
pub mod vmer;
pub mod fastq;
pub mod colors;
pub mod serde;

const BUF: usize = 64*1024;
const BUCKETS: usize = 256;
const ALPHABET_SIZE: usize = 4;
const PROGRESS_STYLE: &str = "{msg} [{elapsed_precise}] {bar:60.cyan/blue} ({pos}/{len})";

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
mod bitops_avx2;

#[cfg(test)]
pub mod test;

/// Convert a 2-bit representation of a base to a char
#[inline]
pub fn bits_to_ascii(c: u8) -> u8 {
    match c {
        0u8 => b'A',
        1u8 => b'C',
        2u8 => b'G',
        3u8 => b'T',
        _ => b'X',
    }
}

/// Convert an ASCII-encoded DNA base to a 2-bit representation,
/// transforming bytes outside of ACGTacgt to A
#[inline]
pub fn base_to_bits(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 0u8,
        b'C' | b'c' => 1u8,
        b'G' | b'g' => 2u8,
        b'T' | b't' => 3u8,
        _ => 0u8,
    }
}

/// Convert an ASCII-encoded DNA base to a 2-bit representation,
/// second value is `false` if the base was ambiguous
#[inline]
pub fn base_to_bits_checked(c: u8) -> (u8, bool) {
    match c {
        b'A' | b'a' => (0u8, true),
        b'C' | b'c' => (1u8, true),
        b'G' | b'g' => (2u8, true),
        b'T' | b't' => (3u8, true),
        _ => (0u8, false)
    }
}

#[inline]
pub fn dna_only_base_to_bits(c: u8) -> Option<u8> {
    match c {
        b'A' | b'a' => Some(0u8),
        b'C' | b'c' => Some(1u8),
        b'G' | b'g' => Some(2u8),
        b'T' | b't' => Some(3u8),
        _ => None,
    }
}

/// Convert an ASCII-encoded DNA base to a 2-bit representation
#[inline]
pub fn is_valid_base(c: u8) -> bool {
    matches!(c, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
}

/// Convert a 2-bit representation of a base to a char
#[inline]
pub fn bits_to_base(c: u8) -> char {
    match c {
        0u8 => 'A',
        1u8 => 'C',
        2u8 => 'G',
        3u8 => 'T',
        _ => 'X',
    }
}

/// The complement of a 2-bit encoded base
#[inline(always)]
pub fn complement(base: u8) -> u8 {
    (!base) & 0x3u8
}

/// Trait for interacting with DNA sequences
pub trait Mer: Sized + fmt::Debug {
    /// Length of DNA sequence
    fn len(&self) -> usize;

    /// True if the sequence is empty.
    fn is_empty(&self) -> bool;

    /// Get 2-bit encoded base at position `pos`
    fn get(&self, pos: usize) -> u8;

    /// Set base at `pos` to 2-bit encoded base `val`
    fn set_mut(&mut self, pos: usize, val: u8);

    /// Set `nbases` positions in the sequence, starting at `pos`.
    /// Values must  be packed into the upper-most bits of `value`.
    fn set_slice_mut(&mut self, pos: usize, nbases: usize, value: u64);

    /// Return a new object containing the reverse complement of the sequence
    fn rc(&self) -> Self;

    /// Iterate over the bases in the sequence
    fn iter(&'_ self) -> MerIter<'_, Self> {
        MerIter {
            sequence: self,
            i: 0,
        }
    }

    /// Count the number of A/T bases in the kmer
    fn at_count(&self) -> u32 {
        let mut count = 0;
        for i in 0..self.len() {
            let base = self.get(i);
            if base == 0 || base == 3 {
                count += 1;
            }
        }
        count
    }

    /// Count the number of G/C bases in the kmer
    fn gc_count(&self) -> u32 {
        let mut count = 0;
        for i in 0..self.len() {
            let base = self.get(i);
            if base == 1 || base == 2 {
                count += 1;
            }
        }
        count
    }
}

/// Iterator over bases of a DNA sequence (bases will be unpacked into bytes).
pub struct MerIter<'a, M: 'a + Mer> {
    sequence: &'a M,
    i: usize,
}

impl<'a, M: 'a + Mer> Iterator for MerIter<'a, M> {
    type Item = u8;

    fn next(&mut self) -> Option<u8> {
        if self.i < self.sequence.len() {
            let value = self.sequence.get(self.i);
            self.i += 1;
            Some(value)
        } else {
            None
        }
    }
}

/// Encapsulates a Kmer sequence with statically known K.
pub trait Kmer: Mer + Sized + Copy + PartialEq + PartialOrd + Eq + Ord + Hash {
    /// Create a Kmer initialized to all A's
    fn empty() -> Self;

    /// K value for this concrete type.
    fn k() -> usize;

    /// Return the rank of this kmer in an lexicographic ordering of all kmers
    /// E.g. 'AAAA' -> 0, 'AAAT' -> 1, etc. This will panic if K > 32.
    fn to_u64(&self) -> u64;

    /// Construct a kmer from the given lexicographic rank of the kmer.
    /// If K > 32, the leads bases will be A's.
    fn from_u64(value: u64) -> Self;

    // Compute Hamming distance between self and other
    fn hamming_dist(&self, other: Self) -> u32;

    /// Add the base `v` to the left side of the sequence, and remove the rightmost base
    fn extend_left(&self, v: u8) -> Self;

    /// Add the base `v` to the right side of the sequence, and remove the leftmost base
    fn extend_right(&self, v: u8) -> Self;

    /// Add the base `v` to the side of sequence given by `dir`, and remove a base at the opposite side
    fn extend(&self, v: u8, dir: Dir) -> Self {
        match dir {
            Dir::Left => self.extend_left(v),
            Dir::Right => self.extend_right(v),
        }
    }

    /// Generate all the extension of this sequence given by `exts` in direction `Dir`
    fn get_extensions(&self, exts: Exts, dir: Dir) -> Vec<Self> {
        let ext_bases = exts.get(dir);
        ext_bases.iter().map(|b| self.extend(*b, dir)).collect()
    }

    /// Return the minimum of the kmer and it's reverse complement, and a flag indicating if sequence was flipped
    fn min_rc_flip(&self) -> (Self, bool) {
        let rc = self.rc();
        //println!("kmer flip: self: {:?}, rc: {:?}, t/f: {}", self, rc, (*self < rc));
        if *self < rc {
            (*self, false)
        } else {
            (rc, true)
        }
    }

    /// Return the minimum of the kmer and it's reverse complement
    fn min_rc(&self) -> Self {
        let rc = self.rc();
        if *self < rc {
            *self
        } else {
            rc
        }
    }

    /// Test if this Kmer and it's reverse complement are the same
    fn is_palindrome(&self) -> bool {
        self.len().is_multiple_of(2) && *self == self.rc()
    }

    /// Create a Kmer from the first K bytes of `bytes`, which must be encoded as the integers 0-4.
    fn from_bytes(bytes: &[u8]) -> Self {
        if bytes.len() < Self::k() {
            panic!("bytes not long enough to form kmer")
        }

        let mut k0 = Self::empty();

        for (i, b) in bytes.iter().take(Self::k()).enumerate() {
            k0.set_mut(i, *b)
        }

        k0
    }

    /// Create a Kmer from the first K bytes of `bytes`, which must be encoded as ASCII letters A,C,G, or T.
    fn from_ascii(bytes: &[u8]) -> Self {
        if bytes.len() < Self::k() {
            panic!("bytes not long enough to form kmer")
        }

        let mut k0 = Self::empty();

        for (i, b) in bytes.iter().take(Self::k()).enumerate() {
            k0.set_mut(i, base_to_bits(*b))
        }

        k0
    }

    /// Return String containing Kmer sequence
    fn to_string(&self) -> String {
        let mut s = String::with_capacity(self.len());
        for pos in 0..self.len() {
            s.push(bits_to_base(self.get(pos)))
        }
        s
    }

    /// Generate vector of all kmers contained in `str` encoded as 0-4.
    fn kmers_from_bytes(str: &[u8]) -> Vec<Self> {
        if str.len() < Self::k() {
            return Vec::default();
        }
        let mut k0 = Self::empty();
        for (i, v) in str.iter().take(Self::k()).enumerate() {
            k0.set_mut(i, *v);
        }

        let mut r = Vec::with_capacity(str.len() - Self::k() + 1);
        r.push(k0);

        for v in str.iter().skip(Self::k()) {
            k0 = k0.extend_right(*v);
            r.push(k0);
        }

        r
    }

    /// Generate vector of all kmers contained in `str`, encoded as ASCII ACGT.
    fn kmers_from_ascii(str: &[u8]) -> Vec<Self> {
        if str.len() < Self::k() {
            return Vec::default();
        }
        let mut k0 = Self::empty();
        for (i, b) in str.iter().take(Self::k()).enumerate() {
            k0.set_mut(i, base_to_bits(*b));
        }

        let mut r = Vec::with_capacity(str.len() - Self::k() + 1);
        r.push(k0);

        for v in str.iter().skip(Self::k()) {
            k0 = k0.extend_right(base_to_bits(*v));
            r.push(k0);
        }

        r
    }

    fn has_low_complexity(&self) -> bool {
        let a = Self::from_u64(0);
        let c = Self::from_u64((0..(Self::k()*2)).filter(|&x| (x % 2 == 0) | (x == 0)).map(|x| 2usize.pow(x as u32)).sum::<usize>() as u64);
        let g = Self::from_u64((0..(Self::k()*2)).filter(|&x| x % 2 != 0).map(|x| 2u64.pow(x as u32)).sum::<u64>());
        let t = Self::from_u64(2u64.pow((Self::k()*2) as u32) - 1);

        (self == &a) | (self == &c) | (self == &g) | (self == &t)
    }
}

/// An immutable interface to a Mer sequence.
pub trait MerImmut: Mer + Clone {
    fn set(&self, pos: usize, val: u8) -> Self {
        let mut new = self.clone();
        new.set_mut(pos, val);
        new
    }

    fn set_slice(&self, pos: usize, nbases: usize, bits: u64) -> Self {
        let mut new = self.clone();
        new.set_slice_mut(pos, nbases, bits);
        new
    }
}

impl<T> MerImmut for T where T: Mer + Clone {}

/// A DNA sequence with run-time variable length, up to a statically known maximum length
pub trait Vmer: Mer + PartialEq + Eq {
    /// Create a new sequence with length `len`, initialized to all A's
    fn new(len: usize) -> Self;

    /// Maximum sequence length that can be stored in this type
    fn max_len() -> usize;

    /// Create a Vmer from a sequence of bytes
    fn from_slice(seq: &[u8]) -> Self {
        let mut vmer = Self::new(seq.len());
        for (i, v) in seq.iter().enumerate() {
            vmer.set_mut(i, *v);
        }

        vmer
    }

    /// Efficiently extract a Kmer from the sequence
    fn get_kmer<K: Kmer>(&self, pos: usize) -> K;

    /// Get the first Kmer from the sequence
    fn first_kmer<K: Kmer>(&self) -> K {
        self.get_kmer(0)
    }

    /// Get the last kmer in the sequence
    fn last_kmer<K: Kmer>(&self) -> K {
        self.get_kmer(self.len() - K::k())
    }

    /// Get the terminal kmer of the sequence, on the both side of the sequence
    fn both_term_kmer<K: Kmer>(&self) -> (K, K) {
        (self.first_kmer(), self.last_kmer())
    }

    /// Get the terminal kmer of the sequence, on the side of the sequence given by dir
    fn term_kmer<K: Kmer>(&self, dir: Dir) -> K {
        match dir {
            Dir::Left => self.first_kmer(),
            Dir::Right => self.last_kmer(),
        }
    }

    /// Iterate over the kmers in the sequence
    fn iter_kmers<K: Kmer>(&self) -> KmerIter<'_, K, Self> {
        let kmer = if self.len() >= K::k() {
            self.first_kmer()
        } else {
            // Default kmer, will not be used
            K::empty()
        };

        KmerIter {
            bases: self,
            kmer,
            pos: K::k(),
        }
    }

    /// Iterate over the kmers and their extensions, given the extensions of the whole sequence
    fn iter_kmer_exts<K: Kmer>(&self, seq_exts: Exts) -> KmerExtsIter<'_, K, Self> {
        let kmer = if self.len() >= K::k() {
            self.first_kmer()
        } else {
            // Default kmer, will not be used
            K::empty()
        };

        KmerExtsIter {
            bases: self,
            exts: seq_exts,
            kmer,
            pos: K::k(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct KmerDataItem<K: Kmer, DI> {
    pub kmer: K,
    pub exts: Exts,
    pub data: DI,
    pub quality: Option<BaseQuality>
}

impl<K: Kmer, DI> KmerDataItem<K, DI> {
    pub fn new(kmer: K, exts: Exts, data: DI, quality: Option<BaseQuality>) -> KmerDataItem<K, DI> {
        KmerDataItem { kmer, exts, data, quality }
    }
}

/// A newtype wrapper around a `Vec<u8>` with implementations
/// of the `Mer` and `Vmer` traits.
#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct DnaBytes(pub Vec<u8>);

impl Mer for DnaBytes {
    fn len(&self) -> usize {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn get(&self, pos: usize) -> u8 {
        self.0[pos]
    }

    /// Set base at `pos` to 2-bit encoded base `val`
    fn set_mut(&mut self, pos: usize, val: u8) {
        self.0[pos] = val
    }

    /// Set `nbases` positions in the sequence, starting at `pos`.
    /// Values must  be packed into the upper-most bits of `value`.
    fn set_slice_mut(&mut self, _pos: usize, _nbases: usize, _value: u64) {
        unimplemented!();
        //for i in pos .. (pos + nbases) {
        //
        //}
    }

    /// Return a new object containing the reverse complement of the sequence
    fn rc(&self) -> Self {
        unimplemented!();
    }
}

impl Vmer for DnaBytes {
    /// Create a new sequence with length `len`, initialized to all A's
    fn new(len: usize) -> Self {
        DnaBytes(vec![0u8; len])
    }

    /// Maximum sequence length that can be stored in this type
    fn max_len() -> usize {
        1 << 48
    }

    /// Efficiently extract a Kmer from the sequence
    fn get_kmer<K: Kmer>(&self, pos: usize) -> K {
        K::from_bytes(&self.0[pos..pos + K::k()])
    }
}

/// A newtype wrapper around a `&[u8]` with implementations
/// of the `Mer` and `Vmer` traits.
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct DnaSlice<'a>(pub &'a [u8]);

impl Mer for DnaSlice<'_> {
    fn len(&self) -> usize {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn get(&self, pos: usize) -> u8 {
        self.0[pos]
    }

    /// Set base at `pos` to 2-bit encoded base `val`
    fn set_mut(&mut self, _pos: usize, _val: u8) {
        unimplemented!()
    }

    /// Set `nbases` positions in the sequence, starting at `pos`.
    /// Values must  be packed into the upper-most bits of `value`.
    fn set_slice_mut(&mut self, _pos: usize, _nbases: usize, _value: u64) {
        unimplemented!();
        //for i in pos .. (pos + nbases) {
        //
        //}
    }

    /// Return a new object containing the reverse complement of the sequence
    fn rc(&self) -> Self {
        unimplemented!();
    }
}

impl Vmer for DnaSlice<'_> {
    /// Create a new sequence with length `len`, initialized to all A's
    fn new(_len: usize) -> Self {
        unimplemented!();
    }

    /// Maximum sequence length that can be stored in this type
    fn max_len() -> usize {
        1 << 48
    }

    /// Efficiently extract a Kmer from the sequence
    fn get_kmer<K: Kmer>(&self, pos: usize) -> K {
        K::from_bytes(&self.0[pos..pos + K::k()])
    }
}

/// Direction of motion in a DeBruijn graph
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq)]
pub enum Dir {
    Left,
    Right,
}

impl Dir {
    /// Return a fresh Dir with the opposite direction
    pub fn flip(&self) -> Dir {
        match *self {
            Dir::Left => Dir::Right,
            Dir::Right => Dir::Left,
        }
    }

    /// Return a fresh Dir opposite direction if do_flip == True
    pub fn cond_flip(&self, do_flip: bool) -> Dir {
        if do_flip {
            self.flip()
        } else {
            *self
        }
    }

    /// Pick between two alternatives, depending on the direction
    pub fn pick<T>(&self, if_left: T, if_right: T) -> T {
        match self {
            Dir::Left => if_left,
            Dir::Right => if_right,
        }
    }

    /// get the index of the base in dir for [`Exts`] and [`EdgeMult`]
    fn index(&self, base: u8) -> u8 {
        match self {
            Self::Right => ALPHABET_SIZE as u8 - 1 - base,
            Self::Left => 2 * ALPHABET_SIZE as u8 - 1 - base,
        }
    }

    /// get the index range of the dir for [`Exts`] and [`EdgeMult`]
    fn index_range(&self) -> Range<usize>{
        match self {
            Self::Right => 0..ALPHABET_SIZE,
            Self::Left => ALPHABET_SIZE..(2 * ALPHABET_SIZE),
        }
    }
}

/// Store single-base extensions for a DNA Debruijn graph.
///
/// 8 bits, 4 higher order ones represent extensions to the right, 4 lower order ones
/// represent extensions to the left. For each direction the bits (from lower order
/// to higher order) represent whether there exists an extension with each of the
/// letters A, C, G, T. So overall the bits are:
///  right   left
/// T G C A T G C A
#[derive(Eq, PartialEq, Copy, Clone, Ord, PartialOrd, Hash, Serialize, Deserialize)]
pub struct Exts {
    pub val: u8,
}

impl Exts {
    pub fn new(val: u8) -> Self {
        Exts { val }
    }

    pub fn empty() -> Exts {
        Exts { val: 0u8 }
    }

    pub fn from_single_dirs(left: Exts, right: Exts) -> Exts {
        Exts {
            val: (right.val << 4) | (left.val & 0xf),
        }
    }

    pub fn merge(left: Exts, right: Exts) -> Exts {
        Exts {
            val: left.val & 0x0f | right.val & 0xf0,
        }
    }

    pub fn add(&self, v: Exts) -> Exts {
        Exts {
            val: self.val | v.val,
        }
    }

    /// subtract an Exts from an Exts
    pub fn subtract(&self, v: Exts) -> Exts {
        Exts { val: self.val & !v.val }
    }

    pub fn set(&self, dir: Dir, pos: u8) -> Exts {
        let shift = pos
            + match dir {
                Dir::Right => 4,
                Dir::Left => 0,
            };

        let new_val = self.val | (1u8 << shift);
        Exts { val: new_val }
    }

    pub fn remove(&self, dir: Dir, pos: u8) -> Exts {
        let shift = pos
            + match dir {
                Dir::Right => 4,
                Dir::Left => 0,
            };

        let new_val = self.val & !(1u8 << shift);
        Exts { val: new_val }
    }

    #[inline]
    fn dir_bits(&self, dir: Dir) -> u8 {
        match dir {
            Dir::Right => self.val >> 4,
            Dir::Left => self.val & 0xf,
        }
    }

    pub fn get(&self, dir: Dir) -> Vec<u8> {
        let bits = self.dir_bits(dir);
        let mut v = Vec::with_capacity(4);
        for i in 0..4 {
            if bits & (1 << i) > 0 {
                v.push(i);
            }
        }

        v
    }

    pub fn has_ext(&self, dir: Dir, base: u8) -> bool {
        let bits = self.dir_bits(dir);
        (bits & (1 << base)) > 0
    }

    pub fn from_slice_bounds(src: &[u8], start: usize, length: usize) -> Exts {
        let l_extend = if start > 0 {
            1u8 << (src[start - 1])
        } else {
            0u8
        };
        let r_extend = if start + length < src.len() {
            1u8 << src[start + length]
        } else {
            0u8
        };

        Exts {
            val: (r_extend << 4) | l_extend,
        }
    }

    pub fn from_dna_string(src: &dna_string::DnaString, start: usize, length: usize) -> Exts {
        let l_extend = if start > 0 {
            1u8 << (src.get(start - 1))
        } else {
            0u8
        };
        let r_extend = if start + length < src.len() {
            1u8 << src.get(start + length)
        } else {
            0u8
        };

        Exts {
            val: (r_extend << 4) | l_extend,
        }
    }

    pub fn num_exts_l(&self) -> u8 {
        self.num_ext_dir(Dir::Left)
    }

    pub fn num_exts_r(&self) -> u8 {
        self.num_ext_dir(Dir::Right)
    }

    pub fn num_ext_dir(&self, dir: Dir) -> u8 {
        let e = self.dir_bits(dir);
        (e & 1u8) + ((e & 2u8) >> 1) + ((e & 4u8) >> 2) + ((e & 8u8) >> 3)
    }

    pub fn mk_left(base: u8) -> Exts {
        Exts::empty().set(Dir::Left, base)
    }

    pub fn mk_right(base: u8) -> Exts {
        Exts::empty().set(Dir::Right, base)
    }

    pub fn mk(left_base: u8, right_base: u8) -> Exts {
        Exts::merge(Exts::mk_left(left_base), Exts::mk_right(right_base))
    }

    pub fn get_unique_extension(&self, dir: Dir) -> Option<u8> {
        if self.num_ext_dir(dir) != 1 {
            None
        } else {
            let e = self.dir_bits(dir);
            for i in 0..4 {
                if (e & (1 << i)) > 0 {
                    return Some(i);
                }
            }

            None
        }
    }

    pub fn single_dir(&self, dir: Dir) -> Exts {
        match dir {
            Dir::Right => Exts { val: self.val >> 4 },
            Dir::Left => Exts {
                val: self.val & 0xfu8,
            },
        }
    }

    /// Complement the extension bases for each direction
    pub fn complement(&self) -> Exts {
        let v = self.val;

        // swap bits
        let mut r = (v & 0x55u8) << 1 | ((v >> 1) & 0x55u8);

        // swap pairs
        r = (r & 0x33u8) << 2 | ((r >> 2) & 0x33u8);
        Exts { val: r }
    }

    pub fn reverse(&self) -> Exts {
        let v = self.val;
        let r = (v & 0xf) << 4 | (v >> 4);
        Exts { val: r }
    }

    pub fn rc(&self) -> Exts {
        self.reverse().complement()
    }
}

impl fmt::Debug for Exts {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();

        for b in self.get(Dir::Left) {
            s.push(bits_to_base(b));
        }
        s.push('|');

        for b in self.get(Dir::Right) {
            s.push(bits_to_base(b));
        }

        write!(f, "{}", s)
    }
}

/// Iterate over the `Kmer`s of a DNA sequence efficiently
pub struct KmerIter<'a, K: Kmer, D>
where
    D: 'a,
{
    bases: &'a D,
    kmer: K,
    pos: usize,
}

impl<K: Kmer, D: Mer> Iterator for KmerIter<'_, K, D> {
    type Item = K;

    #[inline]
    fn next(&mut self) -> Option<K> {
        if self.pos <= self.bases.len() {
            let retval = self.kmer;

            if self.pos < self.bases.len() {
                self.kmer = self.kmer.extend_right(self.bases.get(self.pos));
            }

            self.pos += 1;
            Some(retval)
        } else {
            None
        }
    }
}

/// Iterate over the `(Kmer, Exts)` tuples of a sequence and it's extensions efficiently
pub struct KmerExtsIter<'a, K: Kmer, D>
where
    D: 'a,
{
    bases: &'a D,
    exts: Exts,
    kmer: K,
    pos: usize,
}

impl<K: Kmer, D: Mer> Iterator for KmerExtsIter<'_, K, D> {
    type Item = (K, Exts);

    fn next(&mut self) -> Option<(K, Exts)> {
        if self.pos <= self.bases.len() {
            let next_base = if self.pos < self.bases.len() {
                self.bases.get(self.pos)
            } else {
                0u8
            };

            let cur_left = if self.pos == K::k() {
                self.exts
            } else {
                Exts::mk_left(self.bases.get(self.pos - K::k() - 1))
            };

            let cur_right = if self.pos < self.bases.len() {
                Exts::mk_right(next_base)
            } else {
                self.exts
            };

            let cur_exts = Exts::merge(cur_left, cur_right);

            let retval = self.kmer;
            self.kmer = self.kmer.extend_right(next_base);
            self.pos += 1;
            Some((retval, cur_exts))
        } else {
            None
        }
    }
}


/// Compress up to 64 tags to one `u64` (8 bytes) (or up to 128 tags to a  
/// `u128`(16 bytes) with the feature `sample128` enabled)
#[derive(Clone, PartialEq, Copy, Serialize, Deserialize)]
pub struct Tags {
    pub val: Marker,
}

impl Tags {

    /// Make a new Tags from a `u64` (or `u128` with the feature `sample128` enabled)
    pub fn new(val: Marker) -> Self {
        Tags { val }
    }

    /// get number of labels saved in the Tags
    pub fn len(&self) -> usize {
        self.val.count_ones() as usize
    }

    /// check if the tags are empty
    pub fn is_empty(&self) -> bool {
        self.val == 0
    }

    /// encodes a sorted (!) Vec<Tag> and encodes it as a u64
    pub fn from_tag_vec(vec: &Vec<Tag>) -> Self {
        let mut x = 0;
        
        // if the vector is empty, return an empty Tags
        if vec.is_empty() { return Tags { val: 0 } }

        // panic if Tags would overflow
        if ( *vec.last().expect("vector empty when it shouldn't be") ) / 8 as Tag >= mem::size_of::<Tags>() as Tag { 
            panic!("too many tags - maximum number of supported tags is 64 by default, 128 with compile flag / feature '--feature sample128'") 
        }
        
        // iterate backwards over all elements of the vector
        for i in (1..vec.len()).rev() {
            x += 1;
            x <<= vec[i] - vec[i-1];
        }

        x += 1;
        x <<= vec[0];

        Tags { val: x }
    }

    // turn Tags into Vec<Tag>
    pub fn to_tag_vec(&self) -> Vec<Tag> {
        let mut x = self.val;
        let mut vec: Vec<Tag> = Vec::new();

        // do bit-wise right shifts trough u64
        // each time first digit is 1 (is an odd number), push i to vec
        for i in 0..(mem::size_of::<Tags>()*8) as Tag {
            if !x.is_multiple_of(2) {
                vec.push(i)
            }
            x >>= 1;
        }

        vec
    }

    // directly translate Tags to Vec<&str>
    // str_map is translatror BiMap between Tag and &str 
    pub fn to_string_vec<'a>(&'a self, str_map: &'a BiMap<String, Tag>) -> Vec<&'a str> {
        let mut x = self.val;
        let mut vec: Vec<&str> = Vec::with_capacity(x.count_ones() as usize);

        // iterate through bits of the u64
        for i in 0..(mem::size_of::<Tags>()*8) as Tag {
            // check if odd number: current first bit is 1
            if !x.is_multiple_of(2) {
                match str_map.get_by_right(&{ i }) {
                    Some(label) => vec.push(label),
                    None => panic!("tried to access label that does not exist!"),
                }
            }
            // shift the u64 bitise to rotate though it
            x >>= 1;
        }
        vec    
    }


    /// compares the value of the tags with another value (marker) with a bit-wise and,
    /// returns true if the result is greater than 0:
    /// `00101 & 01000 -> false`
    /// `00101 & 00100 -> true`
    pub fn bit_and(&self, marker: Marker) -> bool {
        (self.val & marker) > 0
    }

    /// compares the value of the tags with another value (marker) with a bit-wise and
    /// counts the overlaps:
    /// `00101 & 01000 -> 0`
    /// `00101 & 00100 -> 1`
    /// `00101 & 00101 -> 2`
    pub fn bit_and_dist(&self, marker: Marker) -> usize {
        (self.val & marker).count_ones() as usize
    }

    /// get an iterator over the tags in the [`Tags`]
    pub fn iter(&self) -> TagsIterator {
        TagsIterator::new(*self)
    }

    /// get the memory of the [`Tags`] (depends on activated features)
    pub fn mem(&self) -> usize {
        mem::size_of::<Marker>()
    }
}

impl fmt::Debug for Tags {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.to_tag_vec())
    }
}

pub struct TagsIterator {
    tags: Tags,
    i: Tag
}

impl TagsIterator  {
    fn new(tags: Tags) -> TagsIterator {
        TagsIterator {tags, i: 0}
    }
}

impl Iterator for TagsIterator {
    type Item = Tag;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.i as usize == mem::size_of::<Tags>()*8 { return None }
            let result = !self.tags.val.is_multiple_of(2);
            self.tags.val >>= 1;
            self.i += 1;
            if result { return Some(self.i - 1); }
        }
    }
}

pub struct TagsFormatter<'a> {
    tags: Tags,
    translator: &'a Translator
}

impl<'a> TagsFormatter<'a> {
    pub fn new(tags: Tags, translator: &'a Translator) -> TagsFormatter<'a> {
        TagsFormatter {
            tags,
            translator
        }
    }
}

impl fmt::Display for TagsFormatter<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(tag_translator) = self.translator.tag_translator() {
            let tag_vec = self.tags.to_string_vec(tag_translator);

            writeln!(f, "samples:")?;

            for label in tag_vec.into_iter() {
                writeln!(f, "{}", label)?
            }
        } else {
            write!(f, "samples: {:?}", self.tags.to_tag_vec())?
        }

        Ok(())
    }
}

pub struct TagsCountsFormatter<'a> {
    tags: Tags,
    counts: &'a [u32],
    translator: &'a Translator
}

impl<'a> TagsCountsFormatter<'a> {
    pub fn new(tags: Tags, counts: &'a [u32], translator: &'a Translator) -> TagsCountsFormatter<'a> {
        TagsCountsFormatter {
            tags,
            counts,
            translator
        }
    }
}

impl fmt::Display for TagsCountsFormatter<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "{:<20} - counts", "samples")?;

        if let Some(tag_translator) = self.translator.tag_translator() {
            let label_vec = self.tags.to_string_vec(tag_translator);

            for (label, count) in label_vec.into_iter().zip(self.counts) {
                writeln!(f, "{:<20} - {}", label, count)?
            }
        } else {
            let tag_vec = self.tags.to_tag_vec();

            for (tag, count) in tag_vec.into_iter().zip(self.counts) {
                writeln!(f, "{:<20} - {}", tag, count)?
            }
        }


        Ok(())
    }
}


// would be more intuitive with left and right switched but Exts were built this way
/// multiplicities or coverage for each of the 8 possible edges
/// indices: 
/// 0: T right
/// 1: G right
/// 2: C right
/// 3: A right
/// 4: T left
/// 5: G left
/// 6: C left
/// 7: A left
#[derive(PartialEq, PartialOrd, Eq, Ord, Serialize, Deserialize, Clone)]
pub struct EdgeMult {
    edge_mults: [u32; 2*ALPHABET_SIZE],
}

impl Debug for EdgeMult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let edge_f = ["A:", ", C:", ", G:", ", T:", " | A:", ", C:", ", G:", ", T:"];
        for (ef, em) in edge_f.iter().zip(self.edge_mults.iter().rev()) {
             write!(f, "{} {}", ef, em)?
        }
        Ok(())
    }
}

impl Display for EdgeMult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let base = ["A", "C", "G", "T"];
        for (i, b) in (0..ALPHABET_SIZE).rev().zip(base) {
            writeln!(f, "{}: {} | {}", 
                b, 
                self.edge_mults[i + ALPHABET_SIZE], 
                self.edge_mults[i]
            )?
        }

        Ok(())
    }
}

impl EdgeMult {
    /// a new, empty `EdgeMult`
    pub fn new() -> Self {
        EdgeMult { edge_mults: [0;  2 * ALPHABET_SIZE] }
    }

    /// a new `EdgeMult` with values
    pub fn new_from(edge_mults: [u32; 2*ALPHABET_SIZE]) -> Self {
        EdgeMult { edge_mults }
    }

    /// add a count to the an edge
    pub fn add(&mut self, base: u8, dir: Dir, count: u32) {
        self.edge_mults[dir.index(base) as usize] += count
    }

    /// remove an edge from the `EdgeMult`
    pub fn remove(&mut self, base: u8, dir: Dir) {
        self.edge_mults[dir.index(base) as usize] = 0;
    }

    /// add an [`Exts`] to the `EdgeMult`
    pub fn add_exts(&mut self, exts: Exts) {
        let mut exts = exts.val;
        for index in (0..(2 * ALPHABET_SIZE)).rev() {
            if !exts.is_multiple_of(2) {
                self.edge_mults[index] += 1
            }
            exts >>= 1;
        }
    }

    /// get the edge multiplicities as an array 
    pub fn edge_mults(&self) -> [u32; 2*ALPHABET_SIZE] {
        self.edge_mults
    }

    /// get the edge multiplicities to the right of the node
    pub fn right(&self) -> &[u32] {
        &self.edge_mults[(Dir::Right).index_range()]
    }

    /// get the edge multiplicities to the left of the node
    pub fn left(&self) -> &[u32] {
        &self.edge_mults[(Dir::Left).index_range()]
    }

    /// get the sum of all edges of the node
    pub fn sum(&self) -> u32 {
        self.edge_mults.iter().sum::<u32>()
    }

    /// get the multiplicity of a certain edge
    pub fn edge_mult(&self, base: u8, dir: Dir) -> u32 {
        self.edge_mults[dir.index(base) as usize]
    }

    /// get the [`Exts`] corresponding to the `EdgeMult`
    pub fn exts(&self) -> Exts {
        let mut exts_val = 0u8;
        for (i, edge) in self.edge_mults.iter().rev().enumerate() {
            if *edge > 0 { exts_val += 2u8.pow(i as u32) }
        }

        Exts::new(exts_val)
    }

    /// clean the edge mults by removing edge counts that led to filtered kmers
    /// based on a correct [`Exts`]
    pub fn clean_edges(&mut self, exts: Exts) {
        let mut exts = exts.val;
        for index in (0..(2 * ALPHABET_SIZE)).rev() {
            if exts.is_multiple_of(2) {
                self.edge_mults[index] = 0;
            }
            exts >>= 1;
        }
        
    }

    pub fn rc(&mut self) {
        self.edge_mults.reverse();
    }

    pub fn combine(left: &EdgeMult, right: &EdgeMult) -> EdgeMult {
        let mut combined = [0u32; 2*ALPHABET_SIZE];
        (0..ALPHABET_SIZE).for_each(|i| combined[i] = right.edge_mults[i]);
        (ALPHABET_SIZE..(2*ALPHABET_SIZE)).for_each(|i| combined[i] = left.edge_mults[i]);

        EdgeMult::new_from(combined)
    }

    pub fn from_single_dirs(left: &Option<SingleDirEdgeMult>, right: &Option<SingleDirEdgeMult>) -> Option<EdgeMult> {
        if let Some(l_em) = left {
            if let Some(r_em) = right {
                let mut combined = [0u32; 2*ALPHABET_SIZE];
                (0..ALPHABET_SIZE).for_each(|i| combined[i] = r_em.edge_mults[i]);
                (0..ALPHABET_SIZE).for_each(|i| combined[i + ALPHABET_SIZE] = l_em.edge_mults[i]);
        
                return Some(EdgeMult::new_from(combined))
            }
        }

        None
    }

    pub fn single_dir(&self, dir: Dir) -> SingleDirEdgeMult {
        match dir {
            Dir::Left => SingleDirEdgeMult::new(self.edge_mults[ALPHABET_SIZE..(2*ALPHABET_SIZE)]
                .try_into().expect("Error: slice has incorrect length")),
            Dir::Right => SingleDirEdgeMult::new(self.edge_mults[0..ALPHABET_SIZE]
                .try_into().expect("Error: slice has incorrect length")),
        }
    }
}

impl Default for EdgeMult {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct SingleDirEdgeMult {
    edge_mults: [u32; ALPHABET_SIZE]
}

impl SingleDirEdgeMult {
    pub fn new(edge_mults: [u32; ALPHABET_SIZE]) -> Self {
        SingleDirEdgeMult { edge_mults }
    }

    pub fn complement(&self) -> Self {
        let mut reverse = self.edge_mults;
        reverse.reverse();
        SingleDirEdgeMult::new(reverse)
    }

    /// get the multiplicity of a certain edge
    pub fn edge_mult(&self, base: u8) -> u32 {
        self.edge_mults[(ALPHABET_SIZE as u8 - 1 - base) as usize]
    }
}

// would be more intuitive with left and right switched but Exts were built this way
/// mapped transcript/gene/chromosome IDs for each of the 8 possible edges
/// indices: 
/// 0: T right
/// 1: G right
/// 2: C right
/// 3: A right
/// 4: T left
/// 5: G left
/// 6: C left
/// 7: A left
#[derive(PartialEq, PartialOrd, Eq, Ord, Serialize, Deserialize, Clone)]
pub struct EdgeMap {
    edge_maps: [Box<[ID]>; 2*ALPHABET_SIZE],
}

impl EdgeMap {
    /// create a new [`EdgeMap`] by supplying the underlying mapped IDs.
    pub fn new(edge_maps: [Box<[ID]>; 2*ALPHABET_SIZE]) -> EdgeMap {
        EdgeMap { edge_maps }
    }

    /// get the IDs mapped to the specified edge
    fn edge_map(&self, base: u8, dir: Dir) -> &[ID] {
        &self.edge_maps[dir.index(base) as usize]
    }

    /// set the IDs mapped to the edge at the specified index
    fn set_edge_map_at_index(&mut self, edge_map: Box<[ID]>, index: usize) {
        self.edge_maps[index] = edge_map;
    }

    /// add an ID to the edge at the specified index
    fn add_id_to_edge_map_at_index(&mut self, id: ID, index: usize) {
        let mut em = self.edge_maps[index].to_vec();
        em.push(id);

        self.set_edge_map_at_index(em.into(), index);
    }

    /// add an ID at an [`Exts`] to the [`EdgeMap`]
    pub fn add_id(&mut self, exts: Exts, id: ID) {
        let mut exts = exts.val;
        for index in (0..(2 * ALPHABET_SIZE)).rev() {
            if !exts.is_multiple_of(2) {
                self.add_id_to_edge_map_at_index(id, index);
            }
            exts >>= 1;
        }
    }

    /// returns true if the [`EdgeMap`] does not contain any IDs
    pub fn is_empty(&self) -> bool {
        let mut empty = true;

        for emap in self.edge_maps.iter() {
            if !emap.is_empty() { empty = false }
        }

        empty
    }

    /// returns the heap memory used by the [`EdgeMap`] - the stack memory size
    /// is always 128 bytes (8 edges * (8 byte pointer + 8 byte length))
    pub fn mem_heap(&self) -> usize {
        let mut heap = 0;

        for emap in self.edge_maps.iter() {
            heap += mem::size_of_val(&**emap);
        }

        heap
    }

    /// create an [`EdgeMap`] from two [`SingleDirEdgeMap`]s
    pub fn from_single_dirs(left: &Option<SingleDirEdgeMap>, right: &Option<SingleDirEdgeMap>) -> Option<EdgeMap> {
        if let Some(l_em) = left {
            if let Some(r_em) = right {
                let mut combined = EdgeMap::default().edge_maps;
                (0..ALPHABET_SIZE).for_each(|i| combined[i] = r_em.edge_maps[i].clone());
                (0..ALPHABET_SIZE).for_each(|i| combined[i + ALPHABET_SIZE] = l_em.edge_maps[i].clone());
        
                return Some(EdgeMap::new(combined))
            }
        }

        None
    }

    /// get the [`SingleDirEdgeMap`] in the specified [`Dir`]
    pub fn single_dir(&self, dir: Dir) -> SingleDirEdgeMap {
        let singe_dir: &[Box<[ID]>; 4] = match dir {
            Dir::Left => self.edge_maps[ALPHABET_SIZE..(2*ALPHABET_SIZE)].try_into().expect("Error: slice has incorrect length"),
            Dir::Right => self.edge_maps[0..ALPHABET_SIZE].try_into().expect("Error: slice has incorrect length"),
        };

        SingleDirEdgeMap::new(singe_dir.clone())
    }

    /// clean the edge maps by removing IDs mapped to edges that led to filtered kmers,
    /// based on a correct [`Exts`]
    pub fn clean_edges(&mut self, exts: Exts) {
        let mut exts = exts.val;
        for index in (0..(2 * ALPHABET_SIZE)).rev() {
            if exts.is_multiple_of(2) {
                self.edge_maps[index] = [].into();
            }
            exts >>= 1;
        }
        
    }
}

impl Default for EdgeMap {
    fn default() -> Self {
        let edge_maps = array::from_fn(|_n| Vec::new().into());
        Self { edge_maps }
    }
}

impl Debug for EdgeMap {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let edge_f = ["A:", ", C:", ", G:", ", T:", " | A:", ", C:", ", G:", ", T:"];
        for (ef, em) in edge_f.iter().zip(self.edge_maps.iter().rev()) {
             write!(f, "{} {:?}", ef, em)?
        }
        Ok(())
    }
}

impl Display for EdgeMap {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let base = ["A", "C", "G", "T"];
        for (i, b) in (0..ALPHABET_SIZE).rev().zip(base) {
            writeln!(f, "{}: {:?} | {:?}", 
                b, 
                self.edge_maps[i + ALPHABET_SIZE], 
                self.edge_maps[i]
            )?
        }

        Ok(())
    }
}

/// mapped transcript/gene/chromosome IDs for each of the 4 possible edges
/// indices in one direction: 
/// 
/// 0: T 
/// 1: G 
/// 2: C 
/// 3: A 
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct SingleDirEdgeMap {
    edge_maps: [Box<[ID]>; ALPHABET_SIZE]
}

impl SingleDirEdgeMap {
    pub fn new(edge_maps: [Box<[ID]>; ALPHABET_SIZE]) -> Self {
        SingleDirEdgeMap { edge_maps }
    }

    pub fn complement(&self) -> Self {
        let mut reverse = self.edge_maps.clone();
        reverse.reverse();
        SingleDirEdgeMap::new(reverse)
    }

    /// get the IDs mapped to the specified edge
    pub fn edge_map(&self, base: u8) -> &[ID] {
        &self.edge_maps[(ALPHABET_SIZE as u8 - 1 - base) as usize]
    }
}

// TODO add methods
#[derive(Debug, Serialize, Deserialize)]
pub struct Label {
    group: char,
    sample_label: String
}

/// category for the 
#[derive(Debug, Deserialize, Serialize, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
#[serde(rename_all = "kebab-case")]
pub enum BaseQuality {
    NoCall,
    Marginal,
    Medium,
    High
}

impl BaseQuality {
    fn from_u64(quality: u64) -> BaseQuality {
        match quality {
            0 => Self::NoCall,
            1 => Self::Marginal,
            2 => Self::Medium,
            3 => Self::High,
            _ => panic!("invalid base quality value")
        }
    }

    fn as_char(&self) -> char {
        match self {
            Self::NoCall => '#',
            Self::Marginal => '-',
            Self::Medium => ';',
            Self::High => 'C',
        }
    }
}

impl fmt::Display for BaseQuality {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoCall => write!(f, "no-call"),
            Self::Marginal => write!(f, "marginal"),
            Self::Medium => write!(f, "medium"),
            Self::High => write!(f, "high"),
        }
    }
}

#[derive(Debug, Deserialize, Serialize, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct QualityBins {
    marginal_top_cutoff: u8,
    high_bottom_cutoff: u8,
}

impl Default for QualityBins  {
    fn default() -> Self {
        Self { marginal_top_cutoff: 15, high_bottom_cutoff: 30 }
    }
}

impl QualityBins {
    pub fn new(marginal_top_cutoff: u8, high_bottom_cutoff: u8) -> QualityBins {
        Self { marginal_top_cutoff, high_bottom_cutoff }
    }

    fn base_quality(&self, score: u8) -> BaseQuality {
        if score <= 2 {
            BaseQuality::NoCall
        } else if score < self.marginal_top_cutoff {
            BaseQuality::Marginal
        } else if score > self.high_bottom_cutoff {
            BaseQuality::High
        } else {
            BaseQuality::Medium
        }
    }

    fn base_quality_from_ascii_bytes(&self, char: u8) -> BaseQuality {
        let score = char - 33;
        self.base_quality(score)
    }
}

#[derive(Deserialize, Serialize, Clone, PartialEq, PartialOrd)]
pub struct QualityVec {
    storage: Vec<BaseQuality>
}

impl QualityVec {
    fn from_vec(quality_vec: Vec<BaseQuality>) -> QualityVec {
        QualityVec { storage: quality_vec }
    }

    fn from_ascii_bytes(quality_scores: &[u8], quality_bins: QualityBins) -> QualityVec {
        let mut vec = Vec::new();
        for score in quality_scores {
            vec.push(quality_bins.base_quality_from_ascii_bytes(*score));
        }

        QualityVec { storage: vec }
    }

    fn iter_k_lowest_q<K: Kmer>(&'_ self) -> KLowestQualityIter<'_, K> {
        KLowestQualityIter { 
            quality_vec: self, 
            start_pos: 0, 
            phantom_data: PhantomData,
        }
    }

    fn iter_k_random_q<K: Kmer>(&'_ self) -> KRandomQualityIter<'_, K> {
        KRandomQualityIter { 
            quality_vec: self, 
            start_pos: 0, 
            phantom_data: PhantomData,
            rng: rand::thread_rng()
        }
    }


    fn len(&self) -> usize {
        self.storage.len()
    }
}

impl fmt::Debug for QualityVec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for q in self.storage.iter() {
            write!(f, "{}", q.as_char())?;
        }

        Ok(())
    }
}

pub struct KLowestQualityIter<'a, K: Kmer> {
    quality_vec: &'a QualityVec,
    start_pos: usize,
    phantom_data: PhantomData<K>
}

impl<K: Kmer> Iterator for KLowestQualityIter<'_, K> {
    type Item = BaseQuality;

    fn next(&mut self) -> Option<Self::Item> {
        let end_pos = self.start_pos + K::k();
        if end_pos <= self.quality_vec.len() {
            let range = self.start_pos..end_pos;

            let quality = self.quality_vec.storage[range]
                .iter()
                .min()
                .expect("missing base quality");

            self.start_pos += 1;
            Some(*quality)
        } else {
            None
        }
    }
}

pub struct KRandomQualityIter<'a, K: Kmer> {
    quality_vec: &'a QualityVec,
    start_pos: usize,
    phantom_data: PhantomData<K>,
    rng: ThreadRng,

}

impl<K: Kmer> Iterator for KRandomQualityIter<'_, K> {
    type Item = BaseQuality;

    fn next(&mut self) -> Option<Self::Item> {
        let end_pos = self.start_pos + K::k();
        if end_pos <= self.quality_vec.len() {
            let quality = BaseQuality::from_u64(self.rng.gen_range(1, 4) as u64);
            Some(quality)
        } else {
            None
        }
    }
}

/// add the alignment buffer to a structure
/// - `size_heap`: contents of boxes/vectors -> are stored separately and do not go into alignment calculation
pub fn size_aligned(size_stack: usize, size_heap: usize, align: usize) -> usize {
    let empty = size_stack % align;
    let buffer = match empty {
        0 => 0,
        _ => align - empty
    };

    buffer + size_heap + size_stack
}

pub fn build_test_graph<K, SD, DI>() -> (SerReads<DI>, SerKmers<K, SD>, SerGraph<K, SD>)
where
    K: Kmer +  Send + Sync,
    SD: SummaryData<DI>,
    DI: ReadData
{
    /*
    transcrips
    GCAGCTAGCTAGCGCGACTACGATCGTAGCGCAGCGAGCAGGGGGGGGGATAGCTGTCGCGGGGACGTATTATTATTAAAATTGCGGCGCGAGCTATTCGAGCGGAGCGAGCGACAGGAGCGGAGTTTGCGGTACGGGATTTTCGGATATCGGC
    GCGATTATTTTGCGGGGGATTTTCGGTAGCGACTGGGGGGGGGTATCGATCGTGACAGCTTTCGACTGGGAGCGCAGCTAGGCAGGACGCATTAATTATATATCATTATTTTTTTCTATAAAAAAAAAAGAGCTAGCGATCGACGCGATCGAC
    TATATTATCGGCTGAGCGAGCGGGGGCAGCTATATTACGCGATAAAGAGCCCCCCGAGGCGAGGCGGACTTACGTAGCGCAGGCACCATGACGAGCTAGCAGTCAGTCGTAGCGATCA
    GCTAGCTAGCTGACTACGATCGACGGGGAGCATTAATTAGAAAAAAGAGAGAGACAGCTTTCGACTGGGAGCGCAGCTAGGCAGGACGCATTACTATCTATTATTATATATCATTATTTGCGATTGGGGTGCTAGCATGCGT
    */

    let raw_reads = [
        ["GCAGCTAGCTAGCGCGACTACGATCGTAGCGCAGCGAGCAGGGGGGGGGA", "gene1", "sample1", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;;CCCCCCCCCCCCCCC"],
        ["CTAGCGCGACTACGATCGTAGCGCAGCGAGCAGGGGGGGGGATAGCTGTC", "gene1", "sample1", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["GCAGCGAGCAGGGGGGGGGATAGCTGTCGCGGGGACGTATTATTATTAAA", "gene1", "sample1", "CCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["CGGGGACGTATTATTATTAAAATTGCGGCGCGAGCTATTCGAGCGGAGCG", "gene1", "sample1", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["TATTCGAGCGGAGCGAGCGACAGGAGAGGAGTTTGCGGTACGGGATTTTC", "gene1", "sample1", "CCCCCCCCCCCCCCCCCCCCCCCCCC--CCCCCCCCCCCCCCCCCCCCCC"],
        ["CGGAGCGAGCGACAGGAGCGGAGTTTGCGGTACGGGATTTTCGGATATCG", "gene1", "sample1", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["GAGCGAGCGACAGGAGCGGAGTTTGCGGTACGGGATTTTCGGATATCGGC", "gene1", "sample1", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["GCAGCTAGCTAGCGCGACTACGATCGTAGCGCAGCGAGCAGGGGGGGGGA", "gene1", "sample1", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCC"],
        ["TATTATTAAAATTGCGGCGCGAGCTATTCGAGCGGAGCGAGCGACAGGAG", "gene1", "sample1", "CCCCCCCCCCCCCCCCCC-CC--CCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["ACTACGATCGTAGCGCAGCGAGCAGGGGGGGGGATAGCTGTCGCGGGGAC", "gene1", "sample1", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["AGCGGAGCGAGCGACAGGAGCGGAGTTTGCGGTACGGGATTTTCGGATAT", "gene1", "sample1", "CCCCCCCCCCCCCCCCCCCCC--CCCC--CCCCCCCCCCCCCCCCCCCCC"],
        ["TAGCTAGCGCGACTACGATCGTAGCGCAGCGAGCAGGGGGGGGGATAGCT", "gene1", "sample1", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["TACGATCGTAGCGCAGCGAGCAGGGGGGGGGATAGCTGTCGCGGGGACGT", "gene1", "sample1", "CCCCCCCCCCCCCCCCC--CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["ATTGCGGCGCGAGCTATTCGAGCGGAGCGAGCGACAGGAGCGGAGTTTGC", "gene1", "sample1", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["CAGCTAGCTAGCGCGACTACGATCGTAGCGCAGCGAGCAGGGGGGGGGAT", "gene1", "sample1", "CCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["GCGATTATTTTGCGGGGGATTTTCGGTAGCGACTGGGGGGGGGTATCGAT", "gene2", "sample2", "CCCCCCCC---CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["GGGGATTTTCGGTAGCGACTGGGGGGGGGTATCGATCGTGACAGCTTTCG", "gene2", "sample2", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCC"],
        ["TTTCGACTGGGAGCGCAGCTAGGCAGGACGCATTACTATCTATTATTATA", "gene2", "sample2", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC.-CCCCCCCCCCCCCCCCCC"],
        ["CAGGACGCATTACTATCTATTATTATATATCATTATTTTTTTCTATAAAA", "gene2", "sample2", "CCCCCCCCCCCCCCCC---CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["TTTCGGTAGCGACTGGGGGGGGGTATCGATCGTGACAGCTTTCGACTGGG", "gene2", "sample2", "CCCCCCCCCCCCCC-CCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["AGCTAGGCAGGACGCATTACTATCTATTATTATATATCATTATTTTTTTC", "gene2", "sample2", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["AGGACGCATTACTATCTATTATTATATATCATTATTTTTTTCTATAAAAA", "gene2", "sample2", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC--;CCCCCCCCCCCCCCCCCC"],
        ["GGGATTTTCGGTAGCGACTGGGGGGGGGTATCGATCGTGACAGCTTTCGA", "gene2", "sample2", "CCCCCCCCCCCC---CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["CTTTCGACTGGGAGCGCAGCTAGGCAGGACGCATTACTATCTATTATTAT", "gene2", "sample2", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["CGACTGGGGGGGGGTATCGATCGTGACAGCTTTCGACTGGGAGCGCAGCT", "gene2", "sample2", "CCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCC"],
        ["ATTATATATCATTATTTTTTTCTATAAAAAAAAAAGAGCTAGCGATCGAC", "gene2", "sample2", "CCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["TAATTATATATCATTATTTTTTTCTATAAAAAAAAAAGAGCTAGCGATCG", "gene2", "sample2", "CCCCCC-CCCCCCCCCCCCCCCCCCCCCCCC--CCCCCCCCCCCCCCCCC"],
        ["ATCATTATTTTTTTCTATAAAAAAAAAAGAGCTAGCGATCGACGCGATCG", "gene2", "sample2", "CCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["CGTGACAGCTTTCGACTGGGAGCGCAGCTAGGCAGGACGCATTAATTATA", "gene2", "sample2", "CCCCCCCCCCCCCCCCCC;;;CCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["GGGGGTATCGATCGTGACAGCTTTCGACTGGGAGCGCAGCTAGGCAGGAC", "gene2", "sample2", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["TATATTATCGGCTGAGCGAGCGGGGGGAGCTATATTACGCGATAAAGAGC", "gene3", "sample3", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["AGCGAGCGGGGGCAGCTATATTACGCGATAAAGAGCCCCCCGAGGCGAGG", "gene3", "sample3", "CCCCCCCCCCCCCCCCCCCCCCCCC---CCCCCCCCCCCCCCCCCCCCCC"],
        ["CTATATTACGCGATAAAGAGCCCCCCGAGGCGAGGCGGACTTACGTAGCG", "gene3", "sample3", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["CGCGATAAAGAGCCCCCCGAGGCGAGGCGGACTTACGTAGCGCAGGCACC", "gene3", "sample3", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["CTTACGTAGCGCAGGCACCATGACGAGCTAGCAGTCAGTCGTAGCGATCA", "gene3", "sample3", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["GCGATAAAGAGCCCCCCGAGGCGAGGCGGACTTACGTAGCGCAGGCACCA", "gene3", "sample3", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["TATATTACGCGATAAAGAGCCCCCCGAGGCGAGGCGGACTTACGTAGCGC", "gene3", "sample3", "CCCCCCCCCCCCCCCCCC-----CCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["GAGCGGGGGCAGCTATATTACGCGATAAAGAGCCCCCCGAGGCGAGGCGG", "gene3", "sample3", "CCCCCCCCCCCCCCCCCCCCCCC;;;CCCCC;;CCCCCCCCCCCCCCCCC"],
        ["AGAGCCCCCCGAGGCGAGGCGGACTTACGTAGCGCAGGCACCATGACGAG", "gene3", "sample3", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["GGCGGACTTACGTAGCGCAGGCACCATGACGAGCTAGCAGTCAGTCGTAG", "gene3", "sample3", "CCCCCCCCCCCCC---CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["ACTTACGTAGCGCAGGCACCATGACGAGCTAGCAGTCAGTCGTAGCGATC", "gene3", "sample3", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["GGGGCAGCTATATTACGCGATAAAGAGCCCCCCGAGGCGAGGCGGACTTA", "gene3", "sample3", "CCCCCCCCCCCCCCCCCCCCCCCCCCCC---CC-CCCCCCCCCCCCCCCC"],
        ["CCGAGGCGAGGCGGACTTACGTAGCGCAGGCACCATGACGAGCTAGCAGT", "gene3", "sample3", "CCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["GGCGAGGCGGACTTACGTAGCGCAGGCACCATGACGAGCTAGCAGTCAGT", "gene3", "sample3", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["GAGCCCCCCGAGGCGAGGCGGACTTACGTAGCGCAGGCACCATGACGAGC", "gene3", "sample3", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["GCTAGCTAGCTGACTACGATCGACGGGGAGCATTAATTAGAAAAAAGAGA", "gene4", "sample4", "CCCCCCCCCCCCCCCCCC--CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["ACTATCTATTATTATATATCATTATTTGCGATTGGGGTGCTAGCATGCGT", "gene4", "sample4", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["AGGCAGGACGCATTACTATCTATTATTATATATCATTATTTGCGATTGGG", "gene4", "sample4", "CCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCC"],
        ["ATCGACGGGGAGCATTAATTAGAAAAAAGAGAGAGACAGCTTTCGACTGG", "gene4", "sample4", "CCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["AGGACGCATTACTATCTATTATTATATATCATTATTTGCGATTGGGGTGC", "gene4", "sample4", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["GAGCATTAATTAGAAAAAAGAGAGAGACAGCTTTCGACTGGGAGCGCAGC", "gene4", "sample4", "CCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCC"],
        ["ACTGGGAGCGCAGCTAGGCAGGACGCATTACTATCTATTATTATATATCA", "gene4", "sample4", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["GATCGACGGGGAGCATTAATTAGAAAAAAGAGAGAGACAGCTTTCGACTG", "gene4", "sample4", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["CGCAGCTAGGCAGGACGCATTACTATCTATTATTATATATCATTATTTGC", "gene4", "sample4", "CCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["CGATCGACGGGGAGCATTAATTAGAAAAAAGAGAGAGACAGCTTTCGACT", "gene4", "sample4", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["ATCGACGGGGAGCATTAATTAGAAAAAAGAGAGAGACAGCTTTCGACTGG", "gene4", "sample4", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC--CCCCCCCCCCCCCC"],
        ["GCAGCTAGGCAGGACGCATTACTATCTATTATTATATATCATTATTTGCG", "gene4", "sample4", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["GACGGGGAGCATTAATTAGAAAAAAGAGAGAGACAGCTTTCGACTGGGAG", "gene4", "sample4", "CCCC----CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"],
        ["CGATCGACGGGGAGCATTAATTAGAAAAAAGAGAGAGACAGCTTTCGACT", "gene4", "sample4", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC---CCCCCCCCCCCCCC"],
        ["TAGGCAGGACGCATTACTATCTATTATTATATATCATTATTTGCGATTGG", "gene4", "sample4", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCC"]
    ];

    let mut reads = Reads::new_with_quality(crate::reads::Strandedness::Forward);
    let mut id_translator = BiMap::new();
    let mut tag_translator = BiMap::new();

    for [seq, gene, sample, quality] in raw_reads {
        let gene = String::from(gene);
        let sample = String::from(sample);

        let id = match id_translator.get_by_left(&gene) {
            Some(id) => *id,
            None =>  {
                let new_id = id_translator.len() as ID;
                id_translator.insert(gene, new_id);
                new_id
            }
        };

        let tag = match tag_translator.get_by_left(&sample) {
            Some(tag) => *tag,
            None =>  {
                let new_tag = tag_translator.len() as Tag;
                tag_translator.insert(sample, new_tag);
                new_tag
            }
        };

        reads.add_read(DnaString::from_acgt_bytes(seq.as_bytes()), None, DI::new(id, tag), Some(quality.as_bytes()));
    }

    let translator = Translator::new(id_translator, tag_translator);
    let reads_paired = ReadsPaired::Unpaired { reads };
    let sample_kmers = reads_paired.tag_kmers_vec(K::k(), 4);
    let ser_reads = SerReads::new(reads_paired, translator.clone());
    
    let sample_info = SampleInfo::new(0b1100, 0b0011, sample_kmers);
    let summary_config = SummaryConfig::new(sample_info);

    let (kmers, _) = filter_kmers::<SD, K, _>(
        ser_reads.reads(), 
        &summary_config, 
        false, 
        5., 
        false
    );

    let ser_kmers = SerKmers::new(kmers.clone(), translator.clone(), summary_config.clone());

    let comp_spec = CheckCompress::new(|d: SD, _| d, |d, d1| d.join_test(d1));
    let graph = compress_kmers_with_hash(true, &comp_spec, kmers, false, false).finish();
    let ser_graph = SerGraph::new(graph, translator, summary_config);

    (ser_reads, ser_kmers, ser_graph)
}

#[cfg(test)]
mod tests {
    use bimap::BiMap;

    use crate::{ALPHABET_SIZE, BaseQuality, Dir, EdgeMap, EdgeMult, Exts, Kmer, QualityBins, Tags, TagsCountsFormatter, TagsFormatter, kmer::{Kmer4, Kmer17}, size_aligned, summarizer::{ID, Marker, Tag, Translator}};

    #[test]
    fn test_dir_index() {
        let bases = [0, 1, 2, 3];

        for (i, (dir, base)) in [Dir::Right, Dir::Left].iter().flat_map(|d| bases.into_iter().rev().map(move |b| (d, b))).enumerate() {
            assert_eq!(i as u8, dir.index(base))
        }

        assert_eq!((Dir::Left).index_range(), 4..8);
        assert_eq!((Dir::Right).index_range(), 0..4);
    }

    #[test]
    fn test_remove_ext() {
        let ext = Exts::new(0b11111011);
        assert_eq!(ext.remove(Dir::Right, 0).val, 0b11101011);
        assert_eq!(ext.remove(Dir::Left, 1).val, 0b11111001);
    }

    #[test]
    fn test_edge_mult() {
        let mut edge_mult = EdgeMult::new();
        assert_eq!(edge_mult.edge_mults, [0; 2*ALPHABET_SIZE]);

        let exts = Exts::new(0b11111111);
        edge_mult.add_exts(exts);
        assert_eq!(edge_mult.edge_mults, [1, 1, 1, 1, 1, 1, 1, 1]);

        edge_mult.add(0, crate::Dir::Right, 2);
        edge_mult.add(2, crate::Dir::Left, 78989);
        assert_eq!(edge_mult.edge_mults, [1, 1, 1, 3, 1, 78990, 1, 1]);

        let exts = Exts::new(0);
        let comp = [1, 1, 1, 3, 1, 78990, 1, 1];
        edge_mult.add_exts(exts);
        assert_eq!(edge_mult.edge_mults, comp);
        assert_eq!(edge_mult.sum(), comp.iter().sum::<u32>());

        let exts = Exts::new(0b10101010);
        edge_mult.add_exts(exts);
        assert_eq!(edge_mult.edge_mults, [2, 1, 2, 3, 2, 78990, 2, 1]);

        let exts = Exts::new(0b01010101);
        edge_mult.add_exts(exts);
        assert_eq!(edge_mult.edge_mults, [2, 2, 2, 4, 2, 78991, 2, 2]);

        let clean_exts = Exts::new(0b10010101);
        edge_mult.clean_edges(clean_exts);
        assert_eq!(edge_mult.edge_mults, [2, 0, 0, 4, 0, 78991, 0, 2]);

        let mut em = EdgeMult::new();
        let exts = Exts::new(0b00011101);
        println!("{:?}", exts);
        em.add_exts(exts);
        println!("{}", em);
        println!("{:?}", em);
        assert_eq!(em.exts(), exts);

        assert_eq!(em.left(), &[1, 1, 0, 1]);
        assert_eq!(em.right(), &[0, 0, 0, 1]);

        // single dir em
        let sdir_em = em.single_dir(Dir::Left);
        assert_eq!(sdir_em.edge_mults, [1, 1, 0, 1]);
        let sdir_em = em.single_dir(Dir::Right);
        assert_eq!(sdir_em.edge_mults, [0, 0, 0, 1]);
        let em_2 = EdgeMult::from_single_dirs(&Some(em.single_dir(Dir::Left)), &Some(em.single_dir(Dir::Right)));
        assert_eq!(em_2.unwrap(), em);
        assert_eq!(sdir_em.complement().edge_mults, [1, 0, 0, 0]);

        // reverse complement
        em.rc();
        assert_eq!(em.edge_mults, [1, 0, 1, 1, 1, 0, 0, 0]);
    }

    #[test]
    fn test_edge_map() {
        let empty_emaps: [Box<[ID]>; 8] = Default::default();
        let default_emap = EdgeMap::default();

        let mut emap = EdgeMap::new(empty_emaps.clone());

        assert!(emap.is_empty());

        assert_eq!(emap.edge_maps, empty_emaps);
        assert_eq!(emap, default_emap);

        let m = vec![1, 2];
        emap.set_edge_map_at_index(m.clone().into(), 1); // right G
        let r = emap.edge_map(2, Dir::Right);

        assert!(!emap.is_empty());

        assert_eq!(&m, r);

        emap.add_id_to_edge_map_at_index(3, 1);
        let r = emap.edge_map(2, Dir::Right);
        assert_eq!(&[1, 2, 3], r);

        let exp_mem = 3*std::mem::size_of::<ID>();
        let r_mem = emap.mem_heap();
        assert_eq!(exp_mem, r_mem);

        emap.add_id(Exts::new(0b01000010), 4); // right G and left C
        assert_eq!(&[1, 2, 3, 4], emap.edge_map(2, Dir::Right));
        assert_eq!(&[4], emap.edge_map(1, Dir::Left));

        let left = emap.single_dir(Dir::Left);
        let right = emap.single_dir(Dir::Right);
        let new_emap = EdgeMap::from_single_dirs(&Some(left), &Some(right)).unwrap();
        assert_eq!(emap, new_emap);

        let mut clean_emap = EdgeMap::default();
        clean_emap.add_id_to_edge_map_at_index(4, 6); // left C

        emap.clean_edges(Exts::new(0b00000010));

        assert_eq!(emap, clean_emap);

        assert_eq!("A: [] | []\nC: [4] | []\nG: [] | []\nT: [] | []\n", &format!("{}", emap));
        assert_eq!("A: [], C: [4], G: [], T: [] | A: [], C: [], G: [], T: []", &format!("{:?}", emap));

    }

    #[test]
    fn test_base_quality() {
        let bq = BaseQuality::from_u64(0);
        assert_eq!(&format!("{bq} {}", bq.as_char()), "no-call #");

        let bq = BaseQuality::from_u64(1);
        assert_eq!(&format!("{bq} {}", bq.as_char()), "marginal -");

        let bq = BaseQuality::from_u64(2);
        assert_eq!(&format!("{bq} {}", bq.as_char()), "medium ;");

        let bq = BaseQuality::from_u64(3);
        assert_eq!(&format!("{bq} {}", bq.as_char()), "high C");

        let bins = QualityBins::new(15, 30);
        assert_eq!(bins, QualityBins::default());
    }

    #[test]
    #[should_panic]
    fn test_base_quality_panic() {
        let _bq = BaseQuality::from_u64(5);
    }

    #[test]
    fn test_bit_and_dist() {
        let marker: Marker = 0b1111000011110000111100001111000011110000111100001111000011110000;
        println!("marker:   {:064b}", marker);

        let tags = Tags::from_tag_vec(&vec![0, 1, 4]);
        println!("tags:     {:064b}", tags.val);
        let dist = tags.bit_and_dist(marker);
        println!("dist: {}", dist);
        assert_eq!(tags.len(), 3);

        let tags = Tags::from_tag_vec(&vec![1, 5, 19, 25, 32]);
        println!("tags:     {:064b}", tags.val);
        let dist = tags.bit_and_dist(marker);
        println!("dist: {}", dist);
        assert_eq!(tags.len(), 5);

        let tags = Tags::from_tag_vec(&vec![0, 1, 2, 3, 4, 5, 6, 7, 63]);
        println!("tags:     {:064b}", tags.val);
        let dist = tags.bit_and_dist(marker);
        println!("dist: {}", dist);
        assert_eq!(tags.len(), 9);

        let tags = Tags::from_tag_vec(&vec![31]);
        println!("tags:     {:064b}", tags.val);
        let dist = tags.bit_and_dist(marker);
        println!("dist: {}", dist);
        assert_eq!(tags.len(), 1);

        let tags = Tags::from_tag_vec(&vec![63]);
        println!("tags:     {:064b}", tags.val);
        let dist = tags.bit_and_dist(marker);
        println!("dist: {}", dist);
        assert_eq!(tags.len(), 1);
    }

    #[test]
    fn test_tag_formatter() {
        let mut tag_translator = BiMap::new();
        let samples = vec!["A", "B", "C", "D", "E", "F", "G"];

        for (i, label) in samples.into_iter().enumerate() {
            tag_translator.insert(label.to_string(), i as Tag);
        }

        let translator = Translator::new_tag_translator(tag_translator);

        let tags = Tags::from_tag_vec(&vec![0, 1, 4]);
        let counts = vec![1, 2, 3].into_boxed_slice();
        print!("{}", TagsCountsFormatter::new(tags, &counts, &translator));

        let tags = Tags::from_tag_vec(&vec![0, 1, 4, 6]);
        let counts = vec![1, 2, 3, 0].into_boxed_slice();
        print!("{}", TagsCountsFormatter::new(tags, &counts, &translator));

        let tags = Tags::from_tag_vec(&vec![0, 1, 4]);
        print!("{}", TagsFormatter::new(tags, &translator));

        let tags = Tags::from_tag_vec(&vec![0, 1, 4, 6]);
        print!("{}", TagsFormatter::new(tags, &translator));

    }

    #[test]
    fn test_iter_tags() {
        let tags = Tags::from_tag_vec(&vec![0, 1, 4, 12, 32, 63]);
        for tag in tags.iter() {
            println!("tag: {tag}")
        }
    }

    #[test]
    fn test_kmer_complexity() {
        let kmer = Kmer4::empty();
        assert!(kmer.has_low_complexity());
        let kmer = Kmer4::from_ascii("TTTTTTTT".as_bytes());
        assert!(kmer.has_low_complexity());
        let kmer = Kmer4::from_ascii("CCCCCCCC".as_bytes());
        assert!(kmer.has_low_complexity());
        let kmer = Kmer4::from_ascii("GGGGGGGG".as_bytes());
        assert!(kmer.has_low_complexity());

        let kmer = Kmer4::from_ascii("ACGATCGA".as_bytes());
        assert!(!kmer.has_low_complexity());
        let kmer = Kmer4::from_ascii("AGCAGCTC".as_bytes());
        assert!(!kmer.has_low_complexity());

        let kmer = Kmer17::empty();
        assert!(kmer.has_low_complexity());
        let kmer = Kmer17::from_ascii("TTTTTTTTTTTTTTTTT".as_bytes());
        assert!(kmer.has_low_complexity());
        let kmer = Kmer17::from_ascii("CCCCCCCCCCCCCCCCC".as_bytes());
        assert!(kmer.has_low_complexity());
        let kmer = Kmer17::from_ascii("GGGGGGGGGGGGGGGGG".as_bytes());
        assert!(kmer.has_low_complexity());

        let kmer = Kmer17::from_ascii("ACGATCGAGACTGACTG".as_bytes());
        assert!(!kmer.has_low_complexity());
        let kmer = Kmer17::from_ascii("AGCAGCTCAGCTAGCTG".as_bytes());
        assert!(!kmer.has_low_complexity());
    }
}


