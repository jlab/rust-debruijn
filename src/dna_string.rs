// Copyright 2014 Johannes Köster and 10x Genomics
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! A 2-bit encoding of arbitrary length DNA sequences.
//!
//! Store arbitrary-length DNA strings in a packed 2-bit encoding. Individual base values are encoded
//! as the integers 0,1,2,3 corresponding to A,C,G,T.
//!
//! # Example
//! ```
//! use debruijn::Kmer;
//! use debruijn::dna_string::*;
//! use debruijn::kmer::Kmer16;
//! use debruijn::Vmer;
//!
//! // Construct a new DNA string
//! let dna_string1 = DnaString::from_dna_string("ACAGCAGCAGCACGTATGACAGATAGTGACAGCAGTTTGTGACCGCAAGAGCAGTAATATGATG");
//!
//! // Get an immutable view into the sequence
//! let slice1 = dna_string1.slice(10, 40);
//!
//! // Get a kmer from the DNA string
//! let first_kmer: Kmer16 = slice1.get_kmer(0);
//! assert_eq!(first_kmer, Kmer16::from_ascii(b"CACGTATGACAGATAG"))

use itertools::Itertools;
use serde_derive::{Deserialize, Serialize};
use std::borrow::Borrow;
use std::cmp::min;
use std::collections::hash_map::DefaultHasher;
use std::error::Error;
use std::fmt::{self, Display};
use std::hash::{Hash, Hasher};

use crate::{base_to_bits, base_to_bits_checked};
use crate::bits_to_ascii;
use crate::bits_to_base;
use crate::dna_only_base_to_bits;

use crate::Kmer;
use crate::Mer;
use crate::MerIter;
use crate::Vmer;

const BLOCK_BITS: usize = 64;
const WIDTH: usize = 2;

const MASK: u64 = 0x3;

/// A container for sequence of DNA bases.
/// ```
/// use debruijn::dna_string::DnaString;
/// use debruijn::kmer::Kmer8;
/// use debruijn::{Mer, Vmer};
///
/// let dna_string = DnaString::from_dna_string("ATCGTACGTACGTAGTC");
///
/// // Iterate over 8-mers
/// for k in dna_string.iter_kmers::<Kmer8>() {
///     println!("{:?}", k);
/// }
///
/// // Get a base, encoded as a byte in 0-3 range
/// assert_eq!(dna_string.get(0), 0);
/// assert_eq!(dna_string.get(1), 3);
///
/// // Make a read-only 'slice' of a DnaString
/// let slc = dna_string.slice(1, 10);
///
///  assert_eq!(slc.iter_kmers::<Kmer8>().next(), dna_string.iter_kmers::<Kmer8>().skip(1).next());
/// ```
#[derive(Ord, PartialOrd, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct DnaString {
    storage: Vec<u64>,
    len: usize,
}

impl Mer for DnaString {
    fn len(&self) -> usize {
        self.len
    }

    fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Get the value at position `i`.
    #[inline(always)]
    fn get(&self, i: usize) -> u8 {
        let (block, bit) = self.addr(i);
        self.get_by_addr(block, bit)
    }

    /// Set the value as position `i`.
    fn set_mut(&mut self, i: usize, value: u8) {
        let (block, bit) = self.addr(i);
        self.set_by_addr(block, bit, value);
    }

    fn set_slice_mut(&mut self, _: usize, _: usize, _: u64) {
        unimplemented!()
    }

    fn rc(&self) -> DnaString {
        let mut dna_string = DnaString::new();
        let rc = (0..self.len()).rev().map(|i| 3 - self.get(i));

        dna_string.extend(rc);
        dna_string
    }
}

impl Vmer for DnaString {
    fn new(len: usize) -> Self {
        Self::blank(len)
    }

    fn max_len() -> usize {
        <usize>::MAX
    }

    /// Get the kmer starting at position pos
    fn get_kmer<K: Kmer>(&self, pos: usize) -> K {
        assert!(self.len() - pos >= K::k());

        // Which block has the first base
        let (mut block, _) = self.addr(pos);

        // Where we are in the kmer
        let mut kmer_pos = 0;

        // Where are in the block
        let mut block_pos = pos % 32;

        let mut kmer = K::empty();

        while kmer_pos < K::k() {
            // get relevent bases for current block
            let nb = min(K::k() - kmer_pos, 32 - block_pos);

            let v = self.storage[block];
            let val = v << (2 * block_pos);
            kmer.set_slice_mut(kmer_pos, nb, val);

            // move to next block, move ahead in kmer.
            block += 1;
            kmer_pos += nb;
            // alway start a beginning of next block
            block_pos = 0;
        }

        kmer
    }
}

impl DnaString {
    /// Create an empty DNA string
    pub fn new() -> DnaString {
        DnaString {
            storage: Vec::new(),
            len: 0,
        }
    }

    /// Length of the sequence
    pub fn len(&self) -> usize {
        self.len
    }

    /// Create a new instance with a given capacity.
    pub fn with_capacity(n: usize) -> Self {
        let blocks = ((n * WIDTH) >> 6) + (if (n * WIDTH) & 0x3F > 0 { 1 } else { 0 });
        let storage = Vec::with_capacity(blocks);

        DnaString { storage, len: 0 }
    }

    /// Create a DnaString of length n initialized to all A's
    pub fn blank(n: usize) -> Self {
        let blocks = ((n * WIDTH) >> 6) + (if (n * WIDTH) & 0x3F > 0 { 1 } else { 0 });
        let storage = vec![0; blocks];

        DnaString { storage, len: n }
    }

    /// Create a DnaString corresponding to an ACGT-encoded str.
    pub fn from_dna_string(dna: &str) -> DnaString {
        let mut dna_string = DnaString {
            storage: Vec::new(),
            len: 0,
        };

        dna_string.extend(dna.chars().map(|c| base_to_bits(c as u8)));
        dna_string
    }

    /// Create a DnaString corresponding to an ACGT-encoded str.
    pub fn from_dna_only_string(dna: &str) -> Vec<DnaString> {
        let mut dna_vector: Vec<DnaString> = Vec::new();
        let mut dna_string = DnaString::new();

        for c in dna.chars() {
            match dna_only_base_to_bits(c as u8) {
                Some(bit) => {
                    dna_string.push(bit);
                }
                None => {
                    if !dna_string.is_empty() {
                        dna_vector.push(dna_string);
                        dna_string = DnaString::new();
                    }
                }
            }
        }
        if !dna_string.is_empty() {
            dna_vector.push(dna_string);
        }

        dna_vector
    }

    /// Create a DnaString from an ASCII ACGT-encoded byte slice.
    /// Non ACGT positions will be converted to 'A'
    pub fn from_acgt_bytes(bytes: &[u8]) -> DnaString {
        let mut dna_string = DnaString::with_capacity(bytes.len());

        // Accelerated avx2 mode. Should run on most machines made since 2013.
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                for chunk in bytes.chunks(32) {
                    if chunk.len() == 32 {
                        let (conv_chunk, _) = unsafe { crate::bitops_avx2::convert_bases(chunk) };
                        let packed = unsafe { crate::bitops_avx2::pack_32_bases(conv_chunk) };
                        dna_string.storage.push(packed);
                    } else {
                        let b = chunk.iter().map(|c| base_to_bits(*c));
                        dna_string.extend(b);
                    }
                }

                dna_string.len = bytes.len();
                return dna_string;
            }
        }

        let b = bytes.iter().map(|c| base_to_bits(*c));
        dna_string.extend(b);
        dna_string
    }

    /// Create a DnaString from an ASCII ACGT-encoded byte slice.
    /// Will return `None` if there are ambiguous bases in the DnaString
    pub fn from_acgt_bytes_checked(bytes: &[u8]) -> Result<DnaString, AmbiguousBasesError> {
        let mut dna_string = DnaString::with_capacity(bytes.len());

        // Accelerated avx2 mode. Should run on most machines made since 2013.
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                for chunk in bytes.chunks(32) {
                    if chunk.len() == 32 {
                        let (conv_chunk, correct) = unsafe { crate::bitops_avx2::convert_bases(chunk) };
                        if !correct { return Err(AmbiguousBasesError {  }) }
                        let packed = unsafe { crate::bitops_avx2::pack_32_bases(conv_chunk) };
                        dna_string.storage.push(packed);
                    } else {
                        let (b, corrects): (Vec<u8>, Vec<bool>) = chunk.iter().map(|c| base_to_bits_checked(*c)).collect();
                        let correct = !corrects.iter().contains(&false);
                        if !correct { return Err(AmbiguousBasesError {  }) }
                        dna_string.extend(b.into_iter());
                    }
                }

                dna_string.len = bytes.len();
                return Ok(dna_string);
            }
        }
        
        let (b, corrects): (Vec<u8>, Vec<bool>) = bytes.iter().map(|c| base_to_bits_checked(*c)).collect();
        let correct = corrects.iter().contains(&false);
        if !correct { return Err(AmbiguousBasesError {  }) }
        dna_string.extend(b.into_iter());

        Ok(dna_string)
    }

    /// Create a DnaString from an ACGT-encoded byte slice,
    /// Non ACGT positions will be converted to repeatable random base determined
    /// by a hash of the read name and the position within the string.
    pub fn from_acgt_bytes_hashn(bytes: &[u8], read_name: &[u8]) -> DnaString {
        let mut hasher = DefaultHasher::new();
        read_name.hash(&mut hasher);

        let mut dna_string = DnaString::with_capacity(bytes.len());

        for (pos, c) in bytes.iter().enumerate() {
            let v = match c {
                b'A' | b'a' => 0u8,
                b'C' | b'c' => 1u8,
                b'G' | b'g' => 2u8,
                b'T' | b't' => 3u8,
                _ => {
                    let mut hasher_clone = hasher.clone();
                    pos.hash(&mut hasher_clone);
                    (hasher_clone.finish() % 4) as u8
                }
            };

            dna_string.push(v);
        }

        dna_string
    }

    /// Create a DnaString from a 0-4 encoded byte slice
    pub fn from_bytes(bytes: &[u8]) -> DnaString {
        let mut dna_string = DnaString {
            storage: Vec::new(),
            len: 0,
        };

        dna_string.extend(bytes.iter().cloned());
        dna_string
    }

    /// Convert sequence to a Vector of 0-4 encoded bytes
    pub fn to_bytes(&self) -> Vec<u8> {
        self.iter().collect()
    }

    /// Convert sequence to a Vector of ascii-encoded bytes
    pub fn to_ascii_vec(&self) -> Vec<u8> {
        self.iter().map(bits_to_ascii).collect()
    }

    /// Append a 0-4 encoded base.
    #[inline]
    pub fn push(&mut self, value: u8) {
        let (block, bit) = self.addr(self.len);
        if bit == 0 && block >= self.storage.len() {
            self.storage.push(0);
        }
        self.set_by_addr(block, bit, value);
        self.len += 1;
    }

    pub fn extend(&mut self, mut bytes: impl Iterator<Item = u8>) {
        // fill the last incomplete u64 block
        while self.len % 32 != 0 {
            match bytes.next() {
                Some(b) => self.push(b),
                None => return,
            }
        }

        let mut bytes = bytes.peekable();

        // chunk the remaining items into groups of at most 32 and handle them together
        while bytes.peek().is_some() {
            let mut val: u64 = 0;
            let mut offset = 62;
            let mut n_added = 0;

            for _ in 0..32 {
                if let Some(b) = bytes.next() {
                    assert!(b < 4);
                    val |= (b as u64) << offset;
                    offset -= 2;
                    n_added += 1;
                } else {
                    break;
                }
            }

            self.storage.push(val);
            self.len += n_added;
        }
    }

    /// Push 0-4 encoded bases from a byte array.
    ///
    /// # Arguments
    /// `bytes`: byte array to read values from
    /// `seq_length`: how many values to read from the byte array. Note that this
    /// is number of values not number of elements of the byte array.
    pub fn push_bytes(&mut self, bytes: &[u8], seq_length: usize) {
        assert!(
            seq_length <= bytes.len() * 8 / WIDTH,
            "Number of elements to push exceeds array length"
        );

        for i in 0..seq_length {
            let byte_index = (i * WIDTH) / 8;
            let byte_slot = (i * WIDTH) % 8;

            let v = bytes[byte_index];
            let bits = (v >> byte_slot) & (MASK as u8);

            self.push(bits);
        }
    }

    /// Iterate over stored values (values will be unpacked into bytes).
    pub fn iter(&self) -> DnaStringIter<'_> {
        DnaStringIter {
            dna_string: self,
            i: 0,
        }
    }

    /// Clear the sequence.
    pub fn clear(&mut self) {
        self.storage.clear();
        self.len = 0;
    }

    #[inline(always)]
    fn get_by_addr(&self, block: usize, bit: usize) -> u8 {
        ((self.storage[block] >> (62 - bit)) & MASK) as u8
    }

    #[inline(always)]
    fn set_by_addr(&mut self, block: usize, bit: usize, value: u8) {
        let mask = MASK << (62 - bit);
        self.storage[block] |= mask;
        self.storage[block] ^= mask;
        self.storage[block] |= (value as u64 & MASK) << (62 - bit);
    }

    #[inline(always)]
    fn addr(&self, i: usize) -> (usize, usize) {
        let k = i * WIDTH;
        (k / BLOCK_BITS, k % BLOCK_BITS)
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Get the length `k` prefix of the DnaString
    pub fn prefix(&self, k: usize) -> DnaStringSlice<'_> {
        assert!(k <= self.len, "Prefix size exceeds number of elements.");
        DnaStringSlice {
            dna_string: self,
            start: 0,
            length: k,
            is_rc: false,
        }
    }

    /// Get the length `k` suffix of the DnaString
    pub fn suffix(&self, k: usize) -> DnaStringSlice<'_> {
        assert!(k <= self.len, "Suffix size exceeds number of elements.");

        DnaStringSlice {
            dna_string: self,
            start: self.len() - k,
            length: k,
            is_rc: false,
        }
    }

    /// Get slice containing the interval [`start`, `end`) of `self`
    pub fn slice(&self, start: usize, end: usize) -> DnaStringSlice<'_> {
        assert!(start <= self.len, "coordinate exceeds number of elements.");
        assert!(end <= self.len, "coordinate exceeds number of elements.");

        DnaStringSlice {
            dna_string: self,
            start,
            length: end - start,
            is_rc: false,
        }
    }

    /// Create a fresh DnaString containing the reverse of `self`
    pub fn reverse(&self) -> DnaString {
        let values: Vec<u8> = self.iter().collect();
        let mut dna_string = DnaString::new();
        for v in values.iter().rev() {
            dna_string.push(*v);
        }
        dna_string
    }

    /// Compute Hamming distance between this DnaString and another DnaString. The two strings must have the same length.
    pub fn hamming_distance(&self, other: &DnaString) -> usize {
        ndiffs(self, other)
    }

    // pub fn complement(&self) -> DnaString {
    //    assert!(self.width == 2, "Complement only supported for 2bit encodings.");
    //    let values: Vec<u32> = Vec::with_capacity(self.len());
    //    for i, v in self.storage.iter() {
    //        values[i] = v;
    //    }
    //    values[values.len() - 1] =
    // }
}

impl fmt::Display for DnaString {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for v in self.iter() {
            write!(f, "{}", bits_to_base(v))?;
        }
        Ok(())
    }
}

impl fmt::Debug for DnaString {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        for pos in 0..self.len() {
            s.push(bits_to_base(self.get(pos)))
        }

        write!(f, "{}", s)
    }
}

impl Default for DnaString {
    fn default() -> Self {
        Self::new()
    }
}

/// Iterator over values of a DnaStringoded sequence (values will be unpacked into bytes).
pub struct DnaStringIter<'a> {
    dna_string: &'a DnaString,
    i: usize,
}

impl Iterator for DnaStringIter<'_> {
    type Item = u8;

    fn next(&mut self) -> Option<u8> {
        if self.i < self.dna_string.len() {
            let value = self.dna_string.get(self.i);
            self.i += 1;
            Some(value)
        } else {
            None
        }
    }
}

impl<'a> IntoIterator for &'a DnaString {
    type Item = u8;
    type IntoIter = DnaStringIter<'a>;

    fn into_iter(self) -> DnaStringIter<'a> {
        self.iter()
    }
}

/// count Hamming distance between 2 2-bit DNA packed u64s
#[inline]
fn count_diff_2_bit_packed(a: u64, b: u64) -> u32 {
    let bit_diffs = a ^ b;
    let two_bit_diffs = (bit_diffs | bit_diffs >> 1) & 0x5555555555555555;
    two_bit_diffs.count_ones()
}

/// Compute the number of base positions at which two DnaStrings differ, assuming
/// that they have the same length.
pub fn ndiffs(b1: &DnaString, b2: &DnaString) -> usize {
    assert_eq!(b1.len(), b2.len());
    let mut diffs = 0;
    let (s1, s2) = (&b1.storage, &b2.storage);
    for i in 0..s1.len() {
        diffs += count_diff_2_bit_packed(s1[i], s2[i])
    }
    diffs as usize
}

/// An immutable slice into a DnaString
#[derive(Clone)]
pub struct DnaStringSlice<'a> {
    pub dna_string: &'a DnaString,
    pub start: usize,
    pub length: usize,
    pub is_rc: bool,
}

impl PartialEq for DnaStringSlice<'_> {
    fn eq(&self, other: &DnaStringSlice) -> bool {
        if other.length != self.length {
            return false;
        }
        for i in 0..self.length {
            if self.get(i) != other.get(i) {
                return false;
            }
        }
        true
    }
}
impl Eq for DnaStringSlice<'_> {}

impl<'a> Mer for DnaStringSlice<'a> {
    #[inline(always)]
    fn len(&self) -> usize {
        self.length
    }

    fn is_empty(&self) -> bool {
        self.length == 0
    }

    /// Get the base at position `i`.
    #[inline(always)]
    fn get(&self, i: usize) -> u8 {
        if !self.is_rc {
            self.dna_string.get(i + self.start)
        } else {
            crate::complement(self.dna_string.get(self.start + self.length - 1 - i))
        }
    }

    /// Set the base as position `i`.
    fn set_mut(&mut self, _: usize, _: u8) {
        unimplemented!()
        //debug_assert!(i < self.length);
        //self.dna_string.set(i + self.start, value);
    }

    fn set_slice_mut(&mut self, _: usize, _: usize, _: u64) {
        unimplemented!();
    }

    fn rc(&self) -> DnaStringSlice<'a> {
        DnaStringSlice {
            dna_string: self.dna_string,
            start: self.start,
            length: self.length,
            is_rc: !self.is_rc,
        }
    }
}

impl Vmer for DnaStringSlice<'_> {
    fn new(_: usize) -> Self {
        unimplemented!()
    }

    fn max_len() -> usize {
        <usize>::MAX
    }

    /// Get the kmer starting at position pos
    fn get_kmer<K: Kmer>(&self, pos: usize) -> K {
        debug_assert!(pos + K::k() <= self.length);
        if !self.is_rc {
            self.dna_string.get_kmer(self.start + pos)
        } else {
            let k = self
                .dna_string
                .get_kmer(self.start + self.length - K::k() - pos);
            K::rc(&k)
        }
    }
}

impl DnaStringSlice<'_> {
    pub fn is_palindrome(&self) -> bool {
        unimplemented!();
    }

    pub fn bytes(&self) -> Vec<u8> {
        let mut v = Vec::with_capacity(self.length);
        for pos in 0..self.length {
            v.push(self.get(pos));
        }
        v
    }

    pub fn ascii(&self) -> Vec<u8> {
        let mut v = Vec::with_capacity(self.length);
        for pos in 0..self.length {
            v.push(bits_to_ascii(self.get(pos)));
        }
        v
    }

    pub fn to_dna_string(&self) -> String {
        let mut dna: String = String::with_capacity(self.length);
        for pos in 0..self.length {
            dna.push(bits_to_base(self.get(pos)));
        }
        dna
    }

    pub fn to_owned(&self) -> DnaString {
        // FIXME make this faster
        let mut be = DnaString::with_capacity(self.length);
        for pos in 0..self.length {
            be.push(self.get(pos));
        }

        be
    }
    /// Get slice containing the interval [`start`, `end`) of `self`
    pub fn slice(&self, start: usize, end: usize) -> DnaStringSlice<'_> {
        assert!(
            start <= self.length,
            "coordinate exceeds number of elements."
        );
        assert!(end <= self.length, "coordinate exceeds number of elements.");
        assert!(end >= start, "invalid interval");

        if !self.is_rc {
            DnaStringSlice {
                dna_string: self.dna_string,
                start: self.start + start,
                length: end - start,
                is_rc: self.is_rc,
            }
        } else {
            // remap coords for RC
            let new_start = self.start + self.length - end;
            let new_length = end - start;

            DnaStringSlice {
                dna_string: self.dna_string,
                start: new_start,
                length: new_length,
                is_rc: self.is_rc,
            }
        }
    }

    /// Compute the Hamming distance between this DNA string and another
    pub fn hamming_dist(&self, other: &DnaStringSlice) -> u32 {
        use crate::kmer::Kmer32;
        assert_eq!(self.len(), other.len());

        let mut ndiffs = 0;

        let whole_blocks = self.len() >> 5;

        // iterate over the whole K=32 blocks
        for block in (0..whole_blocks).step_by(32) {
            let b1: Kmer32 = self.get_kmer(block);
            let b2: Kmer32 = self.get_kmer(block);
            ndiffs += count_diff_2_bit_packed(b1.to_u64(), b2.to_u64());
        }

        // iterate over trailing bases
        for pos in (whole_blocks >> 5)..self.len() {
            if self.get(pos) != other.get(pos) {
                ndiffs += 1;
            }
        }

        ndiffs
    }
}

impl fmt::Display for DnaStringSlice<'_> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for pos in 0..self.length {
            write!(f, "{}", bits_to_base(self.get(pos)))?;
        }
        Ok(())
    }
}

impl fmt::Debug for DnaStringSlice<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut s = String::new();
        if self.length < 256 {
            for pos in self.start..(self.start + self.length) {
                s.push(bits_to_base(self.dna_string.get(pos)))
            }
            write!(f, "{}", s)
        } else {
            write!(
                f,
                "start: {}, len: {}, is_rc: {}",
                self.start, self.length, self.is_rc
            )
        }
    }
}

impl<'a> IntoIterator for &'a DnaStringSlice<'a> {
    type Item = u8;
    type IntoIter = MerIter<'a, DnaStringSlice<'a>>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

/// Container for many distinct sequences, concatenated into a single DnaString.  Each
/// sequence is accessible by index as a DnaStringSlice.
#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct PackedDnaStringSet {
    pub sequence: DnaString,
    pub start: Vec<usize>,
    pub length: Vec<u32>,
}

impl PackedDnaStringSet {
    /// Create an empty `PackedDnaStringSet`
    pub fn new() -> Self {
        PackedDnaStringSet {
            sequence: DnaString::new(),
            start: Vec::new(),
            length: Vec::new(),
        }
    }

    /// Get a `DnaStringSlice` containing `i`th sequence in the set
    pub fn get(&self, i: usize) -> DnaStringSlice {
        DnaStringSlice {
            dna_string: &self.sequence,
            start: self.start[i],
            length: self.length[i] as usize,
            is_rc: false,
        }
    }

    /// Get a `DnaStringSlice` containing `i`th sequence in the set
    pub fn slice(&self, i: usize, start: usize, end: usize) -> DnaStringSlice {
        assert!(start <= self.length[i] as usize);
        assert!(end <= self.length[i] as usize);

        DnaStringSlice {
            dna_string: &self.sequence,
            start: self.start[i] + start,
            length: end - start,
            is_rc: false,
        }
    }

    /// Number of sequences in the set
    pub fn len(&self) -> usize {
        self.start.len()
    }

    pub fn is_empty(&self) -> bool {
        self.start.is_empty()
    }

    pub fn add<R: Borrow<u8>, S: IntoIterator<Item = R>>(&mut self, sequence: S) {
        let start = self.sequence.len();
        self.start.push(start);

        let mut length = 0;
        for b in sequence {
            self.sequence.push(*b.borrow());
            length += 1;
        }
        self.length.push(length as u32);
        //debug!("add to sequence for loop {:?} iterations (pr seq len)", length);
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AmbiguousBasesError {}

impl Error for AmbiguousBasesError {}

impl Display for AmbiguousBasesError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ambiguous base found in read")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::IntKmer;
    use crate::kmer::Kmer4;
    use crate::test;
    use rand::{self, Rng};

    fn hamming_dist_slow(s1: &DnaStringSlice, s2: &DnaStringSlice) -> u32 {
        assert_eq!(s1.len(), s2.len());

        let mut ndiff = 0;
        for pos in 0..s1.len() {
            if s1.get(pos) != s2.get(pos) {
                ndiff += 1;
            }
        }

        ndiff
    }

    /// Randomly mutate each base with probability `p`
    pub fn edit_dna_string(string: &mut DnaString, p: f64, r: &mut impl Rng) {
        for pos in 0..string.len() {
            if r.gen_range(0.0, 1.0) < p {
                let base = test::random_base(r);
                string.set_mut(pos, base)
            }
        }
    }

    fn random_dna_string(len: usize) -> DnaString {
        let bytes = test::random_dna(len);
        DnaString::from_bytes(&bytes)
    }

    fn random_dna_string_pair(len: usize) -> (DnaString, DnaString) {
        let s1 = random_dna_string(len);
        let mut s2 = s1.clone();

        if rand::thread_rng().gen_bool(0.1) {
            s2 = s2.rc();
        }

        edit_dna_string(&mut s2, 0.1, &mut rand::thread_rng());
        (s1, s2)
    }

    fn random_slice(len: usize) -> (usize, usize) {
        if len == 0 {
            return (0, 0);
        }
        let mut r = rand::thread_rng();
        let a = (rand::RngCore::next_u64(&mut r) as usize) % len;
        let b = (rand::RngCore::next_u64(&mut r) as usize) % len;
        (std::cmp::min(a, b), std::cmp::max(a, b))
    }

    #[test]
    fn dnastringslice_get_kmer() {
        let seq = DnaString::from_dna_string("ACGGTAC");
        let seqrc = DnaString::from_dna_string("GTACCGT");
        let rcslice = seq.slice(0, 7).rc();
        let slice = seqrc.slice(0, 7);
        for i in 0..=3 {
            // The kmer in a slice should be the kmer of the sequencing represented
            // by that slice, regardless of whether the backing DnaString is RC or not.
            assert_eq!(slice.get_kmer::<Kmer4>(i), rcslice.get_kmer::<Kmer4>(i));
        }
    }
    #[test]
    fn dnastringslice_slice() {
        let seq = DnaString::from_dna_string("ACGGTAC");
        let seqrc = DnaString::from_dna_string("GTACCGT");
        let rcslice = seq.slice(0, 7).rc();
        let slice = seqrc.slice(0, 7);
        // The fact that a slice is backed by a DnaStringSlice that is
        // the rc of the slice sequence shouldn't matter.
        assert_eq!(rcslice.slice(1, 4), slice.slice(1, 4));
    }

    #[test]
    fn test_slice_hamming_dist() {
        for len in 0..1000 {
            for _ in 0..5 {
                let (s1, s2) = random_dna_string_pair(len);
                let (start, end) = random_slice(len);
                let slc1 = s1.slice(start, end);
                let slc2 = s2.slice(start, end);

                let validation_dist = hamming_dist_slow(&slc1, &slc2);
                let test_dist = slc1.hamming_dist(&slc2);
                assert_eq!(validation_dist, test_dist);
            }
        }
    }

    #[test]
    fn test_dna_string() {
        let mut dna_string = DnaString::with_capacity(1000);
        assert_eq!(dna_string.len(), 0);
        dna_string.push(0);
        dna_string.push(2);
        dna_string.push(1);
        println!("{:?}", dna_string);
        let mut values: Vec<u8> = dna_string.iter().collect();
        assert_eq!(values, [0, 2, 1]);
        dna_string.set_mut(1, 3);
        values = dna_string.iter().collect();
        assert_eq!(values, [0, 3, 1]);
    }

    #[test]
    fn test_push_bytes() {
        let in_values: Vec<u8> = vec![2, 20];

        let mut dna_string = DnaString::new();
        dna_string.push_bytes(&in_values, 8);
        // Contents should be [01000000, 00101000]
        let values: Vec<u8> = dna_string.iter().collect();
        assert_eq!(values, [2, 0, 0, 0, 0, 1, 1, 0]);

        let mut dna_string = DnaString::new();
        dna_string.push_bytes(&in_values, 2);
        // Contents should be 01000000
        let values: Vec<u8> = dna_string.iter().collect();
        assert_eq!(values, [2, 0]);
    }

    #[test]
    fn test_from_dna_string() {
        dna_string_test("");
        dna_string_test("A");
        dna_string_test("C");
        dna_string_test("G");
        dna_string_test("T");

        dna_string_test("GC");
        dna_string_test("ATA");

        dna_string_test("ACGTACGT");
        dna_string_test("ACGTAAAAAAAAAATTATATAACGT");
        dna_string_test("AACGTAAAAAAAAAATTATATAACGT");
    }

    fn dna_string_test(dna: &str) {
        let dna_string_a = DnaString::from_dna_string(dna);
        let dna_string = DnaString::from_acgt_bytes(dna.as_bytes());
        assert_eq!(dna_string_a, dna_string);

        let rc = dna_string_a.rc();
        let rc2 = rc.rc();
        assert_eq!(dna_string_a, rc2);

        assert_eq!(dna_string.iter().count(), dna.len());

        assert_eq!(dna_string.len, dna.len());

        let dna_cp = dna_string.to_string();
        assert_eq!(dna, dna_cp);
    }

    #[test]
    fn test_dna_string_ambig() {
        let dna = [
            "TTTTTTTTTTTTTTTTTTTTTTTT",
            "NAGCGGAGATTATTCACGAGCATCGCGTAC",
            "GATCGATGCATGCTAGN",
            "ACGTAAAAAAAAAATTATATAACGTACGTAAAAAAAAAATTATATAACGTAACGTAAAAANAAAAATTATANTAACGT",
            "AGCTAGCTAGCTGACTGAGCGACTGA",
            "AGCTAGCTAGCTGACTGAGCGACTGACGGATC",
            "GCATCGAGCATGCTACGATGCGACGATCGTACGATCGTACGATC",
            "ACGATCGNATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACTGATCGATCG",
            "ACGATCGATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACTGATCGATCGJHSJDSDHKAJSHDK",
        ];

        let comp_unchecked = [
            "TTTTTTTTTTTTTTTTTTTTTTTT",
            "AAGCGGAGATTATTCACGAGCATCGCGTAC",
            "GATCGATGCATGCTAGA",
            "ACGTAAAAAAAAAATTATATAACGTACGTAAAAAAAAAATTATATAACGTAACGTAAAAAAAAAAATTATAATAACGT",
            "AGCTAGCTAGCTGACTGAGCGACTGA",
            "AGCTAGCTAGCTGACTGAGCGACTGACGGATC",
            "GCATCGAGCATGCTACGATGCGACGATCGTACGATCGTACGATC",
            "ACGATCGAATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACTGATCGATCG",
            "ACGATCGATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACTGATCGATCGAAAAAAAAAAAAAAA",
        ];

        let comp_checked = [
            Ok("TTTTTTTTTTTTTTTTTTTTTTTT".to_string()),
            Err(AmbiguousBasesError {}),
            Err(AmbiguousBasesError {}),
            Err(AmbiguousBasesError {}),
            Ok("AGCTAGCTAGCTGACTGAGCGACTGA".to_string()),
            Ok("AGCTAGCTAGCTGACTGAGCGACTGACGGATC".to_string()),
            Ok("GCATCGAGCATGCTACGATGCGACGATCGTACGATCGTACGATC".to_string()),
            Err(AmbiguousBasesError {}),
            Err(AmbiguousBasesError {}),
        ];

        for (i, seq) in dna.iter().enumerate() {
            assert_eq!(format!("{:?}", DnaString::from_acgt_bytes(seq.as_bytes())), comp_unchecked[i]);
            assert_eq!(DnaString::from_acgt_bytes_checked(seq.as_bytes()).map(|d| format!("{:?}", d)), comp_checked[i]);
        }
    }

    #[test]
    fn test_prefix() {
        let in_values: Vec<u8> = vec![2, 20];
        let mut dna_string = DnaString::new();
        dna_string.push_bytes(&in_values, 8);
        // Contents should be [01000000, 00101000]

        let pref_dna_string = dna_string.prefix(0).to_owned();
        assert_eq!(pref_dna_string.len(), 0);

        let pref_dna_string = dna_string.prefix(8).to_owned();
        assert_eq!(pref_dna_string, dna_string);

        let pref_dna_string = dna_string.prefix(4).to_owned();
        let values: Vec<u8> = pref_dna_string.iter().collect();
        assert_eq!(values, [2, 0, 0, 0]);

        let pref_dna_string = dna_string.prefix(6).to_owned();
        let values: Vec<u8> = pref_dna_string.iter().collect();
        assert_eq!(values, [2, 0, 0, 0, 0, 1]);

        dna_string.push_bytes(&in_values, 8);
        dna_string.push_bytes(&in_values, 8);

        let pref_dna_string = dna_string.prefix(17).to_owned();
        let values: Vec<u8> = pref_dna_string.iter().collect();
        assert_eq!(values, [2, 0, 0, 0, 0, 1, 1, 0, 2, 0, 0, 0, 0, 1, 1, 0, 2]);
    }

    #[test]
    fn test_suffix() {
        let in_values: Vec<u8> = vec![2, 20];
        let mut dna_string = DnaString::new();
        dna_string.push_bytes(&in_values, 8);
        // Contents should be [01000000, 00101000]

        let suf_dna_string = dna_string.suffix(0).to_owned();
        assert_eq!(suf_dna_string.len(), 0);

        let suf_dna_string = dna_string.suffix(8).to_owned();
        assert_eq!(suf_dna_string, dna_string);

        let suf_dna_string = dna_string.suffix(4).to_owned();
        let values: Vec<u8> = suf_dna_string.iter().collect();
        assert_eq!(values, [0, 1, 1, 0]);

        // 000101000000 64+256
        let suf_dna_string = dna_string.suffix(6).to_owned();
        let values: Vec<u8> = suf_dna_string.iter().collect();
        assert_eq!(values, [0, 0, 0, 1, 1, 0]);

        dna_string.push_bytes(&in_values, 8);
        dna_string.push_bytes(&in_values, 8);

        let suf_dna_string = dna_string.suffix(17).to_owned();
        let values: Vec<u8> = suf_dna_string.iter().collect();
        assert_eq!(values, [0, 2, 0, 0, 0, 0, 1, 1, 0, 2, 0, 0, 0, 0, 1, 1, 0]);
    }

    #[test]
    fn test_reverse() {
        let in_values: Vec<u8> = vec![2, 20];
        let mut dna_string = DnaString::new();
        let rev_dna_string = dna_string.reverse();
        assert_eq!(dna_string, rev_dna_string);

        dna_string.push_bytes(&in_values, 8);
        // Contents should be 00010100 00000010

        let rev_dna_string = dna_string.reverse();
        let values: Vec<u8> = rev_dna_string.iter().collect();
        assert_eq!(values, [0, 1, 1, 0, 0, 0, 0, 2]);
    }

    #[test]
    fn test_kmers() {
        const DNA: &str = "TGCATTAGAAAACTCCTTGCCTGTCAGCCCGACAGGTAGAAACTCATTAATCCACACATTGA\
            CTCTATTTCAGGTAAATATGACGTCAACTCCTGCATGTTGAAGGCAGTGAGTGGCTGAAACAGCATCAAGGCGTGAAGGC";
        let dna_string = DnaString::from_dna_string(DNA);

        let kmers: Vec<IntKmer<u64>> = dna_string.iter_kmers().collect();
        kmer_test::<IntKmer<u64>>(&kmers, DNA, &dna_string);
    }

    #[test]
    fn test_ndiffs() {
        let x1 = DnaString::from_dna_string("TGCATTAGAAAACTCCTTGCCTGTCTAGAAACTCATTAATCCACACATTGA");
        let x2 = DnaString::from_dna_string("TGCATTAGTAAACTCCTTCGCTGTCTAGAAAATCATTAAGCCACACATTGA");
        assert_eq!(ndiffs(&x1, &x2), 5);

        let x1 = DnaString::from_dna_string("TGCATT");
        let x2 = DnaString::from_dna_string("TGCATT");
        assert_eq!(ndiffs(&x1, &x2), 0);

        let x1 = DnaString::from_dna_string("");
        let x2 = DnaString::from_dna_string("");
        assert_eq!(ndiffs(&x1, &x2), 0);

        let x1 = DnaString::from_dna_string("TGCATTAGAAAACTCCTTGCCTGTCTAGAAACTCATTAATCCACACATTGA\
            TGCATTAGAAAACTCCTTGCCTGTCTAGAAACTCATTAATCCACACATTGATGCATTAGAAAACTCCTTGCCTGTCTAGAAACTCATTAATCCACACATTGA");
        let x2 = DnaString::from_dna_string("TGCATTAGTAAACTCCTTCGCTGTCTAGAAAATCATTAAGCCACACATTGA\
            TGCATTAGTAAACTCCTTCGCTGTCTAGAAAATCATTAAGCCACACATTGATGCATTAGTAAACTCCTTCGCTGTCTAGAAAATCATTAAGCCACACATTGA");
        assert_eq!(ndiffs(&x1, &x2), 15);
    }

    #[test]
    fn test_kmers_too_short() {
        const DNA: &str = "TGCATTAGAA";
        let dna_string = DnaString::from_dna_string(DNA);

        let kmers: Vec<IntKmer<u64>> = dna_string.iter_kmers().collect();
        assert_eq!(kmers, Vec::default());
    }

    fn kmer_test<K: Kmer>(kmers: &[K], dna: &str, dna_string: &DnaString) {
        for i in 0..(dna.len() - K::k() + 1) {
            assert_eq!(kmers[i].to_string(), &dna[i..(i + K::k())]);
        }

        let last_kmer: K = dna_string.last_kmer();
        assert_eq!(last_kmer.to_string(), &dna[(dna.len() - K::k())..]);

        for (idx, &k) in kmers.iter().enumerate() {
            assert_eq!(k, dna_string.get_kmer(idx));
        }
    }
}
