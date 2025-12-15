use std::collections::HashMap;
use std::mem::take;
use std::ops::Range;
use bimap::BiMap;
use itertools::Itertools;
use serde::de::DeserializeOwned;
use serde_derive::{Deserialize, Serialize};
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::{mem, str};
use crate::dna_string::DnaString;
use crate::summarizer::{IDTag, Tag, ID};
use crate::{base_to_bits, base_to_bits_checked, Exts, Vmer};

#[derive(Debug, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Clone, Hash, Copy)]
pub enum Strandedness {
    Forward,
    Reverse,
    Unstranded
}

#[derive(Debug, Clone, Copy)]
pub enum ReadEnd {
    R1,
    R2
}

/// a sequencing read and additional information
#[derive(Debug, PartialEq, )]
pub struct Read<D> {
    seq: DnaString,
    exts: Exts,
    data: D,
    strand: Strandedness,
}

impl<D: Clone + Copy> Read<D> {
    /// a new `Read`
    pub fn new(seq: DnaString, exts: Exts, data: D, strand: Strandedness) -> Read<D> {
        Read { seq, exts, data, strand }
    }

    /// the sequence of the `Read`
    pub fn seq(&self) -> &DnaString {
        &self.seq
    }

    /// the [`Exts`] of the `Read`
    pub fn exts(&self) -> Exts {
        self.exts
    }

    /// the data of the `Read`
    pub fn data(&self) -> D {
        self.data
    }

    /// the strandedness of the `Read`
    pub fn stranded(&self) -> Strandedness {
        self.strand
    }
}

/// two paired sequencing reads
pub struct PairedRead<D> {
    read1: Read<D>,
    read2: Read<D>
}

impl<D> PairedRead<D> {
    /// a new `PairedRead`
    pub fn new(read1: Read<D>, read2: Read<D>) -> PairedRead<D> {
        PairedRead {
            read1,
            read2
        }
    }

    /// the paired reads
    pub fn reads(&self) -> (&Read<D>, &Read<D>) {
        (&self.read1, &self.read2)
    }
}

/// Store many DNA sequences together with an Exts and data each compactly packed together
/// 
/// #### fields:
/// 
/// * `storage`: `Vec` with 2-bit encoded DNA bases of all sequences
/// * `ends`:  `Vec` with the ends (exclusive) of the separate sequences in the `Reads`
/// * `exts`: `Option<Vec>` with one Exts for each sequence
/// * `data`: `Vec` with data for each sequence
/// * `len`: length of all sequences together
/// * `stranded`: [`Stranded`] conveying the strandedness and direction of the reads
#[derive(Ord, PartialOrd, Clone, PartialEq, Eq, Hash, Serialize, Deserialize, Debug)]
pub struct Reads<D> {
    storage: Vec<u64>,
    ends: Vec<usize>,
    exts: Option<Vec<Exts>>,
    data: Vec<D>,
    len: usize,
    stranded: Strandedness
}

impl<D: Clone + Copy> Reads<D> {

    /// Returns a new `Reads`
    pub fn new(stranded: Strandedness) -> Self {
        Reads {
            storage: Vec::new(),
            ends: Vec::new(),
            exts: None,
            data: Vec::new(),
            len: 0,
            stranded
        }
    }

    /// Returns a new `Reads` with exts
    pub fn new_with_exts(stranded: Strandedness) -> Self {
        Reads {
            storage: Vec::new(),
            ends: Vec::new(),
            exts: Some(Vec::new()),
            data: Vec::new(),
            len: 0,
            stranded
        }
    }

    #[inline(always)]
    pub fn stranded(&self) -> Strandedness {
        self.stranded
    }

    #[inline(always)]
    /// Returns the number of reads stored
    pub fn n_reads(&self) -> usize {
        self.ends.len()
    }

    /// get the memory required for the reads
    pub fn mem(&self) -> usize {
        let exts_size = if let Some(e_vec) = self.exts.as_ref() { size_of_val(&**e_vec) } else { 0 };
        mem::size_of_val(self) + size_of_val(&*self.storage) + size_of_val(&*self.data) + size_of_val(&*self.ends) + exts_size
    }

    /// set the strandedness and the direction of the reads
    pub fn set_stranded(&mut self, stranded: Strandedness) {
        self.stranded = stranded
    }

    /// add exts to `Reads` if needed - use after adding sequence and new end
    fn add_exts(&mut self, exts: Option<Exts>) {
        match exts {
            Some(e) => {
                // only add exts if some exts are not empty
                match self.exts.as_mut() {
                    // already has exts, simply append new exts
                    Some(e_vec) => e_vec.push(e),
                    // no exts so far
                    None => {
                        // check if exts are empty
                        if e != Exts::empty() {
                            // if not, add vector of empty exts and then push new exts
                            self.exts = Some(vec![Exts::empty(); self.n_reads() - 1]);
                            self.exts.as_mut().unwrap().push(e);
                        } // else keep no exts
                    }
                }
            }
            None => {
                if let Some(e_vec) = self.exts.as_mut() {
                    // Reads has exts, but no exts here, so add empty exts
                    e_vec.push(Exts::empty())
                } // else do nothing
            }
        }
        
    }

    /// Adds a new read to the `Reads`
    // maybe push_base until u64 is full and then do extend like in DnaString::extend ? with accellerated mode
    pub fn add_read<V: Vmer>(&mut self, seq: V, exts: Option<Exts>, data: D) {
        for base in seq.iter() {
            self.push_base(base);
        }
        self.ends.push(self.len);
        self.add_exts(exts);
        self.data.push(data);
    }

    /// Transforms a `[(vmer, exts, data)]` into a `Reads` - watch for memory usage
    // TODO test if memory efficient
    pub fn from_vmer_vec<V: Vmer, S: IntoIterator<Item=(V, Exts, D)>>(vec_iter: S, stranded: Strandedness) -> Self {
        let mut reads = Reads::new(stranded);
        for (vmer, exts, data) in vec_iter {
            for base in vmer.iter() {
                reads.push_base(base);
            }
            reads.ends.push(reads.len);
            reads.add_exts(Some(exts));
            reads.data.push(data);
        }

        reads.shrink_to_fit();
        
        reads
    }


    /// add ASCII encoded bases to the `Reads`
    /// 
    /// will transform all ascii characters outside of ACGTacgt into A
    /// 
    /// if `Reads` previously contained no exts, no new exts will be added
    /// see also: [`Reads::add_from_bytes_checked`]
    pub fn add_from_bytes(&mut self, bytes: &[u8], exts: Option<Exts>, data: D) {
        
        // fill the last incomplete u64 block
        let missing = 32 - (self.len % 32);
        if missing != 0 {
            if  missing > bytes.len() {
                let fill = bytes.iter().map(|c| base_to_bits(*c));
                self.extend(fill);
                self.ends.push(self.len);
                self.add_exts(exts);
                self.data.push(data);
                return;
            } else {
                let fill = bytes[0..missing].iter().map(|c| base_to_bits(*c));
                self.extend(fill);
            }
        }
        
        // Accelerated avx2 mode. Should run on most machines made since 2013.
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                for chunk in bytes[missing..bytes.len()].chunks(32) {
                    if chunk.len() == 32 {
                        let (conv_chunk, _) = unsafe { crate::bitops_avx2::convert_bases(chunk) };
                        let packed = unsafe { crate::bitops_avx2::pack_32_bases(conv_chunk) };
                        self.storage.push(packed);
                        self.len += 32;
                    } else {
                        let b = chunk.iter().map(|c| base_to_bits(*c));
                        self.extend(b);
                    }
                }
                self.ends.push(self.len);
                self.add_exts(exts);
                self.data.push(data);
                
                return;
            }
        }

        let b = bytes.iter().map(|c| base_to_bits(*c));
        self.extend(b);
        self.ends.push(self.len);
        self.add_exts(exts);
        self.data.push(data);
        
        
    }

    /// add ASCII encoded bases to the Reads
    /// 
    /// will return `false` if the bytes contained characters outside of `ACGTacgt`, otherwise return true and add the bases
    /// see also: [`Reads::add_from_bytes`]
    pub fn add_from_bytes_checked(&mut self, bytes: &[u8], exts: Option<Exts>, data: D) -> bool {

        let (_, corrects): (Vec<u8>, Vec<bool>) = bytes.iter().map(|c| base_to_bits_checked(*c)).collect();
        if corrects.iter().contains(&false) { return false }

        
        // fill the last incomplete u64 block
        let missing = 32 - (self.len % 32);
        if missing != 0 {
            if  missing > bytes.len() {
                let fill = bytes.iter().map(|c| base_to_bits(*c));
                self.extend(fill);
                self.add_exts(exts);
                self.data.push(data);
                self.ends.push(self.len);
                return true;
            } else {
                let fill = bytes[0..missing].iter().map(|c| base_to_bits(*c));
                self.extend(fill);
            }
        }
        
        // Accelerated avx2 mode. Should run on most machines made since 2013.
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                for chunk in bytes[missing..bytes.len()].chunks(32) {
                    if chunk.len() == 32 {
                        let (conv_chunk, _) = unsafe { crate::bitops_avx2::convert_bases(chunk) };
                        let packed = unsafe { crate::bitops_avx2::pack_32_bases(conv_chunk) };
                        self.storage.push(packed);
                        self.len += 32;
                    } else {
                        let b = chunk.iter().map(|c| base_to_bits(*c));
                        self.extend(b);
                    }
                }

                self.ends.push(self.len);
                self.add_exts(exts);
                self.data.push(data);
                
                return true;
            }
        }

        let b = bytes.iter().map(|c| base_to_bits(*c));
        self.extend(b);
        self.ends.push(self.len);
        self.add_exts(exts);
        self.data.push(data);

        true         
    }


    /// Add new 2-bit encoded base to the `Reads`
    fn push_base(&mut self, base: u8) {
        let bit = (self.len % 32) * 2;
        if bit != 0 {
            match self.storage.pop() {
                Some(last) => {
                    let last = last + ((base as u64) << (64 - bit - 2));
                    self.storage.push(last);
                },
                None => panic!("tried to push base to empty vector (?)")
            }
        } else {
            self.storage.push((base as u64) << 62);
        }
        self.len += 1; 
    }


    /// extend the reads' storage by 2-bit encoded bases
    fn extend(&mut self, mut bytes: impl Iterator<Item = u8>) {
        // fill the last incomplete u64 block
        while self.len % 32 != 0 {
            match bytes.next() {
                Some(b) => self.push_base(b),
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

    #[inline(always)]
    fn addr(&self, i: &usize) -> (usize, usize) {
        (i / 32, (i % 32 ) * 2)
    }

    /// get the `i`th read in a `Reads`
    pub fn get_read(&self, i: usize) -> Option<Read<D>> {
        if i >= self.n_reads() { return None }

        let mut sequence = DnaString::new();
        let end = self.ends[i];
        //let start = if i != 0 { self.ends[i-1] } else { 0 };
        let start = match i {
            0 => 0,
            1.. => self.ends[i-1]
        };

        for b in start..end {
            let (block, bit) = self.addr(&b);
            let base = ((self.storage[block] >> (62 - bit)) & 3u64) as u8;
            sequence.push(base);
        }

        let exts = match self.exts {
            Some(ref e_vec) => e_vec[i],
            None => Exts::empty()
        };

        Some(Read::new(sequence, exts, self.data[i], self.stranded))
    }


    /// shrink the vectors' capacity to fit the length
    /// 
    /// use sparsely
    pub fn shrink_to_fit(&mut self)  {
        self.storage.shrink_to_fit();
        self.data.shrink_to_fit();
        if let Some(e_vec) = self.exts.as_mut() { e_vec.shrink_to_fit(); }
        self.ends.shrink_to_fit();
    }

    /// Iterate over the reads as (DnaString, Exts, D).
    pub fn iter(&self) -> ReadsIter<'_, D> {
        ReadsIter {
            reads: self,
            i: 0,
            end: self.n_reads(),
            length: self.n_reads(),
        }
    }

    /// Iterate over a range start reads as (DnaString, Exts, D).
    pub fn partial_iter(&self, range: Range<usize>) -> ReadsIter<'_, D> {
        assert!(range.end <= self.n_reads());
        assert!(range.start < self.n_reads());
        assert!(range.start < range.end);
        ReadsIter {
            reads: self,
            i: range.start,
            end: range.end,
            length: (range.end - range.start)
        }
    }

    pub fn info(&self) -> String {
        format!("Reads {{ n reads: {}, stranded: {:?} }}", self.n_reads(), self.stranded)
    }
}

impl<D: ReadData> Reads<D> {
    /// get the number of k-mers for each unique data value
    pub fn tag_kmers(&self, k: usize) -> HashMap<Tag, usize> {
        let mut hm = HashMap::new();

        self.iter().for_each(|read| {
            let kmers = read.seq.len().saturating_sub(k - 1);
            if let Some(tag) = read.data.get_tag() {
                if let Some(count) = hm.get_mut(&tag) {
                    *count += kmers;
                } else {
                    hm.insert(tag, kmers);
                }
            }
            
        });

        hm
    }
}

impl<D: Clone + Copy> Default for Reads<D> {
    fn default() -> Self {
        Self::new(Strandedness::Unstranded)
    }
}

/// Iterator over values of a DnaStringoded sequence (values will be unpacked into bytes).
pub struct ReadsIter<'a, D> {
    reads: &'a Reads<D>,
    i: usize,
    end: usize,
    length: usize,
}

impl<D: Clone + Copy> Iterator for ReadsIter<'_, D> {
    type Item = Read<D>;

    fn next(&mut self) -> Option<Self::Item> {
        if (self.i < self.reads.n_reads()) && (self.i < self.end) {
            let value = self.reads.get_read(self.i);
            self.i += 1;
            value
        } else {
            None
        }
    }
}

impl<D: Copy> ExactSizeIterator for ReadsIter<'_, D> {
    fn len(&self) -> usize {
        self.length
    }
}

impl<D: Clone + Copy + Debug> Display for Reads<D> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let vec: Vec<_> = self.iter().collect();
        write!(f, "{:?}", vec)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ReadsPaired<D> {
    Empty,
    Unpaired { reads: Reads<D> },
    Paired { paired1: Reads<D>, paired2: Reads<D> },
    Combined {paired1: Reads<D>, paired2: Reads<D>, unpaired: Reads<D>}
}

impl<D: Clone + Copy> ReadsPaired<D> {
    /// return an iterable element (`Vec`) with all contained `Reads`
    pub fn iterable(&self) -> Vec<&Reads<D>> {
        match self {
            Self::Empty => vec![],
            Self::Unpaired { reads  } => vec![reads],
            Self::Paired { paired1, paired2 } => vec![paired1, paired2],
            Self::Combined { paired1, paired2, unpaired } => vec![paired1, paired2, unpaired],
        }
    }

    /// the overall number of reads
    pub fn n_reads(&self) -> usize {
        match self {
            Self::Empty => 0,
            Self::Unpaired { reads  } => reads.n_reads(),
            Self::Paired { paired1, paired2 } => paired1.n_reads() + paired2.n_reads(),
            Self::Combined { paired1, paired2, unpaired } => paired1.n_reads() + paired2.n_reads() + unpaired.n_reads(),
        }
    }

    /// the number of paired reads
    pub fn n_read_pairs(&self) -> usize {
        match self {
            Self::Empty => 0,
            Self::Unpaired { reads: _ } => 0,
            Self::Paired { paired1, paired2 } => {
                assert_eq!(paired1.n_reads(), paired2.n_reads());
                paired1.n_reads()
            },
            Self::Combined { paired1, paired2, unpaired: _ } => {
                assert_eq!(paired1.n_reads(), paired2.n_reads());
                paired1.n_reads()
            }
        }
    }

    /// the number of unpaired reads
    pub fn n_unpaired_reads(&self) -> usize {
        match self {
            Self::Empty => 0,
            Self::Unpaired { reads } => reads.n_reads(),
            Self::Paired { paired1: _, paired2: _ } => 0,
            Self::Combined { paired1: _, paired2: _, unpaired } => {
                unpaired.n_reads()
            }
        }
    }

    /// get the read with the index `i` from the `ReadsPaired` - with multiple 
    /// underlying `Reads` its is counted linearly trough paired1, paired2, unpaired
    pub fn get_read(&self, i: usize) -> Option<Read<D>> {
        match self {
            ReadsPaired::Empty => None,
            ReadsPaired::Unpaired { reads } => reads.get_read(i),
            ReadsPaired::Paired { paired1, paired2 } => {
                if i < paired1.n_reads() {
                    paired1.get_read(i)
                } else if (i - paired1.n_reads()) < paired2.n_reads() {
                    paired2.get_read(i -  paired1.n_reads())
                } else {
                    None
                }
            },
            ReadsPaired::Combined { paired1, paired2, unpaired } => {
                if i < paired1.n_reads() {
                    paired1.get_read(i)
                } else if (i - paired1.n_reads()) < paired2.n_reads() {
                    paired2.get_read(i -  paired1.n_reads())
                } else if (i - (paired1.n_reads() + paired2.n_reads())) < unpaired.n_reads() {
                    unpaired.get_read(i - (paired1.n_reads() + paired2.n_reads()))
                } else {
                    None
                }
            },
        }
    }

    pub fn get_paired_read(&self, i: usize) -> Option<PairedRead<D>> {
        match self {
            ReadsPaired::Empty => None,
            ReadsPaired::Unpaired { reads: _ } => None,
            ReadsPaired::Paired { paired1, paired2 } => {
                if let (Some(read1), Some(read2)) = (paired1.get_read(i), paired2.get_read(i)) {
                    Some(PairedRead::new(read1, read2))
                } else {
                    None
                }
            }
            ReadsPaired::Combined { paired1, paired2, unpaired: _ } => {
                if let (Some(read1), Some(read2)) = (paired1.get_read(i), paired2.get_read(i)) {
                    Some(PairedRead::new(read1, read2))
                } else {
                    None
                }
            },
        }
    }

    pub fn mem(&self) -> usize {
        match self {
            Self::Empty => 0,
            Self::Unpaired { reads } => reads.mem(),
            Self::Paired { paired1, paired2 } => paired1.mem() + paired2.mem(),
            Self::Combined { paired1, paired2, unpaired } => paired1.mem() + paired2.mem() + unpaired.mem(),
        }
    }

    /// transform a tuple of two paired [`Reads`] and one unpaired [`Reads`] into a `ReadsPaired`
    /// depending on the contents of the [`Reads`]
    pub fn from_reads((paired1, paired2, unpaired): (Reads<D>, Reads<D>, Reads<D>)) -> Self {
        // first two elements should be paired reads and thus have same n
        assert_eq!(paired1.n_reads(), paired2.n_reads(), "Error: R1 read and R2 read counts have to match");

        if (paired1.n_reads() + paired2.n_reads() + unpaired.n_reads()) == 0 {
            // no reads
            ReadsPaired::Empty
        } else if paired1.n_reads() == 0 && unpaired.n_reads() > 0 {
            // only reads in third element -> unpaired
            ReadsPaired::Unpaired { reads: unpaired }
        } else if paired1.n_reads() > 0 && unpaired.n_reads() == 0 {
            // reads in first and second element -> paired
            ReadsPaired::Paired { paired1, paired2 }
        } else if paired1.n_reads() > 0 && unpaired.n_reads() > 0 {
            // reads in all elements: both paired and unpaired reads
            ReadsPaired::Combined { paired1, paired2, unpaired }
        } else {
            panic!("error in transforming Reads into ReadsPaired")
        }
    }

    pub fn iter(&self) -> Box<dyn Iterator<Item = Read<D>> + '_> {
        match self {
            ReadsPaired::Empty => panic!("Error: no reads to process"),
            ReadsPaired::Unpaired { reads } => Box::new(reads.iter()),
            ReadsPaired::Paired { paired1, paired2 } => Box::new(paired1.iter().chain(paired2.iter())),
            ReadsPaired::Combined { paired1, paired2, unpaired } => Box::new(paired1.iter().chain(paired2.iter()).chain(unpaired.iter())),
        }
    }

    pub fn iter_partial(&self, range: Range<usize>) -> Box<dyn Iterator<Item = Read<D>> + '_> {
        match self {
            Self::Empty => panic!("Error: no reads to process"),
            Self::Unpaired { reads } => Box::new(reads.partial_iter(range)),
            Self::Paired { paired1, paired2 } => {
                let n_p1 = paired1.n_reads();
                if range.start >= n_p1 {
                    // range is fully in paired2
                    Box::new(paired2.partial_iter((range.start - n_p1)..(range.end - n_p1)))
                } else if range.end <= n_p1 {
                    // range is fully in paired1
                    Box::new(paired1.partial_iter(range))
                } else {
                    // range is both in paired1 and paired2
                    Box::new(paired1.partial_iter(range.start..n_p1).chain(paired2.partial_iter(0..(range.end - n_p1))))
                }
            },
            Self::Combined { paired1, paired2, unpaired } => {
                let n_p1 = paired1.n_reads();
                let n_p2 = paired2.n_reads();
                let n_p12 = n_p1 + paired2.n_reads();
                if range.end <= n_p1 {
                    // range is only in paired1
                    Box::new(paired1.partial_iter(range))
                } else if range.end >= n_p1 && range.end <= n_p12 && range.start >= n_p1 && range.start <= n_p12 {
                    // range is only in paired2
                    Box::new(paired2.partial_iter((range.start - n_p1)..(range.end - n_p1)))
                } else if range.start >= n_p12 {
                    // range is only in unpaired
                    Box::new(unpaired.partial_iter((range.start - n_p12)..(range.end - n_p12)))
                } else if range.start <= n_p1 && range.end >= n_p1 && range.end <= n_p12 {
                    // range is in paired1 and paired2
                    Box::new(paired1.partial_iter(range.start..n_p1).chain(paired2.partial_iter(0..(range.end - n_p1))))
                } else if range.start >= n_p1 && range.start <= n_p12 && range.end >= n_p12 {
                    // range is in paired2 and unpaired
                    Box::new(paired2.partial_iter((range.start - n_p1)..n_p2).chain(unpaired.partial_iter(0..(range.end - n_p12))))
                } else {
                    // range is in paired1, paired2, and in unpaired
                    Box::new(paired1.partial_iter(range.start..n_p1).chain(paired2.partial_iter(0..n_p2)).chain(unpaired.partial_iter(0..(range.end - n_p12))))
                }
            }
        }
    }

    /// if the `ReadsPaired` is of `Combined` type, remove the unpaired reads,
    /// returns the number of reads that were removed
    pub fn decombine(&mut self) -> usize {
        if let Self::Combined { paired1, paired2, unpaired } = self {
            let rm_reads = unpaired.n_reads();
            *self = ReadsPaired::Paired { paired1: take(paired1), paired2: take(paired2) };
            rm_reads
        } else {
            0
        }
    }
}

impl<DI: ReadData> ReadsPaired<DI> {
    /// get the number of k-mers for each unique data value
    pub fn tag_kmers(&self, k: usize) -> HashMap<Tag, usize> {
        match self {
            Self::Empty => HashMap::new(),
            Self::Unpaired { reads } => reads.tag_kmers(k),
            Self::Paired { paired1, paired2 } => {
                let mut hm_p1 = paired1.tag_kmers(k);
                let hm_p2 = paired2.tag_kmers(k);

                // combine the values for underlying Reads
                hm_p2.into_iter().for_each(|(data, kmers)| {
                   if let Some(count) = hm_p1.get_mut(&data) {
                    *count += kmers;
                   } else {
                    hm_p1.insert(data, kmers);
                   }
                });

                hm_p1
            },
            Self::Combined { paired1, paired2, unpaired } => {
                let mut hm_p1: HashMap<u8, usize> = paired1.tag_kmers(k);
                let hm_p2 = paired2.tag_kmers(k);
                let hm_up = unpaired.tag_kmers(k);

                // combine the values for underlying Reads
                hm_p2.into_iter().for_each(|(data, kmers)| {
                   if let Some(count) = hm_p1.get_mut(&data) {
                    *count += kmers;
                   } else {
                    hm_p1.insert(data, kmers);
                   }
                });

                hm_up.into_iter().for_each(|(data, kmers)| {
                    if let Some(count) = hm_p1.get_mut(&data) {
                     *count += kmers;
                    } else {
                     hm_p1.insert(data, kmers);
                    }
                 });

                hm_p1
            }
        }
    }

    /// return the number of k-mers occuring with each u8-encoded tag, 
    /// with the tag as the index
    /// if there are no tags saved in the Readspauired, it returns a vector of the
    /// length `n_sampeles`, filles with zeroes
    pub fn tag_kmers_vec(&self, k: usize, n_samples: usize) -> Vec<u64> {
        let hashed_kmer_counts = self.tag_kmers(k);

        let mut kmer_counts = vec![0; n_samples];

        for (tag, kmer_count) in hashed_kmer_counts {
            kmer_counts[tag as usize] += kmer_count as u64;
        }

        kmer_counts
    }
}

impl<D: Clone + Copy> Display for ReadsPaired<D> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Empty => write!(f, "empty ReadsPaired"),
            Self::Unpaired { reads } => write!(f, "unpaired ReadsPaired: \n{}", reads.info()),
            Self::Paired { paired1, paired2 } => write!(f, "paired ReadsPaired: \n{}\n{}", paired1.info(), paired2.info()),
            Self::Combined { paired1, paired2, unpaired } => write!(f, "combined ReadsPaired: \n{}\n{}\n{}", paired1.info(), paired2.info(), unpaired.info()),
        }
    }
}

/// Trait for ReadData, [`ID`]s can only be generated from [Marbel](https://github.com/jlab/marbel) reads
pub trait ReadData: PartialEq + Hash + serde::Serialize + DeserializeOwned + Debug + Clone + Copy + Eq + Send + Sync + Ord {
    /// geneate a read data, [`ID`]s and [`IDTag`]s can only be generated from [Marbel](https://github.com/jlab/marbel) reads
    fn read_data(gene_ids: &mut BiMap<String, ID>, read_name: &[u8], tag: Tag) -> Self;
    /// if available, get a tag
    fn get_tag(&self) -> Option<Tag>;
    /// retrun a ReadDatas enum to check which kind of ReadData is present
    fn read_datas() -> ReadDatas;
}

impl ReadData for Tag {
    fn read_data(_: &mut BiMap<String, ID>, _: &[u8], tag: Tag) -> Self {
        tag
    }

    fn get_tag(&self) -> Option<Tag> {
        Some(*self)
    }

    fn read_datas() -> ReadDatas {
        ReadDatas::Tag
    }
}

impl ReadData for ID {
    fn read_data(gene_ids: &mut BiMap<String, ID>, read_name: &[u8], _: Tag) -> Self {

        // read name is e.g. "B7R87_RS28825_2_0/1" -> gene: "B7R87_RS28825"
        // split at '_' and use first two elements and reconnect with '_'
        let read_name_sting = str::from_utf8(read_name).expect("error reading read name").to_string();
        let mut split_iter = read_name_sting.split('_');
        let mut gene = String::new();

        let Some(gene1) = split_iter.next() else {
            panic!("no gene names found in reads - only use id summarizers with marbel data - read name: {}", read_name_sting)
        };
        gene.push_str(gene1);

        gene.push('_');

        let Some(gene2) = split_iter.next() else {
            panic!("no gene names found in reads - only use id summarizers with marbel data - read name: {}", read_name_sting)
        };
        gene.push_str(gene2);

        // if gene not in gene_ids, then add, else get gene id
        let new_id = gene_ids.len() as ID;
        if new_id == ID::MAX { panic!("number of genes has surpassed 65.5k limit of u16") }
        match gene_ids.get_by_left(&gene) {
            Some(id) => *id,
            None => {
                gene_ids.insert(gene, new_id);
                new_id
            },
        }
    }

    fn get_tag(&self) -> Option<Tag> {
        None
    }

    fn read_datas() -> ReadDatas {
        ReadDatas::ID
    }
}

impl ReadData for IDTag {
    fn read_data(gene_ids: &mut BiMap<String, ID>, read_name: &[u8], tag: Tag) -> Self {
        let id = ID::read_data(gene_ids, read_name, tag);
        IDTag::new(id, tag)
    }

    fn get_tag(&self) -> Option<Tag> {
        Some(self.tag())
    }

    fn read_datas() -> ReadDatas {
        ReadDatas::IDTag
    }
}

#[derive(PartialEq, Eq, Debug, Clone, Copy, Serialize, Deserialize)]
pub enum ReadDatas {
    ID,
    Tag,
    IDTag
}

/// storage for which nodes two paired reads map to
#[derive(Debug, Clone, Default)]
pub struct ReadNodesPaired {
    r1_nodes: Vec<usize>,
    r1_node_positions: Vec<u8>, // reads are currently not expected to be longher than 256 bp
    r2_nodes: Vec<usize>,
    r2_node_positions: Vec<u8>,
    
} 

impl ReadNodesPaired {
    /// add node
    pub fn add(&mut self, node: usize, pos: u8, read_end: ReadEnd) {
        match read_end {
            ReadEnd::R1 => {
                self.r1_nodes.push(node);
                self.r1_node_positions.push(pos);
            }
            ReadEnd::R2 => {
                self.r2_nodes.push(node);
                self.r2_node_positions.push(pos);
            }
        }
    }

    /// get the nodes a read pair mapped to
    /// 
    /// nodes might appear twice if, for some reason, both reads mapped to it
    pub fn nodes(&self) -> Vec<usize> {
        let mut nodes = Vec::new();
        nodes.extend_from_slice(&self.r1_nodes);
        nodes.extend_from_slice(&self.r2_nodes);

        nodes
    }
}
/// storage for which reads (read pairs) map to a node
#[derive(Debug, Clone)]
#[derive(Default)]
pub struct NodeReadsPaired {
    reads: Vec<usize>,
    positions: Vec<(u8, ReadEnd)>, // reads are currently not expected to be longher than 256 bp
}

impl NodeReadsPaired {
    pub fn add(&mut self, read: usize, pos: u8, read_end: ReadEnd) {
        self.reads.push(read);
        self.positions.push((pos, read_end));
    }

    pub fn reads(&self) -> &[usize] {
        &self.reads
    }
}

pub struct MappedReads {
    nodes_per_read: Vec<ReadNodesPaired>, 
    reads_per_node: Vec<NodeReadsPaired>
}

impl MappedReads {
    pub fn new(nodes_per_read: Vec<ReadNodesPaired>, reads_per_node: Vec<NodeReadsPaired>) -> MappedReads {
        MappedReads {
            nodes_per_read,
            reads_per_node
        }
    }

    /// get the reads that mapped to a node
    pub fn reads_by_node(&self, node_id: usize) -> &[usize] {
        self.reads_per_node[node_id].reads()
    }

    /// get the nodes a read pair mapped to
    /// 
    /// nodes might appear twice if, for some reason, both reads mapped to it
    pub fn nodes_by_read(&self, read_id: usize) ->  Vec<usize> {
        self.nodes_per_read[read_id].nodes()
    }
}


#[cfg(test)]
mod tests {
    use std::{collections::HashMap, time};

    use bimap::BiMap;
    use itertools::enumerate;
    use rand::random;

    use crate::{dna_string::DnaString, reads::{Read, Strandedness}, summarizer::{IDTag, Tag, ID}, test::random_dna, Exts};
    use crate::reads::ReadData;
    use super::{Reads, ReadsPaired};

    #[test]
    fn test_add() {

        let fastq = vec![
            (DnaString::from_acgt_bytes(str::as_bytes("ACGATCGT")), Exts::empty(), 6u8),
            (DnaString::from_acgt_bytes(str::as_bytes("GGGGGG")), Exts::empty(), 5u8),
            (DnaString::from_acgt_bytes(str::as_bytes("TTGGTT")), Exts::empty(), 7u8),
            (DnaString::from_acgt_bytes(str::as_bytes("ACCAC")), Exts::empty(), 8u8),
            (DnaString::from_acgt_bytes(str::as_bytes("TCCCT")), Exts::empty(), 9u8),
            (DnaString::from_acgt_bytes(str::as_bytes("ACCAC")), Exts::empty(), 8u8),
            (DnaString::from_acgt_bytes(str::as_bytes("TCCCT")), Exts::empty(), 9u8),

        ];


        let mut reads = Reads::new(Strandedness::Unstranded);
        for (read, _, data) in fastq.clone() {
            reads.add_read(read, None, data);
        }

        println!("reads: {:#?}", reads);


        /* for no in reads.storage.iter() {
            println!("{:#b}", no)
        } */

         assert_eq!(reads.storage, vec![1791212948343256433, 5140577499666710528]);

        for (i, _) in fastq.iter().enumerate() {
            //println!("read {}: {:?}", i, reads.get_read(i))
            assert_eq!(Read::new(fastq[i].0.clone(), fastq[i].1, fastq[i].2, Strandedness::Unstranded), reads.get_read(i).unwrap())
        }

        for read in reads.iter() {
            let seq = read.seq;
            println!("{:?}, {}", seq, seq.len())
        }
        println!();

        for read in reads.partial_iter(5..7) {
            println!("{:?}", read)
        }

        println!("memory usage: {}", reads.mem())
    }

    #[test]
    fn test_get_read() {
        let mut reads = Reads::new(Strandedness::Unstranded);
        //reads.add_read(DnaString::from_acgt_bytes("AGCTAGCTAGC".as_bytes()), Exts::empty(), 67u8);
        reads.add_from_bytes("ACGATCGNATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACTGATCGATCG".as_bytes(), None, 67u8);
        let read = reads.get_read(0);
        println!("{:?}", read);
        println!("{:#066b}", reads.storage[0]);
        println!("{:?}", reads)
    }

    #[test]
    fn test_add_from_bytes() {
        let dna = [
            "AAGCGGAGATTATTCACGAGCATCGCGTAC".as_bytes(),
            "GATCGATGCATGCTAGA".as_bytes(),
            "ACGTAAAAAAAAAATTATATAACGTACGTAAAAAAAAAATTATATAACGTAACGTAAAAAAAAAAATTATAATAACGT".as_bytes(),
            "AGCTAGCTAGCTGACTGAGCGACTGA".as_bytes(),
            "AGCTAGCTAGCTGACTGAGCGACTGACGGATC".as_bytes(),
            "TTTTTTTTTTTTTTTTTTTTTTTT".as_bytes(),
            "ACGATCGAATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACTGATCGATCG".as_bytes(),
            "ACGATCGATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACTGATCGATCGAAGGGCAGTTAGGCCGTAAGCGCGAT".as_bytes(),
        ];

        let mut reads: Reads<u8> = Reads::new(Strandedness::Unstranded);
        for seq in dna {
            reads.add_from_bytes(seq, None, random());

        }

        for (i, read) in enumerate(reads.iter()) {
            let sequence = DnaString::from_acgt_bytes(dna[i]);
            assert_eq!(read.seq, sequence);
        }
    }

    #[test]
    fn test_add_from_bytes_checked() {
        let dna = [
            "AAGCGGAGATTATTCACGAGCATCGCGTAC".as_bytes(),
            "GATCGATGCATGCTAGA".as_bytes(),
            "ACGTAAAAAAAAAATTATATAACGTACGTAAAAAAAAAANTTATATAACGTAACGTAAAAAAAAAAATTATAATAACGT".as_bytes(),
            "AGCTAGCTAGCTGACNGAGCGACTGA".as_bytes(),
            "AGCTAGCTAGCTGACTGAGCGACTGACGGATC".as_bytes(),
            "TTTTTTTTTTTTTTTTTTTTTTTT".as_bytes(),
            "ACGATCGAATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACNNNTGATCGATCG".as_bytes(),
            "ACGATCGATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACTGATCGATCGAAGGGCAGTTAGGCCGTAAGCGCGAT".as_bytes(),
            "A".as_bytes(),
            "AAAAN".as_bytes(),
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN".as_bytes(),
        ];

        let mut reads: Reads<u8> = Reads::new(Strandedness::Unstranded);
        let mut corrects = Vec::new();
        for seq in dna {
            corrects.push(reads.add_from_bytes_checked(seq, None, random()));
        }

        let mut read_counter = 0;

        for (i, correct) in enumerate(corrects) {
            let sequence = DnaString::from_acgt_bytes_checked(dna[i]);
            match correct {
                true => {
                    let read = reads.get_read(read_counter).unwrap();
                    assert_eq!(sequence.unwrap(), read.seq);
                    read_counter += 1;

                },
                false => assert!(sequence.is_err()),
            }
        }
    }

    #[test]
    fn test_speed_from_bytes() {
        let dnas = [
            "AAGCGGAGATTATTCACGAGCATCGCGTAC".as_bytes(),
            "GATCGATGCATGCTAGA".as_bytes(),
            "ACGTAAAAAAAAAATTATATAACGTACGTAAAAAAAAAANTTATATAACGTAACGTAAAAAAAAAAATTATAATAACGT".as_bytes(),
            "AGCTAGCTAGCTGACNGAGCGACTGA".as_bytes(),
            "AGCTAGCTAGCTGACTGAGCGACTGACGGATC".as_bytes(),
            "TTTTTTTTTTTTTTTTTTTTTTTT".as_bytes(),
            "ACGATCGAATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACNNNTGATCGATCG".as_bytes(),
            "ACGATCGATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACTGATCGATCGAAGGGCAGTTAGGCCGTAAGCGCGAT".as_bytes(),
            "A".as_bytes(),
            "AAAAN".as_bytes(),
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN".as_bytes(),
        ];

        const REPS: usize = 500000;
        /*
        test with 5 ml:
            through DnaString: 88.09323 s 
            direct to Read: 44.41913 s
         */


        let ds_start= time::Instant::now();
        let mut reads: Reads<u8> = Reads::new(Strandedness::Unstranded);
        for _i in 0..REPS {
            for dna in dnas {
                reads.add_read(DnaString::from_acgt_bytes(dna), None, random());
            }
        }
        let ds_finish = ds_start.elapsed();

        let r_start= time::Instant::now();
        let mut reads: Reads<u8> = Reads::new(Strandedness::Unstranded);
        for _i in 0..REPS {
            for dna in dnas {
                reads.add_from_bytes(dna, None, random());
            }
        }
        let r_finish = r_start.elapsed();

        println!("through DnaString: {} s \n direct to Read: {} s", ds_finish.as_secs_f32(), r_finish.as_secs_f32())


    }


    #[test]
    fn test_reads_stranded() {
        let mut reads: Reads<u8> = Reads::new(Strandedness::Forward);
        assert_eq!(reads.stranded(), Strandedness::Forward);
        reads.set_stranded(Strandedness::Reverse);
        assert_eq!(reads.stranded(), Strandedness::Reverse);

    }

    #[test]
    fn test_reads_data_kmers() {
        let mut reads = Reads::new(Strandedness::Unstranded);
        let seqs = [
            ("ACGATCGTACGTACGTAGCTAGCTGCTAGCTAGCTGACTGACTGA", 0),
            ("CGATGCTATCAGCGAGCGATCGTACGTAGCTACG", 1),
            ("CGATCGACGAGCAGCGTATGCTACGAGCTGACGATCTACGA", 2),
            ("CACACACGGCATCGATCGAGCAGCATCGACTACGTA", 3),
        ];

        seqs.iter().for_each(|(read, tag)| reads.add_from_bytes(read.as_bytes(), None, *tag as u8));
        let data_kmers = reads.tag_kmers(16);
       
        let comp_hm: HashMap<u8, usize> = [(0, 30), (1, 19), (2, 26), (3, 21)].into_iter().collect();

        assert_eq!(comp_hm, data_kmers);
    }

    #[test]
    fn test_reads_add_exts() {
        let mut raw_reads = Vec::new();
        for _i in 0..10 {
            raw_reads.push((DnaString::from_bytes(&random_dna(100)), Exts::new(rand::random::<u8>()), rand::random::<u8>()));
        }
        let reads = Reads::from_vmer_vec(raw_reads.clone(), Strandedness::Unstranded);
        let new_raw_reads = reads.iter().map(|read| (read.seq, read.exts, read.data)).collect::<Vec<_>>();

        assert_eq!(raw_reads, new_raw_reads)
    }

    #[test]
    fn test_reads_info() {
        let mut reads = Reads::new(Strandedness::Unstranded);
        let seqs = [
            ("ACGATCGTACGTACGTAGCTAGCTGCTAGCTAGCTGACTGACTGA", 0),
            ("CGATGCTATCAGCGAGCGATCGTACGTAGCTACG", 1),
            ("CGATCGACGAGCAGCGTATGCTACGAGCTGACGATCTACGA", 2),
            ("CACACACGGCATCGATCGAGCAGCATCGACTACGTA", 3),
        ];

        seqs.iter().for_each(|(read, tag)| reads.add_from_bytes(read.as_bytes(), None, *tag as u8));

        assert_eq!(reads.info(), "Reads { n reads: 4, stranded: Unstranded }".to_string());
    }

    #[test]
    fn test_reads_paired() {
        let mut p1 = Reads::new(Strandedness::Unstranded);
        let mut p2 = Reads::new(Strandedness::Unstranded);
        let mut up = Reads::new(Strandedness::Unstranded);

        let reads_p1 = [
            "ACGATCGTACGTACGTAGCTAGCTGCTAGCTAGCTGACTGACTGA",
            "CGATGCTATCAGCGAGCGATCGTACGTAGCTACG",
            "CGATCGACGAGCAGCGTATGCTACGAGCTGACGATCTACGA",
            "CACACACGGCATCGATCGAGCAGCATCGACTACGTA",
        ];

        let reads_p2 = [
            "AGCTAGCTAGCTACTGATCGTAGCTAGCTGATCGA",
            "AGCGATCGTACGTAGCTAGCTA",
            "CGATCGATCGACTAGCGTAGCTGACTGAC",
            "CAGATGCTCTGCTGACTGACTGATCGTACTGACTAGCATCTAGC",
        ];

        let reads_up = [
            "CGTACTAGCTGACGTAC",
            "CGATGCTAGCTAGCTAGCGATCG",
        ];

        let tags = (0..4).collect::<Vec<u8>>();

        reads_p1.iter().enumerate().for_each(|(i, read)| p1.add_from_bytes(read.as_bytes(), None, tags[i]));
        reads_p2.iter().enumerate().for_each(|(i, read)| p2.add_from_bytes(read.as_bytes(), None, tags[i]));
        reads_up.iter().enumerate().for_each(|(i, read)| up.add_from_bytes(read.as_bytes(), Some(Exts::new(i as u8)), tags[i]));

        let empty: ReadsPaired<u8> = ReadsPaired::from_reads((Reads::new(Strandedness::Unstranded), Reads::new(Strandedness::Unstranded), Reads::new(Strandedness::Unstranded)));
        assert_eq!(empty, ReadsPaired::Empty);
        assert_eq!(empty.mem(), 0);
        assert_eq!(empty.n_reads(), 0);
        assert_eq!(empty.iterable(), Vec::<&Reads<u8>>::new());

        let unpaired = ReadsPaired::from_reads((Reads::new(Strandedness::Unstranded), Reads::new(Strandedness::Unstranded), up.clone()));
        assert_eq!(ReadsPaired::Unpaired { reads: up.clone() }, unpaired);
        println!("exts up {:?}", unpaired);
        assert_eq!(unpaired.mem(), 148);
        assert_eq!(unpaired.n_reads(), 2);
        assert_eq!(unpaired.iterable(), vec![&up]);
      
        let paired = ReadsPaired::from_reads((p1.clone(), p2.clone(), Reads::new(Strandedness::Unstranded)));
        assert_eq!(ReadsPaired::Paired { paired1: p1.clone(), paired2: p2.clone() }, paired);        
        assert_eq!(paired.mem(), 376);
        assert_eq!(paired.n_reads(), 8);
        assert_eq!(paired.iterable(), vec![&p1, &p2]);

        let combined = ReadsPaired::from_reads((p1.clone(), p2.clone(), up.clone()));
        assert_eq!(ReadsPaired::Combined { paired1: p1.clone(), paired2: p2.clone(), unpaired: up.clone() }, combined);
        assert_eq!(combined.mem(), 524);
        assert_eq!(combined.n_reads(), 10);
        assert_eq!(combined.iterable(), vec![&p1, &p2, &up]);


        // test iter

        assert_eq!(unpaired.iter().collect::<Vec<_>>(), up.iter().collect::<Vec<_>>());
        assert_eq!(paired.iter().collect::<Vec<_>>(), p1.iter().chain(p2.iter()).collect::<Vec<_>>());
        assert_eq!(combined.iter().collect::<Vec<_>>(), p1.iter().chain(p2.iter()).chain(up.iter()).collect::<Vec<_>>());

        // test partial iter

        assert_eq!(unpaired.iter_partial(0..1).collect::<Vec<_>>(), up.partial_iter(0..1).collect::<Vec<_>>());

        assert_eq!(paired.iter_partial(0..1).collect::<Vec<_>>(), p1.partial_iter(0..1).collect::<Vec<_>>());
        assert_eq!(paired.iter_partial(5..7).collect::<Vec<_>>(), p2.partial_iter(1..3).collect::<Vec<_>>());
        assert_eq!(paired.iter_partial(1..8).collect::<Vec<_>>(), p1.partial_iter(1..4).chain(p2.partial_iter(0..4)).collect::<Vec<_>>());

        assert_eq!(combined.iter_partial(0..1).collect::<Vec<_>>(), p1.partial_iter(0..1).collect::<Vec<_>>());
        assert_eq!(combined.iter_partial(5..7).collect::<Vec<_>>(), p2.partial_iter(1..3).collect::<Vec<_>>());
        assert_eq!(combined.iter_partial(8..10).collect::<Vec<_>>(), up.partial_iter(0..2).collect::<Vec<_>>());
        assert_eq!(combined.iter_partial(1..8).collect::<Vec<_>>(), p1.partial_iter(1..4).chain(p2.partial_iter(0..4)).collect::<Vec<_>>());
        assert_eq!(combined.iter_partial(6..9).collect::<Vec<_>>(), p2.partial_iter(2..4).chain(up.partial_iter(0..1)).collect::<Vec<_>>());
        assert_eq!(combined.iter_partial(1..9).collect::<Vec<_>>(), p1.partial_iter(1..4).chain(p2.partial_iter(0..4)).chain(up.partial_iter(0..1)).collect::<Vec<_>>());

        // test tag kmers (and data_kmers)
        assert_eq!(unpaired.tag_kmers_vec(16, 2), vec![2, 8]);
        assert_eq!(paired.tag_kmers_vec(16, 4), vec![50, 26, 40, 50]);
        assert_eq!(combined.tag_kmers_vec(16, 4), vec![52, 34, 40, 50]);
        assert_eq!(empty.tag_kmers_vec(16, 0), Vec::<u64>::new());

        // test decombine
        let mut paired_dc = paired.clone();
        let rm_reads_p = paired_dc.decombine();
        let mut combined_dc = combined.clone();
        let rm_reads_up = combined_dc.decombine();
        assert_eq!(paired_dc, paired);
        assert_eq!(rm_reads_p, 0);
        assert_eq!(combined_dc, paired);
        assert_eq!(rm_reads_up, 2);

        // display
        assert_eq!(format!("{}", empty), "empty ReadsPaired".to_string());
        assert_eq!(format!("{}", unpaired), "unpaired ReadsPaired: 
Reads { n reads: 2, stranded: Unstranded }".to_string());
        assert_eq!(format!("{}", paired), "paired ReadsPaired: 
Reads { n reads: 4, stranded: Unstranded }
Reads { n reads: 4, stranded: Unstranded }".to_string());
        assert_eq!(format!("{}", combined), "combined ReadsPaired: 
Reads { n reads: 4, stranded: Unstranded }
Reads { n reads: 4, stranded: Unstranded }
Reads { n reads: 2, stranded: Unstranded }".to_string());

    }

    #[test]
    #[should_panic]
    fn test_reads_paired_panic() {
        let mut p1 = Reads::new(Strandedness::Unstranded);

        let reads_p1 = [
            "ACGATCGTACGTACGTAGCTAGCTGCTAGCTAGCTGACTGACTGA",
            "CGATGCTATCAGCGAGCGATCGTACGTAGCTACG",
            "CGATCGACGAGCAGCGTATGCTACGAGCTGACGATCTACGA",
            "CACACACGGCATCGATCGAGCAGCATCGACTACGTA",
        ];

        reads_p1.iter().for_each(|read| p1.add_from_bytes(read.as_bytes(), None, 0u8));

        let _ = ReadsPaired::from_reads((p1, Reads::new(Strandedness::Unstranded), Reads::new(Strandedness::Unstranded)));
    }

    #[test]
    fn test_read_data() {
        let read_name_1 = "B7R87_RS28825_2_0/1".as_bytes(); // gene B7R87_RS28825
        let read_name_2 = "B7R87_RS21825_2_0/1".as_bytes(); // gene B7R87_RS21825

        let tag = 0 as Tag;
        
        let mut ids = BiMap::new();

        let id = ID::read_data(&mut ids, read_name_1, tag);
        assert_eq!(id, 0);

        let id_tag = IDTag::read_data(&mut ids, read_name_1, tag);
        assert_eq!(id_tag, IDTag::new(id, tag));
        assert_eq!(id_tag.get_tag(), Some(tag));

        let id = ID::read_data(&mut ids, read_name_2, tag);
        assert_eq!(id, 1);
        assert_eq!(id.get_tag(), None);

        assert_eq!(Tag::read_data(&mut ids, read_name_1, tag), tag);
        assert_eq!(tag.get_tag(), Some(tag));
    }

    #[test]
    #[should_panic]
    fn test_read_data_panic() {
        let mut ids = BiMap::new();

        let read_name_1 = "B7R87_RS28825_2_0/1".as_bytes(); // gene B7R87_RS28825
        let read_name_2 = "B7R87_RS21825_2_0/1".as_bytes(); // gene B7R87_RS21825

        let _ = ID::read_data(&mut ids, read_name_1, 0);
        let _ = ID::read_data(&mut ids, read_name_2, 0);

        // trying to add invalid read name -> panic
        let _ = ID::read_data(&mut ids, "AAAAAA".as_bytes(), 0);
    }
}
