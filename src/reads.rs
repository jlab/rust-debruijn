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

/// Store many DNA sequences together with an Exts and data each compactly packed together
/// 
/// #### fields:
/// 
/// * `storage`: `Vec` with 2-bit encoded DNA bases of all sequences
/// * `ends`:  `Vec` with the ends (exclusive) of the separate sequences in the `Reads`
/// * `exts`: `Vec` with one Exts for each sequence
/// * `data`: `Vec` with data for each sequence
/// * `len`: length of all sequences together
/// * `stranded`: [`Stranded`] conveying the strandedness and direction of the reads
#[derive(Ord, PartialOrd, Clone, PartialEq, Eq, Hash, Serialize, Deserialize, Debug)]
pub struct Reads<D> {
    storage: Vec<u64>,
    ends: Vec<usize>,
    exts: Vec<Exts>,
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
            exts: Vec::new(),
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
        mem::size_of_val(self) + size_of_val(&*self.storage) + size_of_val(&*self.data) + size_of_val(&*self.ends) + size_of_val(&*self.exts)
    }

    /// set the strandedness and the direction of the reads
    pub fn set_stranded(&mut self, stranded: Strandedness) {
        self.stranded = stranded
    }

    /// Adds a new read to the `Reads`
    // maybe push_base until u64 is full and then do extend like in DnaString::extend ? with accellerated mode
    pub fn add_read<V: Vmer>(&mut self, read: V, exts: Exts, data: D) {
        for base in read.iter() {
            self.push_base(base);
        }
        self.ends.push(self.len);
        self.data.push(data);
        self.exts.push(exts);
       
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
            reads.data.push(data);
            reads.exts.push(exts);
        }

        reads.storage.shrink_to_fit();
        reads.data.shrink_to_fit();
        reads.exts.shrink_to_fit();
        reads.ends.shrink_to_fit();
        
        reads
    }


    /// add ASCII encoded bases to the Reads
    /// 
    /// will transform all ascii characters outside of ACGTacgt into A
    /// see also: [`Reads::add_from_bytes_checked`]
    pub fn add_from_bytes(&mut self, bytes: &[u8], exts: Exts, data: D) {
        
        // fill the last incomplete u64 block
        let missing = 32 - (self.len % 32);
        if missing != 0 {
            if  missing > bytes.len() {
                let fill = bytes.iter().map(|c| base_to_bits(*c));
                self.extend(fill);
                self.exts.push(exts);
                self.data.push(data);
                self.ends.push(self.len);
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

                self.exts.push(exts);
                self.data.push(data);
                self.ends.push(self.len);
                return;
            }
        }

        let b = bytes.iter().map(|c| base_to_bits(*c));
        self.extend(b);

        self.exts.push(exts);
        self.data.push(data);
        self.ends.push(self.len);
        
    }

    /// add ASCII encoded bases to the Reads
    /// 
    /// will return `false` if the bytes contained characters outside of `ACGTacgt`, otherwise return true and add the bases
    /// see also: [`Reads::add_from_bytes`]
    pub fn add_from_bytes_checked(&mut self, bytes: &[u8], exts: Exts, data: D) -> bool {

        let (_, corrects): (Vec<u8>, Vec<bool>) = bytes.iter().map(|c| base_to_bits_checked(*c)).collect();
        if corrects.iter().contains(&false) { return false }

        
        // fill the last incomplete u64 block
        let missing = 32 - (self.len % 32);
        if missing != 0 {
            if  missing > bytes.len() {
                let fill = bytes.iter().map(|c| base_to_bits(*c));
                self.extend(fill);
                self.exts.push(exts);
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

                self.exts.push(exts);
                self.data.push(data);
                self.ends.push(self.len);
                return true;
            }
        }

        let b = bytes.iter().map(|c| base_to_bits(*c));
        self.extend(b);

        self.exts.push(exts);
        self.data.push(data);
        self.ends.push(self.len);

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
    pub fn get_read(&self, i: usize) -> Option<(DnaString, Exts, D, Strandedness)> {
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

        Some((sequence, self.exts[i], self.data[i], self.stranded))
    }


    /// shrink the vectors' capacity to fit the length
    /// 
    /// use sparsely
    pub fn shrink_to_fit(&mut self)  {
        self.storage.shrink_to_fit();
        self.data.shrink_to_fit();
        self.exts.shrink_to_fit();
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

        self.iter().for_each(|(read, _, data, _)| {
            let kmers = read.len().saturating_sub(k - 1);
            if let Some(tag) = data.get_tag() {
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
    type Item = (DnaString, Exts, D, Strandedness);

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
        let vec: Vec<(DnaString, Exts, D, Strandedness)> = self.iter().collect();
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
    pub fn iterable(&self) -> Vec<&Reads<D>> {
        match self {
            Self::Empty => vec![],
            Self::Unpaired { reads  } => vec![reads],
            Self::Paired { paired1, paired2 } => vec![paired1, paired2],
            Self::Combined { paired1, paired2, unpaired } => vec![paired1, paired2, unpaired],
        }
    }

    pub fn n_reads(&self) -> usize {
        match self {
            Self::Empty => 0,
            Self::Unpaired { reads  } => reads.n_reads(),
            Self::Paired { paired1, paired2 } => paired1.n_reads() + paired2.n_reads(),
            Self::Combined { paired1, paired2, unpaired } => paired1.n_reads() + paired2.n_reads() + unpaired.n_reads(),
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

    pub fn iter(&self) -> Box<dyn Iterator<Item = (DnaString, Exts, D, Strandedness)> + '_> {
        match self {
            ReadsPaired::Empty => panic!("Error: no reads to process"),
            ReadsPaired::Unpaired { reads } => Box::new(reads.iter()),
            ReadsPaired::Paired { paired1, paired2 } => Box::new(paired1.iter().chain(paired2.iter())),
            ReadsPaired::Combined { paired1, paired2, unpaired } => Box::new(paired1.iter().chain(paired2.iter()).chain(unpaired.iter())),
        }
    }

    pub fn iter_partial(&self, range: Range<usize>) -> Box<dyn Iterator<Item = (DnaString, Exts, D, Strandedness)> + '_> {
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

/// Trait for ReadData
pub trait ReadData: PartialEq + Hash + serde::Serialize + DeserializeOwned + Debug + Clone + Copy + Eq + Send + Sync + Ord {
    fn read_data(gene_ids: &mut BiMap<String, ID>, read_name: &[u8], tag: Tag) -> Self;
    fn get_tag(&self) -> Option<Tag>;
}

impl ReadData for Tag {
    fn read_data(_: &mut BiMap<String, ID>, _: &[u8], tag: Tag) -> Self {
        tag
    }

    fn get_tag(&self) -> Option<Tag> {
        Some(*self)
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
            panic!("no gene names found in reads - only use id-sum summarizer with marbel data - read name: {}", read_name_sting)
        };
        gene.push_str(gene1);

        gene.push('_');

        let Some(gene2) = split_iter.next() else {
            panic!("no gene names found in reads - only use id-sum summarizer with marbel data - read name: {}", read_name_sting)
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
}

impl ReadData for IDTag {
    fn read_data(gene_ids: &mut BiMap<String, ID>, read_name: &[u8], tag: Tag) -> Self {
        let id = ID::read_data(gene_ids, read_name, tag);
        IDTag::new(id, tag)
    }

    fn get_tag(&self) -> Option<Tag> {
        Some(self.tag())
    }
}

#[cfg(test)]
mod tests {
    use std::{collections::HashMap, time};

    use itertools::enumerate;
    use rand::random;

    use crate::{dna_string::DnaString, reads::Strandedness, Exts};
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
        for (read, ext, data) in fastq.clone() {
            reads.add_read(read, ext, data);
        }

        println!("reads: {:#?}", reads);


        /* for no in reads.storage.iter() {
            println!("{:#b}", no)
        } */

         assert_eq!(reads.storage, vec![1791212948343256433, 5140577499666710528]);

        for (i, _) in fastq.iter().enumerate() {
            //println!("read {}: {:?}", i, reads.get_read(i))
            assert_eq!((fastq[i].0.clone(), fastq[i].1, fastq[i].2, Strandedness::Unstranded), reads.get_read(i).unwrap())
        }

        for (seq, _, _, _) in reads.iter() {
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
        reads.add_from_bytes("ACGATCGNATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACTGATCGATCG".as_bytes(), Exts::empty(), 67u8);
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
            reads.add_from_bytes(seq, Exts::empty(), random());

        }

        for (i, seq) in enumerate(reads.iter()) {
            let sequence = DnaString::from_acgt_bytes(dna[i]);
            assert_eq!(seq.0, sequence);
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
            corrects.push(reads.add_from_bytes_checked(seq, Exts::empty(), random()));
        }

        let mut read_counter = 0;

        for (i, correct) in enumerate(corrects) {
            let sequence = DnaString::from_acgt_bytes_checked(dna[i]);
            match correct {
                true => {
                    let read = reads.get_read(read_counter).unwrap();
                    assert_eq!(sequence.unwrap(), read.0);
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
                reads.add_read(DnaString::from_acgt_bytes(dna), Exts::empty(), random());
            }
        }
        let ds_finish = ds_start.elapsed();

        let r_start= time::Instant::now();
        let mut reads: Reads<u8> = Reads::new(Strandedness::Unstranded);
        for _i in 0..REPS {
            for dna in dnas {
                reads.add_from_bytes(dna, Exts::empty(), random());
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

        seqs.iter().for_each(|(read, tag)| reads.add_from_bytes(read.as_bytes(), Exts::empty(), *tag as u8));
        let data_kmers = reads.tag_kmers(16);
       
        let comp_hm: HashMap<u8, usize> = [(0, 30), (1, 19), (2, 26), (3, 21)].into_iter().collect();

        assert_eq!(comp_hm, data_kmers);
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

        seqs.iter().for_each(|(read, tag)| reads.add_from_bytes(read.as_bytes(), Exts::empty(), *tag as u8));

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

        reads_p1.iter().enumerate().for_each(|(i, read)| p1.add_from_bytes(read.as_bytes(), Exts::empty(), tags[i]));
        reads_p2.iter().enumerate().for_each(|(i, read)| p2.add_from_bytes(read.as_bytes(), Exts::empty(), tags[i]));
        reads_up.iter().enumerate().for_each(|(i, read)| up.add_from_bytes(read.as_bytes(), Exts::empty(), tags[i]));

        let empty: ReadsPaired<u8> = ReadsPaired::from_reads((Reads::new(Strandedness::Unstranded), Reads::new(Strandedness::Unstranded), Reads::new(Strandedness::Unstranded)));
        assert_eq!(empty, ReadsPaired::Empty);
        assert_eq!(empty.mem(), 0);
        assert_eq!(empty.n_reads(), 0);
        assert_eq!(empty.iterable(), Vec::<&Reads<u8>>::new());

        let unpaired = ReadsPaired::from_reads((Reads::new(Strandedness::Unstranded), Reads::new(Strandedness::Unstranded), up.clone()));
        assert_eq!(ReadsPaired::Unpaired { reads: up.clone() }, unpaired);
        assert_eq!(unpaired.mem(), 148);
        assert_eq!(unpaired.n_reads(), 2);
        assert_eq!(unpaired.iterable(), vec![&up]);
      
        let paired = ReadsPaired::from_reads((p1.clone(), p2.clone(), Reads::new(Strandedness::Unstranded)));
        assert_eq!(ReadsPaired::Paired { paired1: p1.clone(), paired2: p2.clone() }, paired);        
        assert_eq!(paired.mem(), 384);
        assert_eq!(paired.n_reads(), 8);
        assert_eq!(paired.iterable(), vec![&p1, &p2]);

        let combined = ReadsPaired::from_reads((p1.clone(), p2.clone(), up.clone()));
        assert_eq!(ReadsPaired::Combined { paired1: p1.clone(), paired2: p2.clone(), unpaired: up.clone() }, combined);
        assert_eq!(combined.mem(), 532);
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

        reads_p1.iter().for_each(|read| p1.add_from_bytes(read.as_bytes(), Exts::empty(), 0u8));

        let _ = ReadsPaired::from_reads((p1, Reads::new(Strandedness::Unstranded), Reads::new(Strandedness::Unstranded)));
    }
}
