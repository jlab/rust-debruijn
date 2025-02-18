use bimap::BiMap;
use serde::{Deserialize, Serialize};
use crate::{Exts, Kmer, Tags};
use std::{fmt::Debug, marker::PhantomData, mem};

#[cfg(not(feature = "sample128"))]
pub type M = u64;

#[cfg(feature = "sample128")]
pub type M = u128;

#[derive(Copy, Clone, Serialize, Deserialize, Debug)]
pub struct Marker {
    marker0: M,
    marker1: M,
    count0: u8,
    count1: u8,
}

impl Marker {
    pub fn new(marker0: M, marker1: M, count0: u8, count1: u8) -> Self {
        Marker { marker0, marker1, count0, count1 }
    }
}

/// round an unsigned integer to the specified amount of digits,
/// if the integer is shorter than the number if digits, it returns the original integer
pub fn round_digits(number: u32, digits: u32) -> u32 {
    let length = (number as f32).log10() as u32 + 1;
    if digits > length { return number }
    let empty = length - digits;
    ((number as f32/ 10i32.pow(empty) as f32).round() * 10i32.pow(empty) as f32) as u32
}

/// Trait for the output of the KmerSummarizers
pub trait SummaryData<D> {
    /// Make a new `SummaryData<D>`
    fn new(data: D) -> Self;
    /// does not actually print but format
    fn print(&self, tag_translator: &BiMap<String, u8>) -> String;
    /// If the `SummaryData` contains sufficient information, return `Vec<u8>` of the tags and the count
    fn vec_for_color(&self) -> Option<(Vec<u8>, i32)>;
    /// If the `SummaryData` contains sufficient information, return the Tags and the count 
    fn get_tags_sum(&self) -> Option<(Tags, i32)>;
    /// return a score (the sum of the kmer appearances), `Vec<D>` simply returns `1`
    fn score(&self) -> f32;
    /// return the size of the structure, including contents of slices
    fn mem(&self) -> usize;
}

impl<> SummaryData<u32> for u32 {
    fn new(data: u32) -> Self {
        data
    }

    fn print(&self, _: &BiMap<String, u8>) -> String {
        format!("count: {}", self).replace("\"", "\'")
    }
    fn vec_for_color(&self) -> Option<(Vec<u8>, i32)> {
        None
    }

    fn get_tags_sum(&self) -> Option<(Tags, i32)> {
        None
    }

    fn score(&self) -> f32 {
        *self as f32
    }

    fn mem(&self) -> usize {
        mem::align_of_val(&*self)
    }

}

impl<D: Debug> SummaryData<Vec<D>> for Vec<D> {
    fn new(data: Vec<D>) -> Self {
        data
    }

    fn print(&self, _: &BiMap<String, u8>) -> String {
        format!("tags: {:?}", self).replace("\"", "\'")
    }
    
    fn vec_for_color(&self) -> Option<(Vec<u8>, i32)> {
        None
    }
    
    fn get_tags_sum(&self) -> Option<(Tags, i32)> {
        None
    }

    fn score(&self) -> f32 {
        1.
    }

    fn mem(&self) -> usize {
        mem::size_of_val(&**self) + mem::size_of_val(&*self)
    }

}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
// aligned would be 16 Bytes, packed is 12 Bytes
#[repr(packed)]
pub struct TagsSumData {
    tags: Tags,
    sum: i32,
}

impl SummaryData<(Tags, i32)> for TagsSumData {
    fn new(data: (Tags, i32)) -> Self {
        TagsSumData { tags: data.0, sum: data.1 }
    }

    fn print(&self, tag_translator: &BiMap<String, u8>) -> String {
        // need to copy fields to local variable because repr(packed) results in unaligned struct
        let tags = self.tags;
        let sum = self.sum;
        // replace " with ' to avoid conflicts in dot file
        format!("tags: {:?}, sum: {}", tags.to_string_vec(tag_translator), sum).replace("\"", "\'")
    }

    fn vec_for_color(&self) -> Option<(Vec<u8>, i32)> {
        // need to copy fields to local variable because repr(packed) results in unaligned struct
        let tags = self.tags;
        Some((tags.to_u8_vec(), self.sum))
    }

    fn get_tags_sum(&self) -> Option<(Tags, i32)> {
        Some((self.tags, self.sum))
    }

    fn score(&self) -> f32 {
        self.sum as f32
    }

    fn mem(&self) -> usize {
        mem::size_of_val(&*self)
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TagsCountData {
    tags: Tags,
    sum: i32,
    counts: Box<[u32]>,
}

impl SummaryData<(Tags, Box<[u32]>, i32)> for TagsCountData {
    fn new(data: (Tags, Box<[u32]>, i32)) -> Self {
        TagsCountData { tags: data.0, counts: data.1, sum: data.2 }
    }

    fn print(&self, tag_translator: &BiMap<String, u8>) -> String {
        format!("tags: {:?}, counts: {:?}, sum: {}", self.tags.to_string_vec(tag_translator), self.counts, self.sum).replace("\"", "\'")
    }

    fn vec_for_color(&self) -> Option<(Vec<u8>, i32)> {
        Some((self.tags.to_u8_vec(), self.sum))
    }

    fn get_tags_sum(&self) -> Option<(Tags, i32)> {
        Some((self.tags, self.sum))
    }

    fn score(&self) -> f32 {
        self.sum as f32
    }

    fn mem(&self) -> usize {
        mem::size_of_val(&*self) + mem::size_of_val(&*self.counts)
    }

}

/// Structure for Summarizer Data which contains the tags and the count for each tag
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TagsCountsData {
    tags: Tags,
    counts: Box<[u32]>
}

impl TagsCountsData {
    #[inline]
    pub fn sum(&self) -> i32 {
        self.counts.iter().sum::<u32>() as i32
    }
}

impl SummaryData<(Tags, Box<[u32]>)> for TagsCountsData{
    fn new(data: (Tags, Box<[u32]>)) -> Self {
        TagsCountsData { tags: data.0, counts: data.1 }
    }

    fn print(&self, tag_translator: &BiMap<String, u8>) -> String {
        format!("tags: {:?}, counts: {:?}", self.tags.to_string_vec(tag_translator), self.counts).replace("\"", "\'")
    }

    fn vec_for_color(&self) -> Option<(Vec<u8>, i32)> {
        Some((self.tags.to_u8_vec(), self.sum()))
    }

    fn get_tags_sum(&self) -> Option<(Tags, i32)> {
        Some((self.tags, self.sum()))
    }

    fn score(&self) -> f32 {
        self.sum() as f32
    }

    fn mem(&self) -> usize {
        mem::size_of_val(&*self) + mem::size_of_val(&*self.counts)
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GroupCountData {
    group1: u32,
    group2: u32,
}

impl GroupCountData {
    #[inline(always)]
    fn sum(&self) -> u32 {
        self.group1 + self.group2
    }
}

impl SummaryData<(u32, u32)> for GroupCountData {
    fn new(data: (u32, u32)) -> Self {
        GroupCountData { group1: data.0, group2: data.1 }
    }
    
    fn print(&self, _: &BiMap<String, u8>) -> String {
        format!("count 1: {}, count 2: {}", self.group1, self.group2)
    }

    fn vec_for_color(&self) -> Option<(Vec<u8>, i32)> {
        None
    }

    fn get_tags_sum(&self) -> Option<(Tags, i32)> {
        None
    }

    fn score(&self) -> f32 {
        self.sum() as f32
    }

    fn mem(&self) -> usize {
        mem::size_of_val(&*self)
    }
}


/// Data container to store the relative amount of k-mers (in percent) in group 1 and the overall count
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RelCountData {
    percent: u32,
    count: u32
}

impl SummaryData<(u32, u32)> for RelCountData {
    fn new(data: (u32, u32)) -> Self {
        RelCountData { percent: data.0, count: data.1 }
    }
    
    fn print(&self, _: &BiMap<String, u8>) -> String {
        format!("relative amount group 1: {}, count both: {}", self.percent, self.count)
    }

    fn vec_for_color(&self) -> Option<(Vec<u8>, i32)> {
        None
    }

    fn get_tags_sum(&self) -> Option<(Tags, i32)> {
        None
    }

    fn score(&self) -> f32 {
        self.count as f32
    }

    fn mem(&self) -> usize {
        mem::size_of_val(&*self)
    }
}

impl SummaryData<f32> for f32 {
    fn new(data: f32) -> Self {
        data
    }

    fn print(&self, _: &BiMap<String, u8>) -> String {
        format!("relative count group 1: {}", self)
    }

    fn vec_for_color(&self) -> Option<(Vec<u8>, i32)> {
        None
    }

    fn score(&self) -> f32 {
        0.
    }
    
    fn get_tags_sum(&self) -> Option<(Tags, i32)> {
        None
    }

    fn mem(&self) -> usize {
        mem::size_of_val(&*self)
    }
}

/// Implement this trait to control how multiple observations of a kmer
/// are carried forward into a DeBruijn graph.
pub trait KmerSummarizer<DI, DO: SummaryData<SD>, SD> {
    /// The input `items` is an iterator over kmer observations. Input observation
    /// is a tuple of (kmer, extensions, data). The summarize function inspects the
    /// data and returns a tuple indicating:
    /// * whether this kmer passes the filtering criteria (e.g. is there a sufficient number of observation)
    /// * the accumulated Exts of the kmer
    /// * a summary data object of type `DO` that will be used as a color annotation in the DeBruijn graph.
    
    fn new(min_kmer_obs: usize, markers: Marker) -> Self;
    fn summarize<K: Kmer, F: Iterator<Item = (K, Exts, DI)>>(&self, items: F, significant: Option<u32>) -> (bool, Exts, DO);
}

/// A simple KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times. The metadata returned about a Kmer
/// is the number of times it was observed, capped at 2^16.
pub struct CountFilter<D> {
    min_kmer_obs: usize,
    phantom: PhantomData<D>
}

impl<D> KmerSummarizer<D, u32, u32> for CountFilter<D> {
    fn new(min_kmer_obs: usize, _: Marker) -> Self {
        CountFilter {
            min_kmer_obs,
            phantom: PhantomData,
        }
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, D)>>(&self, items: F, significant: Option<u32>) -> (bool, Exts, u32) {
        let mut all_exts = Exts::empty();
        let mut count = 0u32;
        for (_, exts, _) in items {
            count = count.saturating_add(1);
            all_exts = all_exts.add(exts);
        }

        let count = match significant {
            Some(digits) => round_digits(count, digits),
            None => count  
        };

        (count as usize >= self.min_kmer_obs, all_exts, count)
    }
}

/// A simple KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times. The metadata returned about a Kmer
/// is a vector of the unique data values observed for that kmer.
pub struct CountFilterSet<D> {
    min_kmer_obs: usize,
    phantom: PhantomData<D>,
}

impl<D: Ord + Debug> KmerSummarizer<D, Vec<D>, Vec<D>> for CountFilterSet<D> {
    fn new(min_kmer_obs: usize, _: Marker) -> Self {
        CountFilterSet {
            min_kmer_obs,
            phantom: PhantomData,
        }
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, D)>>(&self, items: F, _: Option<u32>) -> (bool, Exts, Vec<D>) {
        let mut all_exts = Exts::empty();

        let mut out_data: Vec<D> = Vec::with_capacity(items.size_hint().0);

        let mut nobs = 0i32;
        for (_, exts, d) in items {
            out_data.push(d);
            all_exts = all_exts.add(exts);
            nobs += 1;
        }

        out_data.sort();
        out_data.dedup();
        out_data.shrink_to_fit();
        
        (nobs as usize >= self.min_kmer_obs, all_exts, out_data)
        
    }
}

/// A simple KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times. The metadata returned about a Kmer
/// is a vector of the unique data values observed for that kmer.
pub struct CountFilterComb {
    min_kmer_obs: usize,
    phantom: PhantomData<u8>,
}

impl KmerSummarizer<u8, TagsSumData, (Tags, i32)> for CountFilterComb {
    fn new(min_kmer_obs: usize, _: Marker) -> Self {
        CountFilterComb {
            min_kmer_obs,
            phantom: PhantomData,
        }
    }

    fn summarize<K: Kmer, F: Iterator<Item = (K, Exts, u8)>>(&self, items: F, _: Option<u32>) -> (bool, Exts, TagsSumData) {
        let mut all_exts = Exts::empty();

        let mut out_data: Vec<u8> = Vec::with_capacity(items.size_hint().0);

        let mut nobs = 0i32;
        for (_, exts, d) in items {
            out_data.push(d); // uses a shit ton of heap memory
            all_exts = all_exts.add(exts);
            nobs += 1;
        }

        out_data.sort();
        out_data.dedup();

        (nobs as usize >= self.min_kmer_obs, all_exts, TagsSumData::new((Tags::from_u8_vec(out_data), nobs)))
    }
}


/// A simple KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times. The metadata returned about a Kmer
/// is a vector of the unique data values observed for that kmer.
pub struct CountFilterStats {
    min_kmer_obs: usize,
    phantom: PhantomData<u8>,
}

impl KmerSummarizer<u8, TagsCountData, (Tags, Box<[u32]>, i32)> for CountFilterStats {
    fn new(min_kmer_obs: usize, _: Marker) -> Self {
        CountFilterStats {
            min_kmer_obs,
            phantom: PhantomData,
        }
    }

    fn summarize<K: Kmer, F: Iterator<Item = (K, Exts, u8)>>(&self, items: F, _: Option<u32>) -> (bool, Exts, TagsCountData) {
        let mut all_exts = Exts::empty();

        let mut out_data: Vec<u8> = Vec::with_capacity(items.size_hint().0);

        let mut nobs = 0i32;
        for (_, exts, d) in items {
            out_data.push(d); 
            all_exts = all_exts.add(exts);
            nobs += 1;
        }

        out_data.sort();

        let mut tag_counter = 1;
        let mut tag_counts: Vec<u32> = Vec::new();

        // count the occurences of the labels
        for i in 1..out_data.len() {
            if out_data[i] == out_data[i-1] {
                tag_counter += 1;
            } else {
                tag_counts.push(tag_counter.clone());
                tag_counter = 1;
            }
        }
        tag_counts.push(tag_counter);

        out_data.dedup();

        let tag_counts: Box<[u32]> = tag_counts.into();

        (nobs as usize >= self.min_kmer_obs, all_exts, TagsCountData::new((Tags::from_u8_vec(out_data), tag_counts, nobs))) 
    }
}


/// A simple KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times. The metadata returned about a Kmer
/// is a vector of the unique data values observed for that kmer.
pub struct CountsFilterStats {
    min_kmer_obs: usize,
    phantom: PhantomData<u8>,
}

impl KmerSummarizer<u8, TagsCountsData, (Tags, Box<[u32]>)> for CountsFilterStats {
    fn new(min_kmer_obs: usize, _: Marker) -> Self {
        CountsFilterStats {
            min_kmer_obs,
            phantom: PhantomData,
        }
    }

    fn summarize<K: Kmer, F: Iterator<Item = (K, Exts, u8)>>(&self, items: F, _: Option<u32>) -> (bool, Exts, TagsCountsData) {
        let mut all_exts = Exts::empty();

        let mut out_data: Vec<u8> = Vec::with_capacity(items.size_hint().0);

        let mut nobs = 0i32;
        for (_, exts, d) in items {
            out_data.push(d); 
            all_exts = all_exts.add(exts);
            nobs += 1;
        }

        out_data.sort();

        let mut tag_counter = 1;
        let mut tag_counts: Vec<u32> = Vec::new();

        // count the occurences of the labels
        for i in 1..out_data.len() {
            if out_data[i] == out_data[i-1] {
                tag_counter += 1;
            } else {
                tag_counts.push(tag_counter.clone());
                tag_counter = 1;
            }
        }
        tag_counts.push(tag_counter);

        out_data.dedup();

        let tag_counts: Box<[u32]> = tag_counts.into();

        (nobs as usize >= self.min_kmer_obs, all_exts, TagsCountsData::new((Tags::from_u8_vec(out_data), tag_counts))) 
    }
}

/// A simple KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times. The metadata returned about a Kmer
/// is a vector of the unique data values observed for that kmer.
pub struct CountsFilterGroups {
    min_kmer_obs: usize,
    marker: Marker,
    phantom: PhantomData<u8>,
}

impl KmerSummarizer<u8, GroupCountData, (u32, u32)> for CountsFilterGroups {
    fn new(min_kmer_obs: usize, marker: Marker) -> Self {
        CountsFilterGroups {
            min_kmer_obs,
            marker,
            phantom: PhantomData,
        }
    }

    fn summarize<K: Kmer, F: Iterator<Item = (K, Exts, u8)>>(&self, items: F, significant: Option<u32>) -> (bool, Exts, GroupCountData) {
        let mut all_exts = Exts::empty();
        let mut count1 = 0;
        let mut count2 = 0;

        let mut nobs = 0i32;
        for (_, exts, d) in items {
            let tag = (2 as M).pow(d as u32);
            let group1 = ((self.marker.marker0 & tag) > 1) as u32;
            let group2 = ((self.marker.marker1 & tag) > 1) as u32;

            if (group1 + group2) != 1 { panic!("should not happen\n tag: {:#066b}\n m1: {:#066b}\n m2: {:#066b}\n g1: {:#066b}\n g2: {:#066b}", tag, self.marker.marker0, self.marker.marker1, group1, group2)}
            count1 += group1;
            count2 += group2;
            nobs += 1;
            all_exts = all_exts.add(exts);
        }

        let counts = match significant {
            Some(digits) => (round_digits(count1, digits), round_digits(count2, digits)),
            None => (count1, count2)
        };

        assert_eq!((count1 + count2),nobs as u32);
        (nobs as usize >= self.min_kmer_obs, all_exts, GroupCountData::new(counts))
    }
}

/// A simple KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times. The metadata returned about a Kmer
/// is a vector of the unique data values observed for that kmer.
pub struct CountsFilterRel{
    min_kmer_obs: usize,
    marker: Marker,
    phantom: PhantomData<u8>,
}

impl KmerSummarizer<u8, RelCountData, (u32, u32)> for CountsFilterRel {
    fn new(min_kmer_obs: usize, marker: Marker) -> Self {
        CountsFilterRel {
            min_kmer_obs,
            marker,
            phantom: PhantomData,
        }
    }

    fn summarize<K: Kmer, F: Iterator<Item = (K, Exts, u8)>>(&self, items: F, significant: Option<u32>) -> (bool, Exts, RelCountData) {
        let mut all_exts = Exts::empty();
        let mut count1 = 0;
        let mut count2 = 0;

        let mut nobs = 0i32;
        for (_, exts, d) in items {
            let tag = (2 as M).pow(d as u32);
            let group1 = ((self.marker.marker0 & tag) > 1) as u32;
            let group2 = ((self.marker.marker1 & tag) > 1) as u32;

            if (group1 + group2) != 1 { panic!("should not happen\n tag: {:#066b}\n m1: {:#066b}\n m2: {:#066b}\n g1: {:#066b}\n g2: {:#066b}", tag, self.marker.marker0, self.marker.marker1, group1, group2)}
            count1 += group1;
            count2 += group2;
            all_exts = all_exts.add(exts);
            nobs += 1;
        }

        assert_eq!(count1 + count2, nobs as u32);


        let percent = (count1 as f64 / nobs as f64 * 100.) as u32;
        let count = match significant {
            Some(digits) => round_digits(nobs as u32, digits),
            None => nobs as u32
        };

        (nobs as usize >= self.min_kmer_obs, all_exts, RelCountData::new((percent, count))) 
    }
}

pub struct CountsFilterMaj {
    min_kmer_obs: usize,
    marker: Marker,
    phantom: PhantomData<u8>,
}


impl KmerSummarizer<u8, TagsCountsData, (Tags, Box<[u32]>)> for CountsFilterMaj {
    fn new(min_kmer_obs: usize, marker: Marker) -> Self {
        CountsFilterMaj {
            min_kmer_obs,
            marker,
            phantom: PhantomData,
        }
    }

    fn summarize<K: Kmer, F: Iterator<Item = (K, Exts, u8)>>(&self, items: F, _: Option<u32>) -> (bool, Exts, TagsCountsData) {
        let mut all_exts = Exts::empty();

        let mut out_data: Vec<u8> = Vec::with_capacity(items.size_hint().0);

        let mut nobs = 0i32;
        for (_, exts, d) in items {
            out_data.push(d); 
            all_exts = all_exts.add(exts);
            nobs += 1;
        }

        out_data.sort();

        let mut tag_counter = 1;
        let mut tag_counts: Vec<u32> = Vec::new();

        // count the occurences of the labels
        for i in 1..out_data.len() {
            if out_data[i] == out_data[i-1] {
                tag_counter += 1;
            } else {
                tag_counts.push(tag_counter.clone());
                tag_counter = 1;
            }
        }
        tag_counts.push(tag_counter);
        out_data.dedup(); 

        let tags = Tags::from_u8_vec(out_data.clone());

        

        let tag_counts: Box<[u32]> = tag_counts.into();

        (nobs as usize >= self.min_kmer_obs, all_exts, TagsCountsData::new((Tags::from_u8_vec(out_data), tag_counts))) 
    }
}

#[cfg(test)]
mod test {
    use std::fmt::Debug;

    use bimap::BiMap;
    use boomphf::hashmap::BoomHashMap2;
    use rand::Rng;

    use crate::{compression::{ compress_kmers_with_hash, CompressionSpec, ScmapCompress}, filter::filter_kmers, graph::DebruijnGraph, kmer::Kmer12, reads::Reads, summarizer::{CountFilter, CountsFilterGroups, CountsFilterRel, CountsFilterStats, GroupCountData, KmerSummarizer, Marker, RelCountData, SummaryData, TagsCountsData}, Exts, Kmer};

    #[test]
    fn test_summarizers() {

        // generate three reads
        let mut reads = Reads::new();

        let dnas = ["CGATCGAGCTACTGCGACGGACGATTTTTCGAGCGGCGATTTCTCGAGGCGAGCGTCAGC".as_bytes(),
            "CGATCGAGCTACTGCGACGGACGATGACTAGCTAGCTTTTCTCGAGGCGAGCGTCAGC".as_bytes(),
            "ACTGCGACGGACGATGACTAGCTAGCTTTTCTCGAGGCGAGCGTCAGCACGATGCTAGCTGACTAGC".as_bytes(),
            "CGATCGAGCTACTGCGACGGACGATTTTTCGAGCGGCGATTTCTCGAGGCGAGCGTCAGCGGACTAGCGAG".as_bytes(),
            "ACGGACGATTTTTCGAGCGGCGATTTCTCGAGGCGAGCGTCAGC".as_bytes(),
        ];

        let mut rng = rand::thread_rng();
        for dna in dnas {
            reads.add_from_bytes(dna, Exts::empty(), rng.gen_range(1, 4));
        }


        /* let repeats = 2000;
        let read_len = 150;

        for _i in 0..repeats {
            let dna =  random_dna((1.5*read_len as f32) as usize);
            let dna_string1 = DnaString::from_bytes(&dna[0..read_len + read_len/2]);
            let dna_string2 = DnaString::from_bytes(&dna[read_len + read_len/2..(1.5*read_len as f32) as usize]);
            reads.add_read(dna_string1.clone(), Exts::empty(), 1u8);
            if rand::random::<bool>() { reads.add_read(dna_string2, Exts::empty(), 2u8) };
        }

        for _i in 0..repeats {
            let dna =  random_dna((1.5*read_len as f32) as usize);
            let dna_string1 = DnaString::from_bytes(&dna[0..read_len + read_len/2]);
            let dna_string2 = DnaString::from_bytes(&dna[read_len + read_len/2..(1.5*read_len as f32) as usize]);            
            reads.add_read(dna_string1.clone(), Exts::empty(), 2u8);
            if rand::random::<bool>() { reads.add_read(dna_string2, Exts::empty(), 3u8) };
        }

        for _i in 0..repeats {
            let dna_string = DnaString::from_bytes(&random_dna(read_len));
            reads.add_read(dna_string, Exts::empty(), 3u8);
        }
 */
        // markers: 
        let markers = Marker { marker0: 2, marker1: 12, count0: 1, count1: 2 };
        // 0010 => 2
        // 1100 => 12

        // tag translator
        let mut tag_translator: bimap::BiHashMap<String, u8> = BiMap::new();
        tag_translator.insert("sample 1".to_string(), 1);
        tag_translator.insert("sample 2".to_string(), 2);
        tag_translator.insert("sample 3".to_string(), 3);

        let significant= Some(4);
        let min_kmer_obs = 1;
        type K = Kmer12;

        let pr = false;

        if pr { println!("{}", reads) };


        //construct and compress graph with CountFilter
        let summarizer= CountFilter::new(min_kmer_obs, markers);
        let spec: ScmapCompress<u32> = ScmapCompress::new();
        let graph: DebruijnGraph<K, u32> = test_summarizer(&reads, summarizer, spec, significant);

        graph.to_gfa_with_tags("test_cf.gfa",|n| n.data().print(&tag_translator)).unwrap();

        println!("cf graph size: {}", graph.len());


        // construct and compress graph with CountsFilterGroups
        let summarizer= CountsFilterGroups::new(min_kmer_obs, markers);
        let spec: ScmapCompress<GroupCountData> = ScmapCompress::new();
        let graph: DebruijnGraph<K, GroupCountData> = test_summarizer(&reads, summarizer, spec, significant);

        //graph.to_gfa_with_tags("test_csfg.gfa",|n| n.data().print(&tag_translator)).unwrap();

        println!("csfg graph size: {}", graph.len());


        // construct and compress graph with CountsFilterRel
        let summarizer= CountsFilterRel::new(min_kmer_obs, markers);
        let spec: ScmapCompress<RelCountData> = ScmapCompress::new();
        let graph: DebruijnGraph<K, RelCountData> = test_summarizer(&reads, summarizer, spec, significant);

        //graph.to_gfa_with_tags("test_csfr.gfa",|n| n.data().print(&tag_translator)).unwrap();

        println!("csfr graph size: {}", graph.len());


        // construct and compress graph with CountsFilterStats
        let summarizer= CountsFilterStats::new(min_kmer_obs, markers);
        let spec: ScmapCompress<TagsCountsData> = ScmapCompress::new();
        let graph: DebruijnGraph<K, TagsCountsData> = test_summarizer(&reads, summarizer, spec, significant);

        graph.to_gfa_with_tags("test_csfs.gfa",|n| n.data().print(&tag_translator)).unwrap();

        println!("csfs graph size: {}", graph.len());



        // same but with less significant digits
        let significant= Some(1);

        //construct and compress graph with CountFilter
        let summarizer= CountFilter::new(min_kmer_obs, markers);
        let spec: ScmapCompress<u32> = ScmapCompress::new();
        let graph: DebruijnGraph<K, u32> = test_summarizer(&reads, summarizer, spec, significant);

        graph.to_gfa_with_tags("test_cf-1.gfa",|n| n.data().print(&tag_translator)).unwrap();

        println!("cf graph size: {}", graph.len());


        // construct and compress graph with CountsFilterGroups
        let summarizer= CountsFilterGroups::new(min_kmer_obs, markers);
        let spec: ScmapCompress<GroupCountData> = ScmapCompress::new();
        let graph: DebruijnGraph<K, GroupCountData> = test_summarizer(&reads, summarizer, spec, significant);

        //graph.to_gfa_with_tags("test_csfg.gfa",|n| n.data().print(&tag_translator)).unwrap();

        println!("csfg graph size: {}", graph.len());


        // construct and compress graph with CountsFilterRel
        let summarizer= CountsFilterRel::new(min_kmer_obs, markers);
        let spec: ScmapCompress<RelCountData> = ScmapCompress::new();
        let graph: DebruijnGraph<K, RelCountData> = test_summarizer(&reads, summarizer, spec, significant);

        //graph.to_gfa_with_tags("test_csfr.gfa",|n| n.data().print(&tag_translator)).unwrap();

        println!("csfr graph size: {}", graph.len());


        // construct and compress graph with CountsFilterStats
        let summarizer= CountsFilterStats::new(min_kmer_obs, markers);
        let spec: ScmapCompress<TagsCountsData> = ScmapCompress::new();
        let graph: DebruijnGraph<K, TagsCountsData> = test_summarizer(&reads, summarizer, spec, significant);

        //graph.to_gfa_with_tags("test_csfs.gfa",|n| n.data().print(&tag_translator)).unwrap();

        println!("csfs graph size: {}", graph.len());


    }


    fn test_summarizer<K: Kmer + Send + Sync, T: KmerSummarizer<u8, D, SD>, CS: CompressionSpec<D> + Send + Sync, D: SummaryData<SD> + Debug + PartialEq + Clone + Send + Sync, SD>(reads: &Reads<u8>, summarizer: T, spec: CS, sig: Option<u32>) -> DebruijnGraph<K, D> {
        // construct and compress graph with CountsFilterStats
        //let summarizer= CountsFilterStats::new(min_kmer_obs, markers);

        let memory = 1;
        let pr = false;

        let (k_mers, _): (BoomHashMap2<K, Exts, D>, _)  = filter_kmers(&reads,
            &Box::new(summarizer),
            false, 
            false, 
            memory, 
            false, 
            false, 
            sig
        );
        if pr { println!("kmers CountsFilterStats: {:?}", k_mers) };

        let graph = compress_kmers_with_hash(false, &spec, &k_mers, false, false, false);
        if pr { println!("graph CountsFilterStats: {:?}\n", graph) };

        //graph.finish().to_gfa_with_tags("test_csfs.gfa",|n| n.data().print(&tt)).unwrap();
        graph.finish()
    }
}


