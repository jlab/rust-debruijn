use bimap::BiMap;
use crate::{Exts, Kmer, Tags};
use std::{fmt::Debug, marker::PhantomData, mem};

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
    fn print(&self, tag_translator: &BiMap<&str, u8>) -> String;
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

    fn print(&self, _: &BiMap<&str, u8>) -> String {
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

    fn print(&self, _: &BiMap<&str, u8>) -> String {
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

#[derive(Debug, Clone, PartialEq)]
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

    fn print(&self, tag_translator: &BiMap<&str, u8>) -> String {
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

#[derive(Debug, Clone, PartialEq)]
pub struct TagsCountData {
    tags: Tags,
    sum: i32,
    counts: Box<[u32]>,
}

impl SummaryData<(Tags, Box<[u32]>, i32)> for TagsCountData {
    fn new(data: (Tags, Box<[u32]>, i32)) -> Self {
        TagsCountData { tags: data.0, counts: data.1, sum: data.2 }
    }

    fn print(&self, tag_translator: &BiMap<&str, u8>) -> String {
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

    fn print(&self, tag_translator: &BiMap<&str, u8>) -> String {
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
    
    fn print(&self, _: &BiMap<&str, u8>) -> String {
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
pub struct RelCountData {
    percent: u32,
    count: u32
}

impl SummaryData<(u32, u32)> for RelCountData {
    fn new(data: (u32, u32)) -> Self {
        RelCountData { percent: data.0, count: data.1 }
    }
    
    fn print(&self, _: &BiMap<&str, u8>) -> String {
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

    fn print(&self, _: &BiMap<&str, u8>) -> String {
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
    
    fn new(min_kmer_obs: usize, markers: (u64, u64)) -> Self;
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
    fn new(min_kmer_obs: usize, _: (u64, u64)) -> Self {
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
    fn new(min_kmer_obs: usize, _: (u64, u64)) -> Self {
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
    fn new(min_kmer_obs: usize, _: (u64, u64)) -> Self {
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
    fn new(min_kmer_obs: usize, _: (u64, u64)) -> Self {
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
    fn new(min_kmer_obs: usize, _: (u64, u64)) -> Self {
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
    marker1: u64,
    marker2: u64,
    phantom: PhantomData<u8>,
}

impl KmerSummarizer<u8, GroupCountData, (u32, u32)> for CountsFilterGroups {
    fn new(min_kmer_obs: usize, (marker1, marker2): (u64, u64)) -> Self {
        CountsFilterGroups {
            min_kmer_obs,
            marker1,
            marker2,
            phantom: PhantomData,
        }
    }

    fn summarize<K: Kmer, F: Iterator<Item = (K, Exts, u8)>>(&self, items: F, significant: Option<u32>) -> (bool, Exts, GroupCountData) {
        let mut all_exts = Exts::empty();
        let mut count1 = 0;
        let mut count2 = 0;

        let mut nobs = 0i32;
        for (_, exts, d) in items {
            let tag = 2u64.pow(d as u32);
            let group1 = self.marker1 & tag;
            let group2 = self.marker2 & tag;
            if (group1 + group2) > 1 { panic!("should not happen")}
            count1 += group1;
            count2 += group2;
            all_exts = all_exts.add(exts);
            nobs += 1;
        }

        let counts = match significant {
            Some(digits) => (round_digits(count1 as u32, digits), round_digits(count2 as u32, digits)),
            None => (count1 as u32, count2 as u32)
        };

        (nobs as usize >= self.min_kmer_obs, all_exts, GroupCountData::new(counts))
    }
}

/// A simple KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times. The metadata returned about a Kmer
/// is a vector of the unique data values observed for that kmer.
pub struct CountsFilterRel{
    min_kmer_obs: usize,
    marker1: u64,
    marker2: u64,
    phantom: PhantomData<u8>,
}

impl KmerSummarizer<u8, RelCountData, (u32, u32)> for CountsFilterRel {
    fn new(min_kmer_obs: usize, (marker1, marker2): (u64, u64)) -> Self {
        CountsFilterRel {
            min_kmer_obs,
            marker1,
            marker2,
            phantom: PhantomData,
        }
    }

    fn summarize<K: Kmer, F: Iterator<Item = (K, Exts, u8)>>(&self, items: F, significant: Option<u32>) -> (bool, Exts, RelCountData) {
        let mut all_exts = Exts::empty();
        let mut count1 = 0;
        let mut count2 = 0;

        let mut nobs = 0i32;
        for (_, exts, d) in items {
            let tag = 2u64.pow(d as u32);
            let group1 = self.marker1 & tag;
            let group2 = self.marker2 & tag;
            if (group1 + group2) > 1 { panic!("should not happen")}
            count1 += group1;
            count2 += group2;
            all_exts = all_exts.add(exts);
            nobs += 1;
        }

        let percent = (count1 as f64 / (count1 + count2) as f64 * 100.) as u32;
        let count = match significant {
            Some(digits) => round_digits((count1 + count2) as u32, digits),
            None => (count1 + count2) as u32
        };

        (nobs as usize >= self.min_kmer_obs, all_exts, RelCountData::new((percent, count))) 
    }
}

