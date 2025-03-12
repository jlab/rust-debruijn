use bimap::BiMap;
use clap::ValueEnum;
use serde::{Deserialize, Serialize};
use statrs::distribution::{ContinuousCDF, Normal, StudentsT};
use crate::{EdgeMult, Exts, Kmer, Tags};
use std::{cmp::min_by, error::Error, fmt::{Debug, Display}, marker::PhantomData, mem};

#[cfg(not(feature = "sample128"))]
pub type M = u64;

#[cfg(feature = "sample128")]
pub type M = u128;

#[derive(Debug, PartialEq)]
struct NotEnoughSamplesError {}

impl Error for NotEnoughSamplesError {}

impl Display for NotEnoughSamplesError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "not enough samples were supplied to perform a statistical test")
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct SummaryConfig {
    min_kmer_obs: usize,
    significant: Option<u32>,
    group_frac: GroupFrac,
    frac_cutoff: f32,
    sample_info: SampleInfo,
    max_p: Option<f32>,
    stat_test: StatTest,
    stat_test_changed: bool,
}

impl SummaryConfig {
    /// make a new `SummaryConfig`
    /// arguments: 
    /// * `min_kmer_obs`: minimum number of times a k-mer has to be observed in the reads to be valid
    /// * `significant`: some summaries round numbers to a certain number of digits
    /// * `third`
    pub fn new(min_kmer_obs: usize, significant: Option<u32>, group_frac: GroupFrac, frac_cutoff: f32, sample_info: SampleInfo, max_p: Option<f32>, stat_test: StatTest) -> Self {
        SummaryConfig { min_kmer_obs, significant, group_frac, frac_cutoff, sample_info, max_p, stat_test, stat_test_changed: false }
    }

    pub fn empty() -> Self {
        SummaryConfig { min_kmer_obs: 0, significant: None, group_frac: GroupFrac::None, frac_cutoff: 0., sample_info: SampleInfo::empty(), max_p: None, stat_test: StatTest::TTest, stat_test_changed: false }
    }

    pub fn set_min_kmer_obs(&mut self, min_kmer_obs: usize) {
        self.min_kmer_obs = min_kmer_obs;
    }

    pub fn set_group_frac(&mut self, group_frac: GroupFrac, frac_cutoff: f32) {
        self.group_frac = group_frac;
        self.frac_cutoff = frac_cutoff;
    }

    pub fn  set_max_p(&mut self, max_p: Option<f32>) {
        self.max_p = max_p;
    }

    pub fn set_stat_test(&mut self, stat_test: StatTest) {
        if stat_test != self.stat_test { self.stat_test_changed = true }
        self.stat_test = stat_test;
    }
}

#[derive(Copy, Clone, PartialEq, PartialOrd, ValueEnum, Debug, Serialize, Deserialize)]
pub enum GroupFrac {
    None, 
    One, 
    Both,
}

impl std::fmt::Display for GroupFrac {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self {
            GroupFrac::None => write!(f, "none"),
            GroupFrac::One => write!(f, "one"),
            GroupFrac::Both => write!(f, "both")            
        }
    }
    
}

#[derive(Copy, Clone, PartialEq, PartialOrd, ValueEnum, Debug, Serialize, Deserialize)]
pub enum StatTest {
    TTest,
    UTest,
}

impl std::fmt::Display for StatTest {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::TTest => write!(f, "t-test"),
            Self::UTest => write!(f, "u-test"),
            
        }
    }
    
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct SampleInfo {
    marker0: M,
    marker1: M,
    count0: u8,
    count1: u8,
    sample_kmers: Vec<u64>,
}

impl SampleInfo {
    pub fn new(marker0: M, marker1: M, count0: u8, count1: u8, sample_kmers: Vec<u64>) -> Self {
        SampleInfo { marker0, marker1, count0, count1, sample_kmers }
    }

    pub fn empty() -> Self {
        SampleInfo { marker0: 0, marker1: 0, count0: 0, count1: 0, sample_kmers: Vec::new() }
    }

    pub fn get_marker(&self) -> (M, M) {
        (self.marker0, self.marker1)
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

fn valid(tags: Tags, nobs: i32, config: &SummaryConfig) -> bool {
    match config.group_frac {
        GroupFrac::None => nobs as usize >= config.min_kmer_obs,
        GroupFrac::Both => {
            // get amount of labels in tags from each group
            let dist0= tags.bit_and_dist(config.sample_info.marker0);
            let dist1= tags.bit_and_dist(config.sample_info.marker1);
    
            assert_eq!(dist0 + dist1, tags.to_u8_vec().len());
    
            // valid if:
            // - n obs >= min obs AND
            // - observed in at least one third of samples in both groups
            (nobs as usize >= config.min_kmer_obs)
                && (dist0 as f32 / config.sample_info.count0 as f32 >= config.frac_cutoff) 
                && (dist1 as f32 / config.sample_info.count1 as f32 >= config.frac_cutoff)
        },
        GroupFrac::One => {
            // get amount of labels in tags from each group
            let dist0= tags.bit_and_dist(config.sample_info.marker0);
            let dist1= tags.bit_and_dist(config.sample_info.marker1);

    
            assert_eq!(dist0 + dist1, tags.to_u8_vec().len());
    
            // valid if:
            // - n obs >= min obs AND
            // - observed in at least one third of samples in one group
            (nobs as usize >= config.min_kmer_obs)
                && ((dist0 as f32 / config.sample_info.count0 as f32 >= config.frac_cutoff) 
                    | (dist1 as f32 / config.sample_info.count1 as f32 >= config.frac_cutoff))
        }
    }
}

fn t_test(out_data: &Vec<u8>, tag_counts: &Vec<u32>, sample_info: &SampleInfo) -> Result<f32, NotEnoughSamplesError> {
    let n0 = sample_info.count0 as f64;
    let n1 = sample_info.count1 as f64;

    if (n0 < 2.) | (n1 < 2.) { return Err(NotEnoughSamplesError {})}

    let mut counts_g0 = Vec::new();
    let mut counts_g1 = Vec::new();

    let (m0, m1) = sample_info.get_marker();

    for (label, count) in out_data.iter().zip(tag_counts) {
        let bin_rep = (2 as M).pow(*label as u32);
        let norm = *count as f64 / sample_info.sample_kmers[*label as usize] as f64;
        if (m0 & bin_rep) > 0 { counts_g0.push(norm); }
        if (m1 & bin_rep) > 0 { counts_g1.push(norm); }
    }

    let mean0 = counts_g0.iter().sum::<f64>() as f64 / n0;
    let mean1 = counts_g1.iter().sum::<f64>() as f64 / n1;

    let df = n0 + n1 - 2.;

    let var0 = (counts_g0.iter().map(|count| (*count - mean0).powf(2.)).sum::<f64>() + (n0 - counts_g0.len() as f64) * mean0.powf(2.)) / (n0 - 1.);
    let var1 = (counts_g1.iter().map(|count| (*count - mean1).powf(2.)).sum::<f64>() + (n1 - counts_g1.len() as f64) * mean1.powf(2.)) / (n1 - 1.);

    let s = ((1./n0 + 1./n1) * ((n0 - 1.) * var0 + (n1 - 1.) * var1) / df).sqrt();

    let t = (mean0 - mean1) / s;

    let t_dist = StudentsT::new(0.0, 1.0, df).expect("error creating student dist: check if you have enough samples (at least 3)");

    let p_value = 2. * (1. - t_dist.cdf(t.abs())) as f32;

    Ok(p_value)
}

fn u_test(out_data: &Vec<u8>, tag_counts: &Vec<u32>, sample_info: &SampleInfo) -> Result<f32, NotEnoughSamplesError> {

    let n0 = sample_info.count0 as f64;
    let n1 = sample_info.count1 as f64;

    if (n0 < 2.) | (n1 < 2.) { return Err(NotEnoughSamplesError {})}


    let mut counts_g0 = Vec::new();
    let mut counts_g1 = Vec::new();

    let (m0, m1) = sample_info.get_marker();

    for (label, count) in out_data.iter().zip(tag_counts) {
        let bin_rep = (2 as M).pow(*label as u32);
        let norm = *count as f64 / sample_info.sample_kmers[*label as usize] as f64;
        if (m0 & bin_rep) > 0 { counts_g0.push(norm); }
        if (m1 & bin_rep) > 0 { counts_g1.push(norm); }
    }

    let n = n0 + n1;
    let m = (n0 * n1) / 2.;

    // TODO check if more efficient way possible

    let mut all_counts = vec![(0u8, 0f64); n0 as usize - counts_g0.len()];
    all_counts.append(&mut vec![(1, 0f64); n1 as usize - counts_g1.len()]);

    all_counts.append(&mut counts_g0.iter().map(|elt| (0, *elt)).collect());
    all_counts.append(&mut counts_g1.iter().map(|elt| (1, *elt)).collect());

    all_counts.sort_by(|(_, x), (_, y)| x.total_cmp(y));

    let chunked = all_counts
        .chunk_by(|(_, x), (_, y)| x == y);

    let mut ranks = Vec::new();
    let mut tie_factor = 0.;

    for chunk in chunked {
        let rank = (chunk.len() + 1) as f64 / 2. + ranks.len() as f64;
        ranks.append(&mut chunk.iter().map(|(g, _)| (*g, rank)).collect());
        if chunk.len() > 0 {
            tie_factor += (chunk.len().pow(3) - chunk.len()) as f64;
        }
    }

    let mut rank_sum0 = 0.;
    let mut rank_sum1 = 0.;
    ranks.iter().for_each(|(group, rank)| match group { 
        0 => rank_sum0 += rank, 
        1 => rank_sum1 += rank, 
        _ => panic!("should not happen"),
    });

    let u0 = rank_sum0 - n0 * (n0 + 1.) / 2.;
    let u1 = rank_sum1 - n1 * (n1 + 1.) / 2.;
    
    let u = min_by(u0, u1, |a, b| a.total_cmp(&b));

    let s = ((n0 * n1 * (n + 1.) / 12.) - (n0 * n1 * tie_factor / (12. * n * (n - 1.)))).sqrt();
    //let s = ((n0 * n1 / 12.) * ((n + 1.) - (tie_factor / n * (n - 1.)))).sqrt();
    // supposedly same term but returns NaN???

    let z = (u - m) / s;

    let dist = Normal::standard();
    let p_value =  2. * (1. - dist.cdf(z.abs())) as f32;

    Ok(p_value)
}

// calculates the log2 of the log change of the two groups
fn log2_fold_change(tags: Tags, counts: Vec<u32>, sample_info: &SampleInfo) -> f32 {
    let mut norm_count_g0 = 0.;
    let mut norm_count_g1 = 0.;

    let (m0, m1) = sample_info.get_marker();

    for (label, count) in tags.to_u8_vec().iter().zip(&counts) {
        let bin_rep = (2 as M).pow(*label as u32);
        // normalize with number of k-mers in the sample
        let norm = *count as f64 / sample_info.sample_kmers[*label as usize] as f64;
        if (m0 & bin_rep) > 0 { norm_count_g0 += norm; }
        if (m1 & bin_rep) > 0 { norm_count_g1 += norm; }
    }

    // normalize with the number of samples in the group
    norm_count_g0 /= sample_info.count0 as f64;
    norm_count_g1 /= sample_info.count1 as f64;

    (norm_count_g0 / norm_count_g1).log2() as f32
}

/// Trait for summarizing k-mers
pub trait SummaryData<DI, DO> {
    /// Make a new `SummaryData<DO>`
    fn new(data: DO) -> Self;
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
    /// If the `SummaryData` contains sufficient information, return the number of observations
    fn count(&self) -> Option<usize>;
    /// If the `SummaryData` contains sufficient information, return the p-value
    fn p_value(&self, config: &SummaryConfig) -> Option<f32>;
    /// If the `SummaryData` contains sufficient information, return the log2(fold change)
    fn fold_change(&self, config: &SummaryConfig) -> Option<f32>;
    /// If the `SummaryData` contains sufficient information, return the number of samples the sequence was observed in
    fn sample_count(&self) -> Option<usize>;
    /// check if node is valid according to: min kmer obs, third, p-value
    fn valid(&self, config: &SummaryConfig) -> bool;
    /// summarize k-mers
    fn summarize<K, F: Iterator<Item = (K, Exts, DI)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self);
}

/// for CountFilter
impl<DI> SummaryData<DI, u32> for u32 {
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

    fn count(&self) -> Option<usize> {
        Some(*self as usize)
    }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> {
        None
    }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> {
        None
    }

    fn sample_count(&self) -> Option<usize> {
        None
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        *self >= config.min_kmer_obs as u32
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, DI)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let mut all_exts = Exts::empty();
        let mut count = 0u32;
        for (_, exts, _) in items {
            count = count.saturating_add(1);
            all_exts = all_exts.add(exts);
        }

        let count = match config.significant {
            Some(digits) => round_digits(count, digits),
            None => count  
        };

        (count as usize >= config.min_kmer_obs, all_exts, count)
    }

}

/// for CountFilterSet
impl<DI: Debug + Ord> SummaryData<DI, Vec<DI>> for Vec<DI> {
    fn new(data: Vec<DI>) -> Self {
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

    fn count(&self) -> Option<usize> {
        None
    }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> {
        None
    }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> {
        None
    }

    fn sample_count(&self) -> Option<usize> {
        Some(self.len())
    }

    fn valid(&self, _: &SummaryConfig) -> bool {
        true
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, DI)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let mut all_exts = Exts::empty();

        let mut out_data: Vec<DI> = Vec::with_capacity(items.size_hint().0);

        let mut nobs = 0i32;
        for (_, exts, d) in items {
            out_data.push(d);
            all_exts = all_exts.add(exts);
            nobs += 1;
        }

        out_data.sort();
        out_data.dedup();
        out_data.shrink_to_fit();
        
        (nobs as usize >= config.min_kmer_obs, all_exts, out_data)
    }

}

/// For CountFilterComb
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
// aligned would be 16 Bytes, packed is 12 Bytes
//#[repr(packed)]
pub struct TagsSumData {
    tags: Tags,
    sum: i32,
}

impl SummaryData<u8, (Tags, i32)> for TagsSumData {
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

    fn count(&self) -> Option<usize> {
        Some(self.sum as usize)
    }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> {
        None
    }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> {
        None
    }

    fn sample_count(&self) -> Option<usize> {
        Some(self.tags.len())
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        valid(self.tags, self.sum, config)
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, u8)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
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

        let tags = Tags::from_u8_vec(out_data);
        
        (valid(tags, nobs, config), all_exts, TagsSumData::new((tags, nobs)))
    }

}

/// For CountFilterStats
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TagsCountsSumData {
    tags: Tags,
    sum: i32,
    counts: Box<[u32]>,
}

impl SummaryData<u8, (Tags, Box<[u32]>, i32)> for TagsCountsSumData {
    fn new(data: (Tags, Box<[u32]>, i32)) -> Self {
        TagsCountsSumData { tags: data.0, counts: data.1, sum: data.2 }
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

    fn count(&self) -> Option<usize> {
        Some(self.sum as usize)
    }

    fn sample_count(&self) -> Option<usize> {
        Some(self.counts.len() as usize)
    }

    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {
        let p = match config.stat_test {
            StatTest::TTest => t_test(&self.tags.to_u8_vec(), &self.counts.to_vec(), &config.sample_info),
            StatTest::UTest => u_test(&self.tags.to_u8_vec(), &self.counts.to_vec(), &config.sample_info),
        };
        
        match p {
            Ok(p_value) => Some(p_value),
            Err(_) => None,
        }
    }

    fn fold_change(&self, config: &SummaryConfig) -> Option<f32> {
        Some(log2_fold_change(self.tags, self.counts.to_vec(), &config.sample_info))
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        let valid_p = match config.max_p {
            Some(p) => self.p_value(config).expect("error calculating p_value") <= p,
            None => true,
        };

        valid(self.tags, self.sum, config) && valid_p
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, u8)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
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

        let valid_p = match config.max_p {
            Some(max_p) => {
                let p_value = match config.stat_test {
                    StatTest::TTest => t_test(&out_data, &tag_counts, &config.sample_info),
                    StatTest::UTest => u_test(&out_data, &tag_counts, &config.sample_info),
                };
                match p_value {
                    Ok(p) => p <= max_p,
                    Err(_) => true
                } 
            },
            None => true,
        };        

        let tag_counts: Box<[u32]> = tag_counts.into();
        let tags = Tags::from_u8_vec(out_data);

        (valid(tags, nobs, config) && valid_p, all_exts, TagsCountsSumData::new((tags, tag_counts, nobs))) 
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

impl SummaryData<u8, (Tags, Box<[u32]>)> for TagsCountsData{
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

    fn count(&self) -> Option<usize> {
        Some(self.counts.iter().sum::<u32>() as usize)
    }

    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {
        let p = match config.stat_test {
            StatTest::TTest => t_test(&self.tags.to_u8_vec(), &self.counts.to_vec(), &config.sample_info),
            StatTest::UTest => u_test(&self.tags.to_u8_vec(), &self.counts.to_vec(), &config.sample_info),
        };
        
        match p {
            Ok(p_value) => Some(p_value),
            Err(_) => None,
        }
    }

    fn fold_change(&self, config: &SummaryConfig) -> Option<f32> {
        Some(log2_fold_change(self.tags, self.counts.to_vec(), &config.sample_info))
    }

    fn sample_count(&self) -> Option<usize> {
        Some(self.counts.len())
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        let valid_p = match config.max_p {
            Some(p) => self.p_value(config).expect("error calculating p-value") <= p,
            None => true,
        }; 

        valid(self.tags, self.sum(), config) && valid_p
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, u8)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
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

        let valid_p = match config.max_p {
            Some(max_p) => {
                let p_value = match config.stat_test {
                    StatTest::TTest => t_test(&out_data, &tag_counts, &config.sample_info),
                    StatTest::UTest => u_test(&out_data, &tag_counts, &config.sample_info),
                };
                match p_value {
                    Ok(p) => p <= max_p,
                    Err(_) => true
                } 
            },
            None => true,
        };        

        let tag_counts: Box<[u32]> = tag_counts.into();
        let tags = Tags::from_u8_vec(out_data);

        (valid(tags, nobs, config) && valid_p, all_exts, TagsCountsData::new((tags, tag_counts))) 
    }

}

/// Structure for Summarizer Data which contains the tags, the count for each tag, and a p-value regarding groups
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TagsCountsPData {
    tags: Tags,
    counts: Box<[u32]>,
    p_value: f32
}

impl TagsCountsPData {
    #[inline]
    pub fn sum(&self) -> i32 {
        self.counts.iter().sum::<u32>() as i32
    }
}

impl SummaryData<u8, (Tags, Box<[u32]>, f32)> for TagsCountsPData{
    fn new(data: (Tags, Box<[u32]>, f32)) -> Self {
        TagsCountsPData { tags: data.0, counts: data.1, p_value: data.2 }
    }

    fn print(&self, tag_translator: &BiMap<String, u8>) -> String {
        format!("tags: {:?}, counts: {:?}, p-value: {}", self.tags.to_string_vec(tag_translator), self.counts, self.p_value).replace("\"", "\'")
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

    fn count(&self) -> Option<usize> {
        Some(self.sum() as usize)
    }

    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {
        if config.stat_test_changed {
            Some(match config.stat_test {
                StatTest::TTest => t_test(&self.tags.to_u8_vec(), &self.counts.to_vec(), &config.sample_info).unwrap(),
                StatTest::UTest => u_test(&self.tags.to_u8_vec(), &self.counts.to_vec(), &config.sample_info).unwrap(),
            })
        } else {
            Some(self.p_value)
        }
        
    }

    fn fold_change(&self, config: &SummaryConfig) -> Option<f32> {
        Some(log2_fold_change(self.tags, self.counts.to_vec(), &config.sample_info))
    }
    
    fn sample_count(&self) -> Option<usize> {
        Some(self.counts.len())
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        let valid_p = match config.max_p {
            Some(p) => self.p_value <= p,
            None => true,
        }; 

        valid(self.tags, self.sum(), config) && valid_p
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, u8)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
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

        // caluclate p-value with t-test
        let p_value = match config.stat_test {
            StatTest::TTest => t_test(&out_data, &tag_counts, &config.sample_info),
            StatTest::UTest => u_test(&out_data, &tag_counts, &config.sample_info),
        }.unwrap();
         

        let tag_counts: Box<[u32]> = tag_counts.into();
        let tags = Tags::from_u8_vec(out_data);

        let valid = match config.max_p {
            Some(p) => valid(tags, nobs, config) && (p_value <= p),
            None => valid(tags, nobs, config),
        };

        (valid, all_exts, TagsCountsPData::new((tags, tag_counts, p_value))) 
    }

}

/// Structure for Summarizer Data which contains the tags and the count for each tag
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TagsCountsEMData {
    tags: Tags,
    counts: Box<[u32]>,
    edge_mults: EdgeMult,
}

impl TagsCountsEMData {
    #[inline]
    pub fn sum(&self) -> i32 {
        self.counts.iter().sum::<u32>() as i32
    }
}

impl SummaryData<u8, (Tags, Box<[u32]>, EdgeMult)> for TagsCountsEMData{
    fn new(data: (Tags, Box<[u32]>, EdgeMult)) -> Self {
        TagsCountsEMData { tags: data.0, counts: data.1, edge_mults: data.2 }
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

    fn count(&self) -> Option<usize> {
        Some(self.counts.iter().sum::<u32>() as usize)
    }

    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {
        let p = match config.stat_test {
            StatTest::TTest => t_test(&self.tags.to_u8_vec(), &self.counts.to_vec(), &config.sample_info),
            StatTest::UTest => u_test(&self.tags.to_u8_vec(), &self.counts.to_vec(), &config.sample_info),
        };
        
        match p {
            Ok(p_value) => Some(p_value),
            Err(_) => None,
        }
    }

    fn fold_change(&self, config: &SummaryConfig) -> Option<f32> {
        Some(log2_fold_change(self.tags, self.counts.to_vec(), &config.sample_info))
    }

    fn sample_count(&self) -> Option<usize> {
        Some(self.counts.len())
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        let valid_p = match config.max_p {
            Some(p) => self.p_value(config).expect("error calculating p-value") <= p,
            None => true,
        }; 

        valid(self.tags, self.sum(), config) && valid_p
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, u8)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let mut all_exts = Exts::empty();
        let mut edge_mults = EdgeMult::new();

        let mut out_data: Vec<u8> = Vec::with_capacity(items.size_hint().0);

        let mut nobs = 0i32;
        for (_, exts, d) in items {
            out_data.push(d); 
            all_exts = all_exts.add(exts);
            edge_mults.add_exts(exts);
            nobs += 1;
        }

        assert_eq!(nobs as u32, edge_mults.sum());

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

        let valid_p = match config.max_p {
            Some(max_p) => {
                let p_value = match config.stat_test {
                    StatTest::TTest => t_test(&out_data, &tag_counts, &config.sample_info),
                    StatTest::UTest => u_test(&out_data, &tag_counts, &config.sample_info),
                };
                match p_value {
                    Ok(p) => p <= max_p,
                    Err(_) => true
                } 
            },
            None => true,
        };        

        let tag_counts: Box<[u32]> = tag_counts.into();
        let tags = Tags::from_u8_vec(out_data);

        (valid(tags, nobs, config) && valid_p, all_exts, TagsCountsEMData::new((tags, tag_counts, edge_mults))) 
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

impl SummaryData<u8, (u32, u32)> for GroupCountData {
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

    fn count(&self) -> Option<usize> {
        Some((self.group1 + self.group2) as usize)
    }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> {
        None
    }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> {
        None
    }

    fn sample_count(&self) -> Option<usize> {
        None
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        self.sum() >= config.min_kmer_obs as u32
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, u8)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let mut all_exts = Exts::empty();
        let mut count1 = 0;
        let mut count2 = 0;

        let mut nobs = 0i32;
        for (_, exts, d) in items {
            let tag = (2 as M).pow(d as u32);
            let group1 = ((config.sample_info.marker0 & tag) > 1) as u32;
            let group2 = ((config.sample_info.marker1 & tag) > 1) as u32;

            if (group1 + group2) != 1 { 
                panic!(
                    "should not happen\n tag: {:#066b}\n m1: {:#066b}\n m2: {:#066b}\n g1: {:#066b}\n g2: {:#066b}", 
                    tag, config.sample_info.marker0, config.sample_info.marker1, group1, group2
                )
            }
            count1 += group1;
            count2 += group2;
            nobs += 1;
            all_exts = all_exts.add(exts);
        }

        let counts = match config.significant {
            Some(digits) => (round_digits(count1, digits), round_digits(count2, digits)),
            None => (count1, count2)
        };

        assert_eq!((count1 + count2),nobs as u32);
        (nobs as usize >= config.min_kmer_obs, all_exts, GroupCountData::new(counts))
    }

}


/// Data container to store the relative amount of k-mers (in percent) in group 1 and the overall count
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RelCountData {
    percent: u32,
    count: u32
}

impl SummaryData<u8, (u32, u32)> for RelCountData {
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

    fn count(&self) -> Option<usize> {
        Some(self.count as usize)
    }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> {
        None
    }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> {
        None
    }

    fn sample_count(&self) -> Option<usize> {
        None
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        self.count >= config.min_kmer_obs as u32
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, u8)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let mut all_exts = Exts::empty();
        let mut count1 = 0;
        let mut count2 = 0;

        let mut nobs = 0i32;
        for (_, exts, d) in items {
            let tag = (2 as M).pow(d as u32);
            let group1 = ((config.sample_info.marker0 & tag) > 1) as u32;
            let group2 = ((config.sample_info.marker1 & tag) > 1) as u32;

            if (group1 + group2) != 1 { 
                panic!(
                    "should not happen\n tag: {:#066b}\n m1: {:#066b}\n m2: {:#066b}\n g1: {:#066b}\n g2: {:#066b}", 
                    tag, config.sample_info.marker0, config.sample_info.marker1, group1, group2
                )
            }
            count1 += group1;
            count2 += group2;
            all_exts = all_exts.add(exts);
            nobs += 1;
        }

        assert_eq!(count1 + count2, nobs as u32);


        let percent = (count1 as f64 / nobs as f64 * 100.) as u32;
        let count = match config.significant {
            Some(digits) => round_digits(nobs as u32, digits),
            None => nobs as u32
        };

        (nobs as usize >= config.min_kmer_obs, all_exts, RelCountData::new((percent, count))) 
    }
    
}

/// Implement this trait to control how multiple observations of a kmer
/// are carried forward into a DeBruijn graph.
pub trait KmerSummarizer<DI, DO: SummaryData<DI, SD>, SD> {
    /// The input `items` is an iterator over kmer observations. Input observation
    /// is a tuple of (kmer, extensions, data). The summarize function inspects the
    /// data and returns a tuple indicating:
    /// * whether this kmer passes the filtering criteria (e.g. is there a sufficient number of observation)
    /// * the accumulated Exts of the kmer
    /// * a summary data object of type `DO` that will be used as a color annotation in the DeBruijn graph.
    /// 
    
    fn new(min_kmer_obs: usize, sample_info: SampleInfo) -> Self;
    #[deprecated]
    fn summarize<K: Kmer, F: Iterator<Item = (K, Exts, DI)>>(&self, items: F, significant: Option<u32>) -> (bool, Exts, DO);
}

/// A simple KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times. The metadata returned about a Kmer
/// is the number of times it was observed, capped at 2^32.
pub struct CountFilter<D> {
    min_kmer_obs: usize,
    phantom: PhantomData<D>
}


impl<D> KmerSummarizer<D, u32, u32> for CountFilter<D> {
    fn new(min_kmer_obs: usize, _: SampleInfo) -> Self {
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
    fn new(min_kmer_obs: usize, _: SampleInfo) -> Self {
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
/// are the tags the k-mer occured with as a `Tags` and its numver if observations.
pub struct CountFilterComb {
    min_kmer_obs: usize,
    phantom: PhantomData<u8>,
}

impl KmerSummarizer<u8, TagsSumData, (Tags, i32)> for CountFilterComb {
    fn new(min_kmer_obs: usize, _: SampleInfo) -> Self {
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
/// are the tags the k-mer was observed with as a `Tags`, how many times 
/// it was observed with each tag, and how many times it was observed overall
pub struct CountFilterStats {
    min_kmer_obs: usize,
    phantom: PhantomData<u8>,
}

impl KmerSummarizer<u8, TagsCountsSumData, (Tags, Box<[u32]>, i32)> for CountFilterStats {
    fn new(min_kmer_obs: usize, _: SampleInfo) -> Self {
        CountFilterStats {
            min_kmer_obs,
            phantom: PhantomData,
        }
    }

    fn summarize<K: Kmer, F: Iterator<Item = (K, Exts, u8)>>(&self, items: F, _: Option<u32>) -> (bool, Exts, TagsCountsSumData) {
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


        (nobs as usize >= self.min_kmer_obs, all_exts, TagsCountsSumData::new((Tags::from_u8_vec(out_data), tag_counts, nobs))) 
    }
}


/// A simple KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times. The metadata returned about a Kmer
/// are the tags the k-mer was observed with as a `Tags`, and how many times 
/// it was observed with each tag.
pub struct CountsFilterStats {
    min_kmer_obs: usize,
    phantom: PhantomData<u8>,
}

impl KmerSummarizer<u8, TagsCountsData, (Tags, Box<[u32]>)> for CountsFilterStats {
    fn new(min_kmer_obs: usize, _: SampleInfo) -> Self {
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
/// is how many times the k-mer was observed with tags from each the two groups
/// saved in the `SampleInfo`.
pub struct CountsFilterGroups {
    min_kmer_obs: usize,
    marker: SampleInfo,
    phantom: PhantomData<u8>,
}

impl KmerSummarizer<u8, GroupCountData, (u32, u32)> for CountsFilterGroups {
    fn new(min_kmer_obs: usize, marker: SampleInfo) -> Self {
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
/// is how many percent of the observations of the k-mer were with tags
/// from group 1 in the `SanpleInfo` and how many times the k-mer was 
/// observed overall.
pub struct CountsFilterRel{
    min_kmer_obs: usize,
    marker: SampleInfo,
    phantom: PhantomData<u8>,
}

impl KmerSummarizer<u8, RelCountData, (u32, u32)> for CountsFilterRel {
    fn new(min_kmer_obs: usize, marker: SampleInfo) -> Self {
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


/// A KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times and that were observed with at least 
/// a third of the tags in at least one of the groups in the `SampleInfo`. 
/// The metadata returned about a Kmer are the tags the k-mer was observed 
/// with as a `Tags`, and how many times it was observed with each tag.
pub struct CountsFilterMaj {
    min_kmer_obs: usize,
    marker: SampleInfo,
    phantom: PhantomData<u8>,
}

impl KmerSummarizer<u8, TagsCountsData, (Tags, Box<[u32]>)> for CountsFilterMaj {
    fn new(min_kmer_obs: usize, marker: SampleInfo) -> Self {
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

        let tags = Tags::from_u8_vec(out_data);

        // get amount of labels in tags from each group
        let dist0= tags.bit_and_dist(self.marker.marker0);
        let dist1= tags.bit_and_dist(self.marker.marker1);

        assert_eq!(dist0 + dist1, tags.to_u8_vec().len());

        // valid if:
        // - n obs >= min obs AND
        // - observed in at least one third of samples in one group
        let valid = (nobs as usize >= self.min_kmer_obs)
            && ((dist0 as f32 / self.marker.count0 as f32 > 0.333333) 
                | (dist1 as f32 / self.marker.count1 as f32 > 0.333333));

        let tag_counts: Box<[u32]> = tag_counts.into();

        (valid, all_exts, TagsCountsData::new((tags, tag_counts))) 
    }
}


/// A KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times and that were observed with at least 
/// a third of the tags in both of the groups in the `SampleInfo`. 
/// The metadata returned about a Kmer are the tags the k-mer was observed 
/// with as a `Tags`, and how many times it was observed with each tag.
pub struct CountsFilterMajB {
    min_kmer_obs: usize,
    marker: SampleInfo,
    phantom: PhantomData<u8>,
}

impl KmerSummarizer<u8, TagsCountsData, (Tags, Box<[u32]>)> for CountsFilterMajB {
    fn new(min_kmer_obs: usize, marker: SampleInfo) -> Self {
        CountsFilterMajB {
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

        let tags = Tags::from_u8_vec(out_data);

        // get amount of labels in tags from each group
        let dist0= tags.bit_and_dist(self.marker.marker0);
        let dist1= tags.bit_and_dist(self.marker.marker1);

        assert_eq!(dist0 + dist1, tags.to_u8_vec().len());

        // valid if:
        // - n obs >= min obs AND
        // - observed in at least one third of samples in one group
        let valid = (nobs as usize >= self.min_kmer_obs)
            && (dist0 as f32 / self.marker.count0 as f32 > 0.333333) 
            && (dist1 as f32 / self.marker.count1 as f32 > 0.333333);

        let tag_counts: Box<[u32]> = tag_counts.into();

        (valid, all_exts, TagsCountsData::new((tags, tag_counts))) 
    }
}


/// A KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times. 
/// The metadata returned about a Kmer are the tags the k-mer was observed 
/// with as a `Tags`, how many times it was observed with each tag, and a 
/// p-value from an unpaired t-test comparing the groups in the `SampleInfo`.
pub struct CountsFilterStat {
    min_kmer_obs: usize,
    sample_info: SampleInfo,
    phantom: PhantomData<u8>,
}

impl KmerSummarizer<u8, TagsCountsPData, (Tags, Box<[u32]>, f32)> for CountsFilterStat {
    fn new(min_kmer_obs: usize, sample_info: SampleInfo) -> Self {
        CountsFilterStat {
            min_kmer_obs,
            sample_info,
            phantom: PhantomData,
        }
    }

    fn summarize<K: Kmer, F: Iterator<Item = (K, Exts, u8)>>(&self, items: F, _: Option<u32>) -> (bool, Exts, TagsCountsPData) {
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

        // caluclate p-value

        let mut counts_g0 = Vec::new();
        let mut counts_g1 = Vec::new();


        let (m0, m1) = self.sample_info.get_marker();

        for (label, count) in out_data.iter().zip(&tag_counts) {
            let bin_rep = (2 as M).pow(*label as u32);
            let norm = *count as f64 / self.sample_info.sample_kmers[*label as usize] as f64;
            if (m0 & bin_rep) > 0 { counts_g0.push(norm); }
            if (m1 & bin_rep) > 0 { counts_g1.push(norm); }
        }

        let n0 = self.sample_info.count0 as f64;
        let n1 = self.sample_info.count1 as f64;

        let mean0 = counts_g0.iter().sum::<f64>() as f64 / n0;
        let mean1 = counts_g1.iter().sum::<f64>() as f64 / n1;

        let df = n0 + n1 - 2.;

        let var0 = (counts_g0.iter().map(|count| (*count - mean0).powf(2.)).sum::<f64>() + (n0 - counts_g0.len() as f64) * mean0.powf(2.)) / (n0 - 1.);
        let var1 = (counts_g1.iter().map(|count| (*count - mean1).powf(2.)).sum::<f64>() + (n1 - counts_g1.len() as f64) * mean1.powf(2.)) / (n1 - 1.);

        let s = ((1./n0 + 1./n1) * ((n0 - 1.) * var0 + (n1 - 1.) * var1) / df).sqrt();

        let t = (mean0 - mean1) / s;

        let t_dist = StudentsT::new(0.0, 1.0, df).expect("error creating student dist");

        let p_value = 2. * (1. - t_dist.cdf(t.abs())) as f32;

/*         println!("n:    {}      {}", n0, n1);
        println!("mean: {}      {}", mean0, mean1);
        println!("var:  {}      {}", var0, var1);
        println!("s:    {}", s);
        println!("df:   {}", df);
        println!("t:    {}", t); */
        

        let tag_counts: Box<[u32]> = tag_counts.into();

        (nobs as usize >= self.min_kmer_obs, all_exts, TagsCountsPData::new((Tags::from_u8_vec(out_data), tag_counts, p_value))) 
    }
}

#[cfg(test)]
mod test {
    use std::{fmt::Debug, fs::File, io::BufReader};

    use bimap::BiMap;
    use boomphf::hashmap::BoomHashMap2;
    use rand::Rng;

    use crate::{clean_graph::CleanGraph, compression::{ compress_graph, compress_kmers_with_hash, CompressionSpec, ScmapCompress}, dna_string::DnaString, filter::filter_kmers, graph::{BaseGraph, DebruijnGraph, Node}, kmer::{Kmer12, Kmer16, Kmer8}, reads::Reads, summarizer::{self, t_test, u_test, GroupCountData, NotEnoughSamplesError, RelCountData, SampleInfo, SummaryData, TagsCountsData, TagsCountsPData, GroupFrac}, Exts, Kmer, Tags};

    use super::{log2_fold_change, SummaryConfig, TagsCountsSumData};

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

        // markers: 
        let sample_info = SampleInfo { marker0: 2, marker1: 12, count0: 1, count1: 2, sample_kmers: Vec::new() };

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

        let config = SummaryConfig::new(min_kmer_obs, significant, GroupFrac::None, 0.33, sample_info.clone(), None, summarizer::StatTest::TTest);


        let pr = false;

        if pr { println!("{}", reads) };


        //construct and compress graph with CountFilter
        let spec: ScmapCompress<u32> = ScmapCompress::new();
        let graph: DebruijnGraph<K, u32> = test_summarizer(&reads, &config, spec);

        graph.to_gfa_with_tags("test_cf.gfa",|n| <u32 as summarizer::SummaryData<u8, u32>>::print(n.data(), &tag_translator)).unwrap();

        println!("cf graph size: {}", graph.len());


        // construct and compress graph with CountsFilterGroups
        let spec: ScmapCompress<GroupCountData> = ScmapCompress::new();
        let graph: DebruijnGraph<K, GroupCountData> = test_summarizer(&reads, &config, spec);

        //graph.to_gfa_with_tags("test_csfg.gfa",|n| n.data().print(&tag_translator)).unwrap();

        println!("csfg graph size: {}", graph.len());


        // construct and compress graph with CountsFilterRel
        let spec: ScmapCompress<RelCountData> = ScmapCompress::new();
        let graph: DebruijnGraph<K, RelCountData> = test_summarizer(&reads, &config, spec);

        //graph.to_gfa_with_tags("test_csfr.gfa",|n| n.data().print(&tag_translator)).unwrap();

        println!("csfr graph size: {}", graph.len());


        // construct and compress graph with CountsFilterStats
        let spec: ScmapCompress<TagsCountsData> = ScmapCompress::new();
        let graph: DebruijnGraph<K, TagsCountsData> = test_summarizer(&reads, &config, spec);

        graph.to_gfa_with_tags("test_csfs.gfa",|n| n.data().print(&tag_translator)).unwrap();

        println!("csfs graph size: {}", graph.len());



        // same but with less significant digits
        let significant= Some(1);

        let config = SummaryConfig::new(min_kmer_obs, significant, GroupFrac::None, 0.33, sample_info.clone(), None, summarizer::StatTest::TTest);


        //construct and compress graph with CountFilter
        let spec: ScmapCompress<u32> = ScmapCompress::new();
        let graph: DebruijnGraph<K, u32> = test_summarizer(&reads, &config, spec);

        graph.to_gfa_with_tags("test_cf-1.gfa",|n| <u32 as summarizer::SummaryData<u8, u32>>::print(n.data(), &tag_translator)).unwrap();

        println!("cf graph size: {}", graph.len());


        // construct and compress graph with CountsFilterGroups
        let spec: ScmapCompress<GroupCountData> = ScmapCompress::new();
        let graph: DebruijnGraph<K, GroupCountData> = test_summarizer(&reads, &config, spec);

        //graph.to_gfa_with_tags("test_csfg.gfa",|n| n.data().print(&tag_translator)).unwrap();

        println!("csfg graph size: {}", graph.len());


        // construct and compress graph with CountsFilterRel
        let spec: ScmapCompress<RelCountData> = ScmapCompress::new();
        let graph: DebruijnGraph<K, RelCountData> = test_summarizer(&reads, &config, spec);

        //graph.to_gfa_with_tags("test_csfr.gfa",|n| n.data().print(&tag_translator)).unwrap();

        println!("csfr graph size: {}", graph.len());


        // construct and compress graph with CountsFilterStats
        let spec: ScmapCompress<TagsCountsData> = ScmapCompress::new();
        let graph: DebruijnGraph<K, TagsCountsData> = test_summarizer(&reads, &config, spec);

        //graph.to_gfa_with_tags("test_csfs.gfa",|n| n.data().print(&tag_translator)).unwrap();

        println!("csfs graph size: {}", graph.len());


    }


    fn test_summarizer<K, CS, D, SD>(reads: &Reads<u8>, config: &SummaryConfig, spec: CS) -> DebruijnGraph<K, D> 
    where  
        K: Kmer + Send + Sync, 
        CS: CompressionSpec<D> + Send + Sync, 
        D: SummaryData<u8, SD> + Debug + PartialEq + Clone + Send + Sync,
    {
        // construct and compress graph with CountsFilterStats
        //let summarizer= CountsFilterStats::new(min_kmer_obs, markers);

        let memory = 1;
        let pr = false;

        let (k_mers, _): (BoomHashMap2<K, Exts, D>, _)  = filter_kmers(&reads,
            config,
            false, 
            false, 
            memory, 
            false, 
            false, 
        );
        if pr { println!("kmers CountsFilterStats: {:?}", k_mers) };

        let graph = compress_kmers_with_hash(false, &spec, &k_mers, false, false, false);
        if pr { println!("graph CountsFilterStats: {:?}\n", graph) };

        //graph.finish().to_gfa_with_tags("test_csfs.gfa",|n| n.data().print(&tt)).unwrap();
        graph.finish()
    }

    #[test]
    fn test_maj_summarizers() {
        // generate three reads
        /* let mut reads = Reads::new();

        let dnas = ["CGATCGAGCTACTGCGACGGACGATTTTTCGAGCGGCGATTTCTCGAGGCGAGCGTCAGC".as_bytes(),
            "CGATCGAGCTACTGCGACGGACGATGACTAGCTAGCTTTTCTCGAGGCGAGCGTCAGC".as_bytes(),
            "ACTGCGACGGACGATGACTAGCTAGCTTTTCTCGAGGCGAGCGTCAGCACGATGCTAGCTGACTAGC".as_bytes(),
            "CGATCGAGCTACTGCGACGGACGATTTTTCGAGCGGCGATTTCTCGAGGCGAGCGTCAGCGGACTAGCGAG".as_bytes(),
            "ACGGACGATTTTTCGAGCGGCGATTTCTCGAGGCGAGCGTCAGC".as_bytes(),
        ];

        let repeats = 4;

        let mut rng = rand::thread_rng();
        for _i in 0..repeats {
            for dna in dnas {
            reads.add_from_bytes(dna, Exts::empty(), rng.gen_range(0, 12));
            }
        }

        let dnas = ["ACGATGCTAGCTAGCTGACTGACGATCGTAGCTAGCTGATCGGATCGATC".as_bytes(),
            "ACGATGCTAGCTAGCTGACTGACGATCGTAGCTAGCTGATCGGATCGATCGACTGATCGATGCATCGACGATC".as_bytes(),
            "GACGACTGATCGAGCTACGAGCACGATGCTAGCTAGCTGACTGACGATCGTAGCTAGCTGATCGGATGCATCGACGATC".as_bytes(),
        ];

        let mut rng = rand::thread_rng();
        for _i in 0..repeats {
            for dna in dnas {
            reads.add_from_bytes(dna, Exts::empty(), rng.gen_range(0, 5));
            }
        }

        let dnas = ["GACGCCGAGATCATTATCGCGCGCGTATTATATATCGCGGAATATATCCG".as_bytes(),
            "GACGCCGAGATCATTATCGCGCGCGTATTATATGCGTAATATATATTCGGCATTAATCGCGGAATATATCCG".as_bytes(),
            "GACGCCGAGATCATTATCGCGCGCGTATTATATATCGCGGAATATATCCGAGCTGACTGATCGA".as_bytes(),
        ];

        let mut rng = rand::thread_rng();
        for _i in 0..repeats {
            for dna in dnas {
            reads.add_from_bytes(dna, Exts::empty(), rng.gen_range(5, 12));
            }
        } */

        type K = Kmer12;


        let mut rng = rand::thread_rng();
        let repeats = 8;
        let mut reads = Reads::new();

        let dnas = [
            "AGCTAGCTAGCTACGA".as_bytes(), 
            "GCATCGATGCACTGACGACT".as_bytes(),
            "CGATCGATGCTACGACTAGCGACTG".as_bytes()
        ];

        let mut sample_kmers = vec![0; 12];

        let mut dna_kmers = Vec::new();

        for dna in dnas {
            let kmers = dna.len().saturating_sub(K::k() -1) as u64;
            dna_kmers.push(kmers);
        }


        for _i in 0..repeats {
            let group0 = rng.gen_range(0, 12);
            reads.add_from_bytes(dnas[0], Exts::empty(), group0 as u8);
            sample_kmers[group0] += dna_kmers[0];

            let group1 = rng.gen_range(0, 5);
            reads.add_from_bytes(dnas[1], Exts::empty(), group1 as u8);
            sample_kmers[group1] += dna_kmers[1];
            
            let group2 = rng.gen_range(5, 12);
            reads.add_from_bytes(dnas[2], Exts::empty(), group2 as u8);
            sample_kmers[group2] += dna_kmers[2];
        
        }

        
        println!("sample kmers: {:?}", sample_kmers);
        

        /*
        markers: 
        - 111111100000 = 4064
        - 000000011111 = 31*/
        



        let sample_info = SampleInfo { marker0: 0b111111100000, marker1: 31, count0: 7, count1: 5, sample_kmers};
        let config = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), None, summarizer::StatTest::TTest);

        
        println!("markers: {:?}", sample_info);

        // construct and compress graph with CountsFilterMaj
        let spec: ScmapCompress<u32> = ScmapCompress::new();
        let graph: DebruijnGraph<K, u32> = test_summarizer(&reads, &config, spec);

        graph.print();

        // construct and compress graph with CountsFilterMaj
        let spec: ScmapCompress<TagsCountsData> = ScmapCompress::new();
        let graph: DebruijnGraph<K, TagsCountsData> = test_summarizer(&reads, &config, spec);

        graph.print();

        // construct and compress graph with CountsFilterMaj
        let spec: ScmapCompress<TagsCountsData> = ScmapCompress::new();
        let graph: DebruijnGraph<K, TagsCountsData> = test_summarizer(&reads, &config, spec);

        graph.print();

        // construct and compress graph with CountsFilterMajB
        let spec: ScmapCompress<TagsCountsData> = ScmapCompress::new();
        let graph: DebruijnGraph<K, TagsCountsData> = test_summarizer(&reads, &config, spec);

        graph.print();

        // construct and compress graph with CountsFilterStat
        let spec: ScmapCompress<TagsCountsPData> = ScmapCompress::new();
        let graph: DebruijnGraph<K, TagsCountsPData> = test_summarizer(&reads, &config, spec);

        graph.print();

    }

    #[test]
    fn test_p_value() {

        /*
        group 1: 111111100000 = 4064
        group 2: 000000011111 = 31
         */

        let sample_kmers = vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
        let sample_info = SampleInfo::new(31, 4064, 5, 7, sample_kmers);
        let config_t = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), Some(0.1), summarizer::StatTest::TTest);
        let config_u = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), Some(0.1), summarizer::StatTest::UTest);


        let input = [
            (Kmer8::from_u64(12), Exts::new(1), 1u8),
            (Kmer8::from_u64(12), Exts::new(1), 2u8),
            (Kmer8::from_u64(12), Exts::new(1), 3u8),
            (Kmer8::from_u64(12), Exts::new(1), 4u8),
            (Kmer8::from_u64(12), Exts::new(1), 0u8),
            (Kmer8::from_u64(12), Exts::new(1), 8u8),           
        ];

        let summarized_t = TagsCountsPData::summarize(input.into_iter(), &config_t);
        println!("{:?}", summarized_t);

        assert_eq!((summarized_t.2.p_value * 10000.).round() as u32, 5);

        let summarized_u = TagsCountsPData::summarize(input.into_iter(), &config_u);
        println!("{:?}", summarized_u);

        assert_eq!((summarized_u.2.p_value * 1000.).round() as u32, 5);

        let summarized_t = TagsCountsData::summarize(input.into_iter(), &config_t);
        println!("{:?}", summarized_t);

        assert_eq!((summarized_t.2.p_value(&config_t).unwrap() * 10000.).round() as u32, 5);

        let summarized_u = TagsCountsData::summarize(input.into_iter(), &config_u);
        println!("{:?}", summarized_u);

        assert_eq!((summarized_u.2.p_value(&config_u).unwrap() * 1000.).round() as u32, 5);

        let summarized_t = TagsCountsSumData::summarize(input.into_iter(), &config_t);
        println!("{:?}", summarized_t);

        assert_eq!((summarized_t.2.p_value(&config_t).unwrap() * 10000.).round() as u32, 5);

        let summarized_u = TagsCountsSumData::summarize(input.into_iter(), &config_u);
        println!("{:?}", summarized_u);

        assert_eq!((summarized_u.2.p_value(&config_u).unwrap() * 1000.).round() as u32, 5);


        let input = [
            (Kmer8::from_u64(12), Exts::new(1), 7u8),
            (Kmer8::from_u64(12), Exts::new(1), 8u8),
            (Kmer8::from_u64(12), Exts::new(1), 9u8),
            (Kmer8::from_u64(12), Exts::new(1), 10u8),
            (Kmer8::from_u64(12), Exts::new(1), 0u8),
        ];

        let summarized_t = TagsCountsPData::summarize(input.into_iter(), &config_t);
        println!("{:?}", summarized_t);

        assert_eq!((summarized_t.2.p_value * 10000.).round() as u32, 2345);


        let input = [
            (Kmer8::from_u64(12), Exts::new(1), 7u8),
            (Kmer8::from_u64(12), Exts::new(1), 8u8),
            (Kmer8::from_u64(12), Exts::new(1), 9u8),
            (Kmer8::from_u64(12), Exts::new(1), 10u8),
            (Kmer8::from_u64(12), Exts::new(1), 1u8),
            (Kmer8::from_u64(12), Exts::new(1), 0u8),
        ];

        let summarized_t = TagsCountsPData::summarize(input.into_iter(), &config_t);
        println!("{:?}", summarized_t);

        assert_eq!((summarized_t.2.p_value * 10000.).round() as u32, 5995);

        let summarized_u = TagsCountsPData::summarize(input.into_iter(), &config_u);
        println!("{:?}", summarized_u);

        assert_eq!((summarized_u.2.p_value * 1000.).round() as u32, 575);


        let input = [
            (Kmer8::from_u64(12), Exts::new(1), 7u8),
            (Kmer8::from_u64(12), Exts::new(1), 8u8),
            (Kmer8::from_u64(12), Exts::new(1), 9u8),
            (Kmer8::from_u64(12), Exts::new(1), 10u8),
            (Kmer8::from_u64(12), Exts::new(1), 11u8),
            (Kmer8::from_u64(12), Exts::new(1), 0u8),
        ];

        let summarized_t = TagsCountsPData::summarize(input.into_iter(), &config_t);
        println!("{:?}", summarized_t);

        assert_eq!((summarized_t.2.p_value * 10000.).round() as u32, 924);

        let summarized_u = TagsCountsPData::summarize(input.into_iter(), &config_u);
        println!("{:?}", summarized_u);

        assert_eq!((summarized_u.2.p_value * 1000.).round() as u32, 93);


        // with different kmer counts

        let sample_kmers = vec![12, 3345, 3478, 87, 1, 2, 666, 98111, 23982938, 555, 122, 7238];
        /*
        0: 0.08333333333333333
        1: 0.00029895366218236175
        2: 0.0002875215641173088
        3: 0.011494252873563218
        4: 1
        5: 0.5
        6: 0.0015015015015015015
        7: 0.000010192537024390741
        8: 0.00000004169630926786368
        9: 0.0018018018018018018
        10: 0.00819672131147541
        11: 0.00013815971262779773
         */

        let sample_info = SampleInfo::new(31, 4064, 5, 7, sample_kmers);
        let config_t = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), None, summarizer::StatTest::TTest);
        let config_u = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info, None, summarizer::StatTest::UTest);


        let input = [
            (Kmer8::from_u64(12), Exts::new(1), 7u8), 
            (Kmer8::from_u64(12), Exts::new(1), 8u8),
            (Kmer8::from_u64(12), Exts::new(1), 9u8),
            (Kmer8::from_u64(12), Exts::new(1), 10u8),
            (Kmer8::from_u64(12), Exts::new(1), 1u8),
            (Kmer8::from_u64(12), Exts::new(1), 0u8),
        ];

        let summarized_t = TagsCountsPData::summarize(input.into_iter(), &config_t);
        println!("{:?}", summarized_t);

        assert_eq!((summarized_t.2.p_value * 10000.).round() as u32, 2955);

        let summarized_u = TagsCountsPData::summarize(input.into_iter(), &config_u);
        println!("{:?}", summarized_u);

        assert_eq!((summarized_u.2.p_value * 1000.).round() as u32, 862);

        let summarized_t = TagsCountsData::summarize(input.into_iter(), &config_t);
        println!("{:?}", summarized_t);

        assert_eq!((summarized_t.2.p_value(&config_t).unwrap() * 10000.).round() as u32, 2955);

        let summarized_u = TagsCountsData::summarize(input.into_iter(), &config_u);
        println!("{:?}", summarized_u);

        assert_eq!((summarized_u.2.p_value(&config_u).unwrap() * 1000.).round() as u32, 862);

        /*
        group 1: 000000100 = 4
        group 2: 000000011 = 3
         */

        let sample_kmers = vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
        let sample_info = SampleInfo::new(4, 4, 1, 2, sample_kmers);
        let config_t = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), Some(0.1), summarizer::StatTest::TTest);
        let config_u = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), Some(0.1), summarizer::StatTest::UTest);


        let input = [
            (Kmer8::from_u64(12), Exts::new(1), 1u8),
            (Kmer8::from_u64(12), Exts::new(1), 0u8),
            (Kmer8::from_u64(12), Exts::new(1), 1u8),
            (Kmer8::from_u64(12), Exts::new(1), 2u8),
            (Kmer8::from_u64(12), Exts::new(1), 0u8),
            (Kmer8::from_u64(12), Exts::new(1), 1u8),           
        ];

        let summarized_t = TagsCountsData::summarize(input.into_iter(), &config_t);
        assert_eq!(summarized_t.2.p_value(&config_t), None);

        let summarized_u = TagsCountsData::summarize(input.into_iter(), &config_u);
        assert_eq!(summarized_u.2.p_value(&config_u), None);

        let summarized_t = TagsCountsSumData::summarize(input.into_iter(), &config_t);
        assert_eq!(summarized_t.2.p_value(&config_t), None);

        let summarized_u = TagsCountsSumData::summarize(input.into_iter(), &config_u);
        assert_eq!(summarized_u.2.p_value(&config_u), None);

    }

    #[test]
    fn test_mwu() {
        let tags: Vec<u8> = vec![0, 1, 2, 3, 4, 7, 8, 9];
        let test_counts: Vec<u32> = vec![3, 1, 4, 6, 9, 7, 3, 8];

        let m0: u64 = 0b1111111111110000000;
        let m1: u64 = 0b0000000000001111111;
        let c0: u32 = 12;
        let c1: u32 = 7;

        let sample_counts = vec![1; 19];

        assert_eq!(c0, m0.count_ones());
        assert_eq!(c1, m1.count_ones());

        let sample_info = SampleInfo::new(m0, m1, c0 as u8, c1 as u8, sample_counts);

        let p_value_t = t_test(&tags, &test_counts, &sample_info);
        println!("p-value t test: {}", p_value_t.unwrap());

        // soll: U = 24.5, p = 0.0995, z = -1.65, exaxt p: 0.142

        let p_value_u = u_test(&tags, &test_counts, &sample_info);
        println!("p-value u test: {}", p_value_u.unwrap());

        // test Error
        let tags: Vec<u8> = vec![0, 1];
        let test_counts: Vec<u32> = vec![3, 1];

        let m0: u64 = 0b100;
        let m1: u64 = 0b011;
        let c0: u32 = 1;
        let c1: u32 = 2;

        let sample_counts = vec![1; 3];

        assert_eq!(c0, m0.count_ones());
        assert_eq!(c1, m1.count_ones());

        let sample_info = SampleInfo::new(m0, m1, c0 as u8, c1 as u8, sample_counts);

        let p_value_t = t_test(&tags, &test_counts, &sample_info);
        assert_eq!(p_value_t, Err(NotEnoughSamplesError {}));
        println!("p-value t test: {:?}", p_value_t);

        // soll: U = 24.5, p = 0.0995, z = -1.65, exaxt p: 0.142

        let p_value_u = u_test(&tags, &test_counts, &sample_info);
        assert_eq!(p_value_u, Err(NotEnoughSamplesError {}));
        println!("p-value u test: {:?}", p_value_u);

    }

    #[test]
    fn test_data_valid() {
        let mut graph: BaseGraph<Kmer8, TagsCountsSumData> = BaseGraph::new(false);

        let tags = Tags::from_u8_vec(vec![0, 2, 6]);
        let counts: Box<[u32]> = [1, 3, 5].into();
        let sum = counts.iter().sum::<u32>() as i32;
        graph.add(DnaString::from_acgt_bytes("AAAAAAAA".as_bytes()).into_iter(), Exts::empty(), TagsCountsSumData::new((tags, counts, sum)));

        let tags = Tags::from_u8_vec(vec![0]);
        let counts: Box<[u32]> = [1].into();
        let sum = counts.iter().sum::<u32>() as i32;
        graph.add(DnaString::from_acgt_bytes("CCCCCCCC".as_bytes()).into_iter(), Exts::empty(), TagsCountsSumData::new((tags, counts, sum)));
                
        
        let graph = graph.finish();

        graph.print();

        let sample_kmers = vec![123, 234, 12334, 34];
        let sample_info = SampleInfo::new(0b00100101, 0b11011010, 3, 5, sample_kmers);
        let config = SummaryConfig::new(3, None, GroupFrac::None, 0.33,  sample_info, None, summarizer::StatTest::TTest);

        let censor_nodes = CleanGraph::new(|node: &Node<'_, Kmer8, TagsCountsSumData>| !node.data().valid(&config))
                    .find_bad_nodes(&graph);
        println!("censor nodes: {:?}", censor_nodes);
        let filter_graph = compress_graph(false, &ScmapCompress::new(), graph, Some(censor_nodes));

        filter_graph.print();

        // larger test

        let reader = BufReader::with_capacity(64*1024, File::open("400.graph.dbg").unwrap());
        let (graph, _, mut config): (DebruijnGraph<Kmer16, TagsCountsSumData>, Vec<String>, SummaryConfig) = bincode::deserialize_from(reader)
            .expect("error deserializing graph, hashed labels, and config");

        config.set_min_kmer_obs(3);

        let node38 = graph.get_node(38);
        assert!(!node38.data().valid(&config));

        let bad_nodes = graph.find_bad_nodes(|node: &Node<'_, Kmer16, TagsCountsSumData>| node.data().valid(&config));
        let bad_node_correct = vec![0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22,
            23, 24, 25, 26, 27, 28, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51,
            52, 53, 54, 55, 56, 57, 58, 60, 61, 62, 63, 64, 65, 66, 69, 71, 72, 73, 74, 76, 77, 78, 80, 82, 83, 84, 85, 
            86, 87, 89, 90, 94, 95, 98, 99, 100, 101, 102, 103, 104, 105, 107, 109, 111, 114, 116, 117, 118, 120, 122, 
            126, 127, 131, 135, 136, 138, 141, 143, 144];
        assert_eq!(bad_nodes, bad_node_correct);
        let _filtered_graph = compress_graph(false, &ScmapCompress::new(), graph, Some(bad_nodes));

    }

    #[test]
    fn test_fold_change() {    
        let marker0 = 0b0000000001111;
        let marker1 = 0b1111111110000;
        let count0 = 4;
        let count1 = 9;
        let sample_kmers = vec![2834, 2343, 12, 1234, 345345, 122, 234, 23455, 231, 2, 3564, 12344, 34555];
        let sample_info = SampleInfo::new(marker0, marker1, count0, count1, sample_kmers);

        let labels = vec![0, 1, 2, 3, 7, 8];
        let tags = Tags::from_u8_vec(labels);
        let counts = vec![1, 6, 9, 3, 6, 10];
        let fold_change = log2_fold_change(tags, counts, &sample_info);
        assert_eq!(fold_change, 5.28645313556333);

        let labels = vec![0, 6, 7, 8, 10, 11];
        let tags = Tags::from_u8_vec(labels);
        let counts = vec![12, 3, 7, 1, 22, 6];
        let fold_change = log2_fold_change(tags, counts, &sample_info);
        assert_eq!(fold_change, -1.3393244925933);
    }
}




