use bimap::BiMap;
use clap::ValueEnum;
use serde::{Deserialize, Serialize};
use statrs::distribution::{ContinuousCDF, Normal, StudentsT};
use crate::{EdgeMult, Exts, Tags, TagsCountsFormatter, TagsFormatter};
use std::{cmp::min_by, error::Error, fmt::{Debug, Display}, mem};

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

/// Configuration for summary processes
#[derive(Debug, Serialize, Deserialize, Clone)]
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
    /// 
    /// arguments: 
    /// * `min_kmer_obs`: minimum number of times a k-mer has to be observed in the reads to be valid
    /// * `significant`: some summaries round numbers to a certain number of digits ([`u32`], counts in [`GroupCountData`] and [`RelCountData`])
    /// * `group_frac`: a [`GroupFrac`] determining if the k-mers are going to be filtered out
    ///   based on if they are observed in certain percentage of the samples of each group
    /// * `frac_cutoff`: the cutoff for `group_frac`
    /// * `sample_info`: a [`SampleInfo`] with information about the sample groups
    /// * `max_p`: a maximum p-value which will be used for filtering if applicable
    /// * `stat_test`: a [`StatTest`], determining which statistical test will be used
    ///    for calculation of p-values
    pub fn new(min_kmer_obs: usize, significant: Option<u32>, group_frac: GroupFrac, frac_cutoff: f32, sample_info: SampleInfo, max_p: Option<f32>, stat_test: StatTest) -> Self {
        SummaryConfig { min_kmer_obs, significant, group_frac, frac_cutoff, sample_info, max_p, stat_test, stat_test_changed: false }
    }

    /// make an empty `SummaryConfig`
    pub fn empty() -> Self {
        SummaryConfig { min_kmer_obs: 0, significant: None, group_frac: GroupFrac::None, frac_cutoff: 0., sample_info: SampleInfo::empty(), max_p: None, stat_test: StatTest::StudentsTTest, stat_test_changed: false }
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

    /// get the binary encoded group affiliation of tags
    pub fn get_markers(&self) -> (M, M) {
        self.sample_info.get_markers()
    }

    /// get the [`SampleInfo`] stored in the `SummaryConfig`
    pub fn sample_info(&self) -> &SampleInfo {
        &self.sample_info
    }
}

/// In how many of the two sample groups does a specified percentage of the samples 
/// have to be present
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

/// Statistical test for calculation of p-values
#[derive(Copy, Clone, PartialEq, PartialOrd, ValueEnum, Debug, Serialize, Deserialize)]
pub enum StatTest {
    StudentsTTest,
    WelchsTTest,
    UTest,
}

impl std::fmt::Display for StatTest {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::StudentsTTest => write!(f, "students-t-test"),
            Self::WelchsTTest => write!(f, "welchs-t-test"),
            Self::UTest => write!(f, "u-test"),
            
        }
    }
    
}

/// contains information about the samples required for graph construction
/// 
/// ### Example:
/// 
/// - Sample IDs in group 1: 0, 1, 2
/// - Sample IDs in group 2: 3, 4, 5, 6
/// 
/// ```
/// use debruijn::summarizer::{SampleInfo, M};
/// 
/// let marker0: M = 0b0000111; // = 7
/// let marker1: M = 0b1111000; // = 120
/// 
/// let count0: u8 = 3;
/// let count1: u8 = 4;
/// 
/// let sample_kmers = vec![1232, 12323, 24342, 24234, 345456, 21234, 546456];
/// 
/// assert_eq!(marker0.count_ones(), count0 as u32);
/// assert_eq!(marker0.count_ones(), count0 as u32);
/// assert_eq!(count0 + count1, sample_kmers.len() as u8);
/// 
/// let sample_info = SampleInfo::new(marker0, marker1, count0, count1, sample_kmers);
/// 
/// ```
/// 
#[derive(Clone, Serialize, Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct SampleInfo {
    marker0: M,
    marker1: M,
    count0: u8,
    count1: u8,
    sample_kmers: Vec<u64>,
}

impl SampleInfo {
    /// make a new [`SampleInfo`]
    /// 
    /// ### Arguments
    /// * `marker0`: a [`M`] which binary-encodes the affiliation of the tags to a group
    /// * `marker1`: same as `marker0`, for a second group
    /// * `count0`: the number of samples in group 1
    /// * `count1`: the number of samples in group 2
    /// * `sample_kmers`: a [`Vec<u64>`] containing numbers of non-unique k-mers for each sample
    ///    at the index of the sample-id
    pub fn new(marker0: M, marker1: M, count0: u8, count1: u8, sample_kmers: Vec<u64>) -> Self {
        assert_eq!(marker0.count_ones(), count0 as u32);
        assert_eq!(marker0.count_ones(), count0 as u32);
        assert_eq!(count0+count1, sample_kmers.len() as u8);

        SampleInfo { marker0, marker1, count0, count1, sample_kmers }
    }

    /// make a new, empty [`SampleInfo`]
    pub fn empty() -> Self {
        SampleInfo { marker0: 0, marker1: 0, count0: 0, count1: 0, sample_kmers: Vec::new() }
    }

    /// get the binary encoded group affiliation of tags
    pub fn get_markers(&self) -> (M, M) {
        (self.marker0, self.marker1)
    }
}

/// summarize the k-mers, exts and labels
fn summarize<K, F: Iterator<Item = (K, Exts, u8)>>(items: F) -> (Exts, Vec<u8>, Vec<u32>, u32) {
    let mut all_exts = Exts::empty();

    let mut out_data: Vec<u8> = Vec::with_capacity(items.size_hint().0);

    let mut nobs = 0;
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
            tag_counts.push(tag_counter);
            tag_counter = 1;
        }
    }
    tag_counts.push(tag_counter);

    out_data.dedup();

    (all_exts, out_data, tag_counts, nobs)
}

/// summarize the k-mers, exts and labels, also include an [`EdgeMult`]
fn summarize_with_em<K, F: Iterator<Item = (K, Exts, u8)>>(items: F) -> (Exts, Vec<u8>, Vec<u32>, u32, EdgeMult) {
    let mut all_exts = Exts::empty();

    let mut out_data: Vec<u8> = Vec::with_capacity(items.size_hint().0);
    let mut edge_mults = EdgeMult::new();

    let mut nobs = 0;
    for (_, exts, d) in items {
        out_data.push(d); 
        all_exts = all_exts.add(exts);
        edge_mults.add_exts(exts);
        nobs += 1;
    }

    assert_eq!(all_exts, edge_mults.exts());

    out_data.sort();

    let mut tag_counter = 1;
    let mut tag_counts: Vec<u32> = Vec::new();

    // count the occurences of the labels
    for i in 1..out_data.len() {
        if out_data[i] == out_data[i-1] {
            tag_counter += 1;
        } else {
            tag_counts.push(tag_counter);
            tag_counter = 1;
        }
    }
    tag_counts.push(tag_counter);

    out_data.dedup();

    (all_exts, out_data, tag_counts, nobs, edge_mults)
}

/// round an unsigned integer to the specified amount of digits,
/// if the integer is shorter than the number if digits, it returns the original integer
pub fn round_digits(number: u32, digits: u32) -> u32 {
    let length = (number as f32).log10() as u32 + 1;
    if digits > length { return number }
    let empty = length - digits;
    ((number as f32/ 10i32.pow(empty) as f32).round() * 10i32.pow(empty) as f32) as u32
}

// check if the k-mer is valid according to the GroupFrac rule and its n obs
fn valid_counts(tags: Tags, nobs: u32, config: &SummaryConfig) -> bool {
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

enum PInfo<'a> {
    PValue { p: f32 },
    Calculate { out_data: &'a [u8], tag_counts: &'a Vec<u32> }
}

fn valid_p(p_info: PInfo, config: &SummaryConfig) -> bool {
    match config.max_p {
        Some(max_p) => {
            match p_info {
                PInfo::PValue { p } => p <= max_p,
                PInfo::Calculate { out_data, tag_counts } => {
                    match p_value(out_data, tag_counts, config) {
                        Ok(p) => p <= max_p,
                        Err(_) => true
                    } 
                }
            }
        },
        None => true
    }
}

fn p_value(out_data: &[u8], tag_counts: &Vec<u32>, config: &SummaryConfig) -> Result<f32, NotEnoughSamplesError> {
    match config.stat_test {
        StatTest::StudentsTTest => students_t_test(out_data, tag_counts, &config.sample_info),
        StatTest::WelchsTTest => welchs_t_test(out_data, tag_counts, &config.sample_info),
        StatTest::UTest => u_test(out_data, tag_counts, &config.sample_info),
    }
}

// perform a student's t-test
fn students_t_test(out_data: &[u8], tag_counts: &Vec<u32>, sample_info: &SampleInfo) -> Result<f32, NotEnoughSamplesError> {
    let n0 = sample_info.count0 as f64;
    let n1 = sample_info.count1 as f64;

    if (n0 < 2.) | (n1 < 2.) { return Err(NotEnoughSamplesError {})}

    let mut counts_g0 = Vec::new();
    let mut counts_g1 = Vec::new();

    let (m0, m1) = sample_info.get_markers();

    for (label, count) in out_data.iter().zip(tag_counts) {
        let bin_rep = (2 as M).pow(*label as u32);
        let norm = *count as f64 / sample_info.sample_kmers[*label as usize] as f64;
        if (m0 & bin_rep) > 0 { counts_g0.push(norm); }
        if (m1 & bin_rep) > 0 { counts_g1.push(norm); }
    }

    let mean0 = counts_g0.iter().sum::<f64>() / n0;
    let mean1 = counts_g1.iter().sum::<f64>() / n1;

    let df = n0 + n1 - 2.;

    let var0 = (counts_g0.iter().map(|count| (*count - mean0).powi(2)).sum::<f64>() + (n0 - counts_g0.len() as f64) * mean0.powi(2)) / (n0 - 1.);
    let var1 = (counts_g1.iter().map(|count| (*count - mean1).powi(2)).sum::<f64>() + (n1 - counts_g1.len() as f64) * mean1.powi(2)) / (n1 - 1.);

    let s = ((1./n0 + 1./n1) * ((n0 - 1.) * var0 + (n1 - 1.) * var1) / df).sqrt();

    let t = (mean0 - mean1) / s;

    let t_dist = StudentsT::new(0.0, 1.0, df).expect("error creating student dist: check if you have enough samples (at least 3)");

    let p_value = 2. * (1. - t_dist.cdf(t.abs())) as f32;

    Ok(p_value)
}

// perform a welch's t-test
fn welchs_t_test(out_data: &[u8], tag_counts: &Vec<u32>, sample_info: &SampleInfo) -> Result<f32, NotEnoughSamplesError> {
    let n0 = sample_info.count0 as f64;
    let n1 = sample_info.count1 as f64;

    if (n0 < 2.) | (n1 < 2.) { return Err(NotEnoughSamplesError {})}

    let mut counts_g0 = Vec::new();
    let mut counts_g1 = Vec::new();

    let (m0, m1) = sample_info.get_markers();

    for (label, count) in out_data.iter().zip(tag_counts) {
        let bin_rep = (2 as M).pow(*label as u32);
        let norm = *count as f64 / sample_info.sample_kmers[*label as usize] as f64;
        if (m0 & bin_rep) > 0 { counts_g0.push(norm); }
        if (m1 & bin_rep) > 0 { counts_g1.push(norm); }
    }

    let mean0 = counts_g0.iter().sum::<f64>() / n0;
    let mean1 = counts_g1.iter().sum::<f64>() / n1;

    let s0 = (counts_g0.iter().map(|count| (*count - mean0).powi(2)).sum::<f64>() + (n0 - counts_g0.len() as f64) * mean0.powi(2)) / (n0 - 1.);
    let s1 = (counts_g1.iter().map(|count| (*count - mean1).powi(2)).sum::<f64>() + (n1 - counts_g1.len() as f64) * mean1.powi(2)) / (n1 - 1.);

    let s0 = s0.sqrt();
    let s1 = s1.sqrt();

    let t = (mean0 - mean1) / (s0.powi(2)/n0 + s1.powi(2)/n1).sqrt();

    let df0 = n0 - 1.;
    let df1 = n1 - 1.;

    let df = ((s0.powi(2) / n0 + s1.powi(2) / n1).powi(2) / 
        (s0.powi(4) / (n0.powi(2) * df0) 
            + s1.powi(4) / (n1.powi(2) * df1))).floor();

    let t_dist = StudentsT::new(0.0, 1.0, df).expect("error creating student dist: check if you have enough samples");

    let p_value = 2. * (1. - t_dist.cdf(t.abs())) as f32;

    Ok(p_value)
}

// perform a mann-whitney-u-test
fn u_test(out_data: &[u8], tag_counts: &Vec<u32>, sample_info: &SampleInfo) -> Result<f32, NotEnoughSamplesError> {

    let n0 = sample_info.count0 as f64;
    let n1 = sample_info.count1 as f64;

    if (n0 < 2.) | (n1 < 2.) { return Err(NotEnoughSamplesError {})}


    let mut counts_g0 = Vec::new();
    let mut counts_g1 = Vec::new();

    let (m0, m1) = sample_info.get_markers();

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
        if !chunk.is_empty() {
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
    
    let u = min_by(u0, u1, |a, b| a.total_cmp(b));

    let s = ((n0 * n1 * (n + 1.) / 12.) - (n0 * n1 * tie_factor / (12. * n * (n - 1.)))).sqrt();
    //let s = ((n0 * n1 / 12.) * ((n + 1.) - (tie_factor / n * (n - 1.)))).sqrt();
    // supposedly same term but returns NaN???

    let z = (u - m) / s;

    let dist = Normal::standard();
    let p_value =  2. * (1. - dist.cdf(z.abs())) as f32;

    Ok(p_value)
}

// calculate the log2 of the log change of the two groups
fn log2_fold_change(tags: Tags, counts: Vec<u32>, sample_info: &SampleInfo) -> f32 {
    let mut norm_count_g0 = 0.;
    let mut norm_count_g1 = 0.;

    let (m0, m1) = sample_info.get_markers();

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

/// Trait for summarizing k-mers, determines the data saved in the graph nodes
pub trait SummaryData<DI> {
    /// format the noda data 
    fn print(&self, tag_translator: &BiMap<String, DI>, config: &SummaryConfig) -> String;
    /// format the noda data in one line
    fn print_ol(&self, tag_translator: &BiMap<String, DI>, config: &SummaryConfig) -> String;
    /// get `Tags` and the overall count, returns `None` if data is insufficient
    fn tags_sum(&self) -> Option<(Tags, u32)>;
    /// get "score" (the sum of the kmer appearances), `Vec<D>` simply returns `1.`
    fn score(&self) -> f32;
    /// get the size of the structure, including contents of boxed slices
    fn mem(&self) -> usize;
    /// get the number of observations, returns `None` if data is insufficient
    fn count(&self) -> Option<usize>;
    /// get the p-value, returns `None` if data is insufficient
    fn p_value(&self, config: &SummaryConfig) -> Option<f32>;
    /// get the log2(fold change), returns `None` if data is insufficient
    fn fold_change(&self, config: &SummaryConfig) -> Option<f32>;
    /// get the number of samples the sequence was observed in, returns `None` if data is insufficient
    fn sample_count(&self) -> Option<usize>;
    /// get edge multiplicities
    fn edge_mults(&self) -> Option<&EdgeMult>;
    /// fix the [`EdgeMult`] by removing hanging edges
    fn fix_edge_mults(&mut self, exts: Exts);
    /// check if node is valid according to: min kmer obs, group fraction, p-value
    fn valid(&self, config: &SummaryConfig) -> bool;
    /// summarize k-mers
    fn summarize<K, F: Iterator<Item = (K, Exts, DI)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self);
}
// TODO: move SummaryData::print functionality to Display trait?

/// Number of observations for the k-mer
impl<DI> SummaryData<DI> for u32 {
    fn print(&self, _: &BiMap<String, DI>, _: &SummaryConfig) -> String {
        format!("count: {}", self).replace("\"", "\'")
    }

    fn print_ol(&self, _: &BiMap<String, DI>, _: &SummaryConfig) -> String {
        format!("count: {}", self).replace("\"", "\'")
    }

    fn tags_sum(&self) -> Option<(Tags, u32)> { None }

    fn score(&self) -> f32 {
        *self as f32
    }

    fn mem(&self) -> usize {
        mem::align_of_val(self)
    }

    fn count(&self) -> Option<usize> {
        Some(*self as usize)
    }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn sample_count(&self) -> Option<usize> { None }

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn fix_edge_mults(&mut self, _: Exts) { }

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

/// data the k-mer was observed with
impl<DI: Debug + Ord + std::hash::Hash> SummaryData<DI> for Vec<DI> {
    fn print(&self, tag_translator: &BiMap<String, DI>, _: &SummaryConfig) -> String {
        let samples = self.iter().map(|sample_id| tag_translator.get_by_right(sample_id)).collect::<Vec<_>>();
        format!("samples: {:?}", samples).replace("\"", "\'")
    }

    fn print_ol(&self, tag_translator: &BiMap<String, DI>, _: &SummaryConfig) -> String {
        let samples = self.iter().map(|sample_id| tag_translator.get_by_right(sample_id)).collect::<Vec<_>>();
        format!("samples: {:?}", samples).replace("\"", "\'")
    }

    fn tags_sum(&self) -> Option<(Tags, u32)> { None }

    fn score(&self) -> f32 {
        1.
    }

    fn mem(&self) -> usize {
        mem::size_of_val(&**self) + mem::size_of_val(self)
    }

    fn count(&self) -> Option<usize> { None }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn sample_count(&self) -> Option<usize> {
        Some(self.len())
    }

    fn edge_mults(&self) -> Option<&EdgeMult> { None }
    
    fn fix_edge_mults(&mut self, _: Exts) { }

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

/// the u8-labels the k-mer was observed with and its number of observations
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
// aligned would be 16 Bytes, packed would be 12 Bytes
pub struct TagsSumData {
    tags: Tags,
    sum: u32,
}

impl SummaryData<u8> for TagsSumData {
    fn print(&self, tag_translator: &BiMap<String, u8>, _: &SummaryConfig) -> String {
        // replace " with ' to avoid conflicts in dot file
        format!("{}sum: {}", TagsFormatter::new(self.tags, tag_translator), self.sum).replace("\"", "\'")
    }

    fn print_ol(&self, tag_translator: &BiMap<String, u8>, _: &SummaryConfig) -> String {
        // replace " with ' to avoid conflicts in dot file
        format!("samples: {:?}, sum: {}", self.tags.to_string_vec(tag_translator), self.sum).replace("\"", "\'")
    }

    fn tags_sum(&self) -> Option<(Tags, u32)> {
        Some((self.tags, self.sum))
    }

    fn score(&self) -> f32 {
        self.sum as f32
    }

    fn mem(&self) -> usize {
        mem::size_of_val(self)
    }

    fn count(&self) -> Option<usize> {
        Some(self.sum as usize)
    }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn sample_count(&self) -> Option<usize> {
        Some(self.tags.len())
    }

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn fix_edge_mults(&mut self, _: Exts) { }

    fn valid(&self, config: &SummaryConfig) -> bool {
        valid_counts(self.tags, self.sum, config)
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, u8)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let mut all_exts = Exts::empty();

        let mut out_data: Vec<u8> = Vec::with_capacity(items.size_hint().0);

        let mut sum = 0u32;
        for (_, exts, d) in items {
            out_data.push(d); 
            all_exts = all_exts.add(exts);
            sum += 1;
        }

        out_data.sort();
        out_data.dedup();

        let tags = Tags::from_u8_vec(out_data);
        
        (valid_counts(tags, sum, config), all_exts, TagsSumData { tags, sum })
    }

}

/// Implementation of [`SummaryData<u8>`]
/// 
/// Contains the `u8`-labels the k-mer was observed with, how many times it 
/// was observed with each label, and how many times it was observed overall
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TagsCountsSumData {
    tags: Tags,
    sum: u32,
    counts: Box<[u32]>,
}

impl SummaryData<u8> for TagsCountsSumData {
    fn print(&self, tag_translator: &BiMap<String, u8>, config: &SummaryConfig) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        format!("{}sum: {}{}{}", TagsCountsFormatter::new(self.tags, &self.counts, tag_translator), self.sum, p, fc).replace("\"", "\'")
    }

    fn print_ol(&self, tag_translator: &BiMap<String, u8>, config: &SummaryConfig) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        format!("samples: {:?}, counts: {:?}, sum: {}{}{}", self.tags.to_string_vec(tag_translator), self.counts, self.sum, p, fc).replace("\"", "\'")
    }

    fn tags_sum(&self) -> Option<(Tags, u32)> {
        Some((self.tags, self.sum))
    }

    fn score(&self) -> f32 {
        self.sum as f32
    }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.counts)
    }

    fn count(&self) -> Option<usize> {
        Some(self.sum as usize)
    }

    fn sample_count(&self) -> Option<usize> {
        Some(self.counts.len())
    }

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn fix_edge_mults(&mut self, _: Exts) { }

    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {      
        match p_value(&self.tags.to_u8_vec(), &self.counts.to_vec(), config) {
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

        valid_counts(self.tags, self.sum, config) && valid_p
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, u8)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let (all_exts, out_data, tag_counts, sum) = summarize(items);

        let valid_p = valid_p(PInfo::Calculate { out_data: &out_data, tag_counts: &tag_counts}, config);     

        let counts: Box<[u32]> = tag_counts.into();
        let tags = Tags::from_u8_vec(out_data);

        (valid_counts(tags, sum, config) && valid_p, all_exts, TagsCountsSumData { tags, counts, sum }) 
    }

}

/// Implementation of [`SummaryData<u8>`]
/// 
/// Contains the `u8`-labels the k-mer was observed with and how many times it 
/// was observed with each label
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TagsCountsData {
    tags: Tags,
    counts: Box<[u32]>
}

impl TagsCountsData {
    #[inline]
    pub fn sum(&self) -> u32 {
        self.counts.iter().sum::<u32>()
    }
}

impl SummaryData<u8> for TagsCountsData {
    fn print(&self, tag_translator: &BiMap<String, u8>, config: &SummaryConfig) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        format!("{}sum: {}{}{}", TagsCountsFormatter::new(self.tags, &self.counts, tag_translator), self.sum(), p, fc).replace("\"", "\'")
    }

    fn print_ol(&self, tag_translator: &BiMap<String, u8>, config: &SummaryConfig) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        format!("samples: {:?}, counts: {:?}, sum: {}{}{}", self.tags.to_string_vec(tag_translator), self.counts, self.sum(), p, fc).replace("\"", "\'")
    }

    fn tags_sum(&self) -> Option<(Tags, u32)> {
        Some((self.tags, self.sum()))
    }

    fn score(&self) -> f32 {
        self.sum() as f32
    }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.counts)
    }

    fn count(&self) -> Option<usize> {
        Some(self.counts.iter().sum::<u32>() as usize)
    }

    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {      
        match p_value(&self.tags.to_u8_vec(), &self.counts.to_vec(), config) {
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

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn fix_edge_mults(&mut self, _: Exts) { }

    fn valid(&self, config: &SummaryConfig) -> bool {
        let valid_p = match config.max_p {
            Some(p) => self.p_value(config).expect("error calculating p-value") <= p,
            None => true,
        }; 

        valid_counts(self.tags, self.sum(), config) && valid_p
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, u8)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let (all_exts, out_data, tag_counts, sum) = summarize(items);

        let valid_p = valid_p(PInfo::Calculate { out_data: &out_data, tag_counts: &tag_counts}, config);            

        let counts: Box<[u32]> = tag_counts.into();
        let tags = Tags::from_u8_vec(out_data);

        (valid_counts(tags, sum, config) && valid_p, all_exts, TagsCountsData { tags, counts }) 
    }

}

/// Implementation of [`SummaryData<u8>`]
/// 
/// Contains the `u8`-labels the k-mer was observed with, how many times it 
/// was observed with each label, and a p-value
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TagsCountsPData {
    tags: Tags,
    counts: Box<[u32]>,
    p_value: f32
}

impl TagsCountsPData {
    #[inline]
    pub fn sum(&self) -> u32 {
        self.counts.iter().sum::<u32>()
    }
}

impl SummaryData<u8> for TagsCountsPData {
    fn print(&self, tag_translator: &BiMap<String, u8>, config: &SummaryConfig) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        format!("{}sum: {}{}{}", TagsCountsFormatter::new(self.tags, &self.counts, tag_translator), self.sum(), p, fc).replace("\"", "\'")
    }

    fn print_ol(&self, tag_translator: &BiMap<String, u8>, config: &SummaryConfig) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        format!("samples: {:?}, counts: {:?}, sum: {}{}{}", self.tags.to_string_vec(tag_translator), self.counts, self.sum(), p, fc).replace("\"", "\'")
    }

    fn tags_sum(&self) -> Option<(Tags, u32)> {
        Some((self.tags, self.sum()))
    }

    fn score(&self) -> f32 {
        self.sum() as f32
    }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.counts)
    }

    fn count(&self) -> Option<usize> {
        Some(self.sum() as usize)
    }

    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {
        if config.stat_test_changed {
            Some(p_value(&self.tags.to_u8_vec(), &self.counts.to_vec(), config).unwrap())
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

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn fix_edge_mults(&mut self, _: Exts) { }

    fn valid(&self, config: &SummaryConfig) -> bool {
        valid_counts(self.tags, self.sum(), config) 
            && valid_p(PInfo::PValue { p: self.p_value(config).expect("error getting p-values") }, config)
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, u8)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let (all_exts, out_data, tag_counts, sum) = summarize(items);

        // caluclate p-value
        let p_value = p_value(&out_data, &tag_counts, config).unwrap();

        let counts: Box<[u32]> = tag_counts.into();
        let tags = Tags::from_u8_vec(out_data);

        let valid = valid_counts(tags, sum, config) && valid_p(PInfo::PValue { p: p_value }, config);

        (valid, all_exts, TagsCountsPData { tags, counts, p_value }) 
    }

}

/// Implementation of [`SummaryData<u8>`]
/// 
/// Contains the `u8`-labels the k-mer was observed with, how many times it 
/// was observed with each label, and the edge multiplicites
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TagsCountsEMData {
    tags: Tags,
    counts: Box<[u32]>,
    edge_mults: EdgeMult,
}

impl TagsCountsEMData {
    #[inline]
    pub fn sum(&self) -> u32 {
        self.counts.iter().sum::<u32>()
    }
}

impl SummaryData<u8> for TagsCountsEMData {
    fn print(&self, tag_translator: &BiMap<String, u8>, config: &SummaryConfig) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        format!("{}sum: {}{}{}, edge multiplicities: \n{}", TagsCountsFormatter::new(self.tags, &self.counts, tag_translator), self.sum(), p, fc, self.edge_mults).replace("\"", "\'")
    }

    fn print_ol(&self, tag_translator: &BiMap<String, u8>, config: &SummaryConfig) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        format!("samples: {:?}, counts: {:?}, sum: {}{}{}, edge multiplicities: {:?}", self.tags.to_string_vec(tag_translator), self.counts, self.sum(), p, fc, self.edge_mults).replace("\"", "\'")
    }

    fn tags_sum(&self) -> Option<(Tags, u32)> {
        Some((self.tags, self.sum()))
    }

    fn score(&self) -> f32 {
        self.sum() as f32
    }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.counts)
    }

    fn count(&self) -> Option<usize> {
        Some(self.counts.iter().sum::<u32>() as usize)
    }

    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {      
        match p_value(&self.tags.to_u8_vec(), &self.counts.to_vec(), config) {
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

    fn edge_mults(&self) -> Option<&EdgeMult> {
        Some(&self.edge_mults)
    }

    fn fix_edge_mults(&mut self, exts: Exts) {
        self.edge_mults.clean_edges(exts);
    }


    fn valid(&self, config: &SummaryConfig) -> bool {
        let valid_p = match config.max_p {
            Some(p) => self.p_value(config).expect("error calculating p-value") <= p,
            None => true,
        }; 

        valid_counts(self.tags, self.sum(), config) && valid_p
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, u8)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let (all_exts, out_data, tag_counts, sum, edge_mults) = summarize_with_em(items);

        let valid_p = valid_p(PInfo::Calculate { out_data: &out_data, tag_counts: &tag_counts}, config);            

        let counts: Box<[u32]> = tag_counts.into();
        let tags = Tags::from_u8_vec(out_data);

        (valid_counts(tags, sum, config) && valid_p, all_exts, TagsCountsEMData { tags, counts, edge_mults }) 
    }

}

/// Implementation of [`SummaryData<u8>`]
/// 
/// Contains the `u8`-labels the k-mer was observed with, how many times it 
/// was observed with each label, a p-value, and the edge multiplicites
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TagsCountsPEMData {
    tags: Tags,
    counts: Box<[u32]>,
    p_value: f32,
    edge_mults: EdgeMult,
}

impl TagsCountsPEMData {
    #[inline]
    pub fn sum(&self) -> u32 {
        self.counts.iter().sum::<u32>()
    }
}

impl SummaryData<u8> for TagsCountsPEMData{
    fn print(&self, tag_translator: &BiMap<String, u8>, config: &SummaryConfig) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        format!("{}sum: {}{}{}, edge multiplicities: \n{}", TagsCountsFormatter::new(self.tags, &self.counts, tag_translator), self.sum(), p, fc, self.edge_mults).replace("\"", "\'")
    }

    fn print_ol(&self, tag_translator: &BiMap<String, u8>, config: &SummaryConfig) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        format!("samples: {:?}, counts: {:?}, sum: {}{}{}, edge multiplicities: {:?}", self.tags.to_string_vec(tag_translator), self.counts, self.sum(), p, fc, self.edge_mults).replace("\"", "\'")
    }

    fn tags_sum(&self) -> Option<(Tags, u32)> {
        Some((self.tags, self.sum()))
    }

    fn score(&self) -> f32 {
        self.sum() as f32
    }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.counts)
    }

    fn count(&self) -> Option<usize> {
        Some(self.counts.iter().sum::<u32>() as usize)
    }

    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {
        if config.stat_test_changed {
            Some(p_value(&self.tags.to_u8_vec(), &self.counts.to_vec(), config).unwrap())
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

    fn edge_mults(&self) -> Option<&EdgeMult> {
        Some(&self.edge_mults)
    }

    fn fix_edge_mults(&mut self, exts: Exts) {
        self.edge_mults.clean_edges(exts);
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        valid_counts(self.tags, self.sum(), config) 
            && valid_p(PInfo::PValue { p: self.p_value(config).expect("error getting p-values") }, config)
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, u8)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let (all_exts, out_data, tag_counts, sum, edge_mults) = summarize_with_em(items);

        // caluclate p-value with chosen test
        let p_value = p_value(&out_data, &tag_counts, config).unwrap();         

        let counts: Box<[u32]> = tag_counts.into();
        let tags = Tags::from_u8_vec(out_data);

        let valid = valid_counts(tags, sum, config) && valid_p(PInfo::PValue { p: p_value }, config);

        (valid, all_exts, TagsCountsPEMData { tags, counts, p_value, edge_mults }) 
    }

}

/// Implementation of [`SummaryData<u8>`]
/// 
/// Contains how many times the k-mer was observed in each group
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

impl SummaryData<u8> for GroupCountData {
    fn print(&self, _: &BiMap<String, u8>, _: &SummaryConfig) -> String {
        format!("count 1: {}\ncount 2: {}", self.group1, self.group2)
    }

    fn print_ol(&self, _: &BiMap<String, u8>, _: &SummaryConfig) -> String {
        format!("count 1: {}, count 2: {}", self.group1, self.group2)
    }

    fn tags_sum(&self) -> Option<(Tags, u32)> { None }

    fn score(&self) -> f32 {
        self.sum() as f32
    }

    fn mem(&self) -> usize {
        mem::size_of_val(self)
    }

    fn count(&self) -> Option<usize> {
        Some((self.group1 + self.group2) as usize)
    }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn sample_count(&self) -> Option<usize> { None }

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn fix_edge_mults(&mut self, _: Exts) { }

    fn valid(&self, config: &SummaryConfig) -> bool {
        self.sum() >= config.min_kmer_obs as u32
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, u8)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let mut all_exts = Exts::empty();
        let mut count1 = 0;
        let mut count2 = 0;

        let mut nobs = 0u32;
        for (_, exts, d) in items {
            let tag = (2 as M).pow(d as u32);
            let group1 = ((config.sample_info.marker0 & tag) > 0) as u32;
            let group2 = ((config.sample_info.marker1 & tag) > 0) as u32;

            if (group1 + group2) != 1 { 
                panic!(
                    "should not happen\n tag: {:#066b}\n m1:  {:#066b}\n m2:  {:#066b}\n g1:  {}\n g2:  {}", 
                    tag, config.sample_info.marker0, config.sample_info.marker1, group1, group2
                )
            }
            count1 += group1;
            count2 += group2;
            nobs += 1;
            all_exts = all_exts.add(exts);
        }

        let (group1, group2) = match config.significant {
            Some(digits) => (round_digits(count1, digits), round_digits(count2, digits)),
            None => (count1, count2)
        };

        assert_eq!((count1 + count2),nobs);
        (nobs as usize >= config.min_kmer_obs, all_exts, GroupCountData { group1, group2 })
    }

}

/// Implementation of [`SummaryData<u8>`]
/// 
/// Contains the relative number of observations for the k-mer (in percent) 
/// in group 1 and the absolute overall count 
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RelCountData {
    percent: u32,
    count: u32
}

impl SummaryData<u8> for RelCountData {    
    fn print(&self, _: &BiMap<String, u8>, _: &SummaryConfig) -> String {
        format!("relative amount group 1: {}\ncount both: {}", self.percent, self.count)
    }

    fn print_ol(&self, _: &BiMap<String, u8>, _: &SummaryConfig) -> String {
        format!("relative amount group 1: {}, count both: {}", self.percent, self.count)
    }

    fn tags_sum(&self) -> Option<(Tags, u32)> { None }

    fn score(&self) -> f32 {
        self.count as f32
    }

    fn mem(&self) -> usize {
        mem::size_of_val(self)
    }

    fn count(&self) -> Option<usize> {
        Some(self.count as usize)
    }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn sample_count(&self) -> Option<usize> { None }

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn fix_edge_mults(&mut self, _: Exts) { }

    fn valid(&self, config: &SummaryConfig) -> bool {
        self.count >= config.min_kmer_obs as u32
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, u8)>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let mut all_exts = Exts::empty();
        let mut count1 = 0;
        let mut count2 = 0;

        let mut nobs = 0u32;
        for (_, exts, d) in items {
            let tag = (2 as M).pow(d as u32);
            let group1 = ((config.sample_info.marker0 & tag) > 0) as u32;
            let group2 = ((config.sample_info.marker1 & tag) > 0) as u32;

            if (group1 + group2) != 1 { 
                panic!(
                    "should not happen\n tag: {:#066b}\n m1:  {:#066b}\n m2:  {:#066b}\n g1:  {}\n g2:  {}", 
                    tag, config.sample_info.marker0, config.sample_info.marker1, group1, group2
                )
            }
            count1 += group1;
            count2 += group2;
            all_exts = all_exts.add(exts);
            nobs += 1;
        }

        assert_eq!(count1 + count2, nobs);


        let percent = (count1 as f64 / nobs as f64 * 100.) as u32;
        let count = match config.significant {
            Some(digits) => round_digits(nobs, digits),
            None => nobs
        };

        (nobs as usize >= config.min_kmer_obs, all_exts, RelCountData { percent, count }) 
    }
    
}

#[cfg(test)]
mod test {
    use std::{fs::File, io::BufReader};
    

    use crate::{clean_graph::CleanGraph, compression::{ compress_graph, ScmapCompress}, dna_string::DnaString, graph::{BaseGraph, DebruijnGraph, Node}, kmer::{Kmer16, Kmer8}, summarizer::{self, p_value, students_t_test, u_test, valid_p, welchs_t_test, GroupFrac, NotEnoughSamplesError, SampleInfo, SummaryData}, Exts, Tags};

    use super::{log2_fold_change, round_digits, SummaryConfig, TagsCountsSumData};

    #[test]
    fn test_p_value() -> Result<(), NotEnoughSamplesError> {

        /*
        group 1: 111111100000 = 4064
        group 2: 000000011111 = 31
         */

        let sample_kmers = vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
        let sample_info = SampleInfo::new(31, 4064, 5, 7, sample_kmers);
        let summary_config_w = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), None, summarizer::StatTest::WelchsTTest);
        let summary_config_t = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), None, summarizer::StatTest::StudentsTTest);
        let summary_config_u = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), None, summarizer::StatTest::UTest);

        let out_data = [0, 1, 2, 3, 4, 8];
        let tag_counts = vec![1; 6];
        
        let p = welchs_t_test(&out_data, &tag_counts, &sample_info)?;
        assert_eq!((p * 1000.).round(), 1.);
        let p = students_t_test(&out_data, &tag_counts, &sample_info)?;
        assert_eq!((p * 10000.).round(), 5.);
        let p = u_test(&out_data, &tag_counts, &sample_info)?;
        assert_eq!((p * 1000.).round(), 5.);

        let p: f32 = p_value(&out_data, &tag_counts, &summary_config_w)?;
        assert_eq!((p * 1000.).round(), 1.);
        let p = p_value(&out_data, &tag_counts, &summary_config_t)?;
        assert_eq!((p * 10000.).round(), 5.);
        let p = p_value(&out_data, &tag_counts, &summary_config_u)?;
        assert_eq!((p * 1000.).round(), 5.);


        // test with different kmer counts
        let sample_kmers = vec![12, 3345, 3478, 87, 1, 2, 666, 98111, 23982938, 555, 122, 7238];

        let sample_info = SampleInfo::new(31, 4064, 5, 7, sample_kmers);
        let summary_config_w = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), None, summarizer::StatTest::WelchsTTest);
        let summary_config_t = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), None, summarizer::StatTest::StudentsTTest);
        let summary_config_u = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), None, summarizer::StatTest::UTest);

        let out_data = [0, 1, 7, 8, 9, 10];
        let tag_counts = vec![1; 6];

        let p = welchs_t_test(&out_data, &tag_counts, &sample_info)?;
        assert_eq!((p * 10000.).round(), 4113.);
        let p = students_t_test(&out_data, &tag_counts, &sample_info)?;
        assert_eq!((p * 10000.).round(), 2955.);
        let p = u_test(&out_data, &tag_counts, &sample_info)?;
        assert_eq!((p * 1000.).round(), 862.);

        let p: f32 = p_value(&out_data, &tag_counts, &summary_config_w)?;
        assert_eq!((p * 10000.).round(), 4113.);
        let p: f32 = p_value(&out_data, &tag_counts, &summary_config_t)?;
        assert_eq!((p * 10000.).round(), 2955.);
        let p: f32 = p_value(&out_data, &tag_counts, &summary_config_u)?;
        assert_eq!((p * 1000.).round(), 862.);


        // test: not enough samples for statistical test

        /*
        group 1: 000000100 = 4
        group 2: 000000011 = 3
         */

        let sample_kmers = vec![3, 3, 3];
        let sample_info = SampleInfo::new(4, 4, 1, 2, sample_kmers);
        let summary_config_w = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), None, summarizer::StatTest::WelchsTTest);
        let summary_config_t = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), None, summarizer::StatTest::StudentsTTest);
        let summary_config_u = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), None, summarizer::StatTest::UTest);

        let out_data = [0, 1, 2];
        let tag_counts = vec![2, 3, 1];

        if welchs_t_test(&out_data, &tag_counts, &sample_info).is_ok() { panic!("should throw err") }
        if students_t_test(&out_data, &tag_counts, &sample_info).is_ok() { panic!("should throw err") }
        if u_test(&out_data, &tag_counts, &sample_info).is_ok() { panic!("should throw err") }

        if p_value(&out_data, &tag_counts, &summary_config_w).is_ok() { panic!("should throw err") }
        if p_value(&out_data, &tag_counts, &summary_config_t).is_ok() { panic!("should throw err") }
        if p_value(&out_data, &tag_counts, &summary_config_u).is_ok() { panic!("should throw err") }

        Ok(())
    }

    #[test]
    fn test_valid_p() {
        let sample_kmers = vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
        let sample_info = SampleInfo::new(31, 4064, 5, 7, sample_kmers);
        let summary_config_m = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), None, summarizer::StatTest::WelchsTTest);
        let summary_config_p = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), Some(0.05), summarizer::StatTest::WelchsTTest);

        let out_data = [0, 1, 2, 3, 4, 8];  // p should be 0.001
        let tag_counts = vec![1; 6];
        let vp = valid_p(summarizer::PInfo::Calculate { out_data: &out_data, tag_counts: &tag_counts }, &summary_config_m);
        assert!(vp);
        let vp = valid_p(summarizer::PInfo::Calculate { out_data: &out_data, tag_counts: &tag_counts }, &summary_config_p);
        assert!(vp);


        let out_data = [0, 7, 8, 9, 10]; // p should be 0.2238
        let tag_counts = vec![1; 5];

        let vp = valid_p(summarizer::PInfo::Calculate { out_data: &out_data, tag_counts: &tag_counts }, &summary_config_m);
        assert!(vp);
        let vp = valid_p(summarizer::PInfo::Calculate { out_data: &out_data, tag_counts: &tag_counts }, &summary_config_p);
        assert!(!vp);
    }

    #[test]
    #[cfg(not(feature = "sample128"))]
    fn test_data_valid() {
        use crate::BUF;

        let mut graph: BaseGraph<Kmer8, TagsCountsSumData> = BaseGraph::new(false);

        let tags = Tags::from_u8_vec(vec![0, 2, 6]);
        let counts: Box<[u32]> = [1, 3, 5].into();
        let sum = counts.iter().sum::<u32>();
        graph.add(&DnaString::from_acgt_bytes("AAAAAAAA".as_bytes()), Exts::empty(), TagsCountsSumData { tags, counts, sum });

        let tags = Tags::from_u8_vec(vec![0]);
        let counts: Box<[u32]> = [1].into();
        let sum = counts.iter().sum::<u32>();
        graph.add(&DnaString::from_acgt_bytes("CCCCCCCC".as_bytes()), Exts::empty(), TagsCountsSumData { tags, counts, sum });
                
        
        let graph = graph.finish();

        graph.print();

        let sample_kmers = vec![123, 234, 12334, 34, 1232, 123, 123, 34];
        let sample_info = SampleInfo::new(0b00100101, 0b11011010, 3, 5, sample_kmers);
        let config = SummaryConfig::new(3, None, GroupFrac::None, 0.33,  sample_info, None, summarizer::StatTest::StudentsTTest);

        let censor_nodes = CleanGraph::new(|node: &Node<'_, Kmer8, TagsCountsSumData>| !node.data().valid(&config))
                    .find_bad_nodes(&graph);
        println!("censor nodes: {:?}", censor_nodes);
        let filter_graph = compress_graph(false, &ScmapCompress::new(), graph, Some(censor_nodes));

        filter_graph.print();

        // larger test

        let reader = BufReader::with_capacity(BUF, File::open("test_data/400.graph.dbg").unwrap());
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
        //let summary_config = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), None, summarizer::StatTest::WelchsTTest);

        let labels = vec![0, 1, 2, 3, 7, 8];
        let tags = Tags::from_u8_vec(labels);
        let counts = vec![1, 6, 9, 3, 6, 10];
        let fold_change = log2_fold_change(tags, counts, &sample_info);
        assert_eq!(fold_change, 5.286_453_2);

        let labels = vec![0, 6, 7, 8, 10, 11];
        let tags = Tags::from_u8_vec(labels);
        let counts = vec![12, 3, 7, 1, 22, 6];
        let fold_change = log2_fold_change(tags, counts, &sample_info);
        assert_eq!(fold_change, -1.339_324_5);

        // x/0 = inf -> log(inf) = inf
        let labels = vec![0, 1];
        let tags = Tags::from_u8_vec(labels);
        let counts = vec![12, 3];
        let fold_change = log2_fold_change(tags, counts, &sample_info);
        assert_eq!(fold_change, f32::INFINITY);

        // 0/x = 0 -> log2(0) = -inf
        let labels = vec![7, 8];
        let tags = Tags::from_u8_vec(labels);
        let counts = vec![12, 3];
        let fold_change = log2_fold_change(tags, counts, &sample_info);
        assert_eq!(fold_change, f32::NEG_INFINITY);       
    }

    #[test]
    fn test_round_digits() {
        assert_eq!(round_digits(18293092, 4), 18290000);
        assert_eq!(round_digits(333, 4), 333);
        assert_eq!(round_digits(129552, 2), 130000);
        assert_eq!(round_digits(1829399, 6), 1829400);
        assert_eq!(round_digits(1829399, 0), 0);
    }


}




