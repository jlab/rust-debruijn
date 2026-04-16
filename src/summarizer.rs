use bimap::BiMap;
use clap::ValueEnum;
use serde::{de::DeserializeOwned, Deserialize, Serialize};
use statrs::distribution::{ContinuousCDF, Normal, StudentsT};
use crate::{BaseQuality, EdgeMult, Exts, Kmer, KmerDataItem, Tags, TagsCountsFormatter, TagsFormatter};
use std::{cmp::min_by, collections::HashMap, error::Error, fmt::{Debug, Display}, mem};

/// inner type for [`Tags`] and group markers
#[cfg(not(feature = "sample128"))]
pub type Marker = u64;

/// inner type for [`Tags`] and group markers
#[cfg(feature = "sample128")]
pub type Marker = u128;

/// type for IDs (e.g. gene IDs)
#[cfg(not(feature = "id4b"))]
pub type ID = u16;

/// type for IDs (e.g. gene IDs)
#[cfg(feature = "id4b")]
pub type ID = u32;

/// type for tags
pub type Tag = u8;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy, Serialize, Deserialize, Hash)]
/// type for IDs and tags together
pub struct IDTag {
    id: ID,
    tag: Tag
}

impl IDTag {
    pub fn new(id: ID, tag: Tag) -> IDTag {
        IDTag { id, tag }
    }

    pub fn tag(&self) -> Tag {
        self.tag
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
/// translate tags and IDs (e.g. into sample labels and gene names)
pub struct Translator {
    ids: Option<BiMap<String, ID>>,
    tags: Option<BiMap<String,Tag>>,
}

impl Translator {
    /// make a new [`Translator`] for tags and IDs
    pub fn new(ids: BiMap<String, ID>, tags: BiMap<String,Tag>) -> Translator {
        Translator { ids: Some(ids), tags: Some(tags) }
    }

    /// make an empty [`Translator`]
    pub fn empty() -> Translator {
        Translator { ids: None, tags: None }
    }

    /// make a new [`Translator`] for tags
    pub fn new_tag_translator(hashed_tags: BiMap<String, Tag>) -> Translator {
        Translator { ids: None, tags: Some(hashed_tags) }
    }

    /// make a new [`Translator`] for IDs
    pub fn new_id_translator(hashed_ids: BiMap<String, ID>) -> Translator {
        Translator { ids: Some(hashed_ids), tags: None }
    }

    /// get the tag translator, returns None if the `Translator` does not contain a tag translator
    pub fn tag_translator(&self) -> &Option<BiMap<String, Tag>> {
        &self.tags
    }

    /// get the tag translator, returns None if the `Translator` does not contain a tag translator
    pub fn id_translator(&self) -> &Option<BiMap<String, ID>> {
        &self.ids
    }

    /// get a mutable reference to the tag translator, returns None if the `Translator` does not contain a tag translator
    pub fn mut_id_translator(&mut self) -> &mut Option<BiMap<String, ID>> {
        &mut self.ids
    }

    /// dissolve the `Translator` into its underlying [`BiMap`]s
    pub fn dissolve(self) -> (Option<BiMap<String, ID>>, Option<BiMap<String,Tag>>) {
        (self.ids, self.tags)
    }
}

fn id_format(ids: &[ID], translator: &Translator, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
    if let Some(id_gr_tr) = id_group_translator {
        // translate the ids (genes) into id groups (orthogroups)
        let mut t_ids = ids
            .iter()
            .map(|id| id_gr_tr.get(id).unwrap_or_else(|| panic!("ID does not exist - ids {:?}", ids)))
            .collect::<Vec<_>>();
        t_ids.sort();
        t_ids.dedup();
        format!("{:?}", t_ids)
    } else if let Some(id_translator) = translator.id_translator() {
        // translate the ids (genes) into their names
        let t_ids = ids
            .iter()
            .map(|id| id_translator.get_by_right(id).unwrap_or_else(|| panic!("ID does not exist - ids {:?}", ids)))
            .collect::<Vec<_>>();
        format!("{:?}", t_ids)
    } else {
        // do not translate
        format!("{:?}", ids)
    }
}

#[derive(Debug, PartialEq)]
struct NotEnoughSamplesError {}

impl Error for NotEnoughSamplesError {}

impl Display for NotEnoughSamplesError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "not enough samples were supplied to perform a statistical test")
    }
}

/// Configuration for summary processes. It it used to filter the k-mers during
/// and after graph construction and to make statistical analyses based on k-mer
/// occurrence in the sample groups. 
/// 
/// For any statistical analyses regarding the sample groups, the `SummaryData` 
/// requires a [`SampleInfo`].
/// 
/// The available options are:
/// - filter the k-mers by 
///     - their number of occurrences
///     - their quality based on the phred scores from the reads
///     - their p-values regarding sample groups
///     - if they occurr in at least a specific fraction of one or both of
///       the sample groups 
/// - filter the k-mer occurrences by the quality - should this disconnect the 
///   k-mer, alll occurrences will be used
/// - set the number of occurrences which is stored with the k-mers to be rounded
///   to a number of significant digits
/// - choose a statistical test on which the p-value calculation is based, by default 
///   this is set to Welch's t-test; Student's t-test and the Mann-Whitney U test 
///   are also available
/// 
/// By default, all filter options are turned off.
/// 
/// ```
/// use debruijn::summarizer::{GroupFrac, SummaryConfig, SampleInfo, StatTest};
/// use debruijn::BaseQuality;
/// 
/// let sample_info = SampleInfo::new(0b1100, 0b0011, vec![100, 100, 100, 100]);
/// let summary_config = SummaryConfig::new(sample_info.clone())
///     .with_min_kmer_obs(2)
///     .with_min_quality(BaseQuality::Medium)
///     .with_max_p(Some(0.05)) // can also be none which can avoid p-value calculations and save time
///     .with_group_frac(GroupFrac::One, 0.3);
/// 
/// let summary_config2 = SummaryConfig::new(sample_info.clone())
///     .with_min_kmer_obs(4)
///     .with_min_quality_for_edge(BaseQuality::Marginal);
/// 
/// let summary_config3 = SummaryConfig::new(sample_info)
///     .with_significant(Some(5))
///     .with_stat_test(StatTest::StudentsTTest);
/// ```
/// 
/// To filter an already constructed graph, please use the original `SummaryConfig`
/// and change the settings in plase with the `set_...` methods. This is so the program
/// knows that some settings have been changed and corresponding values have
/// to be re-calculated.
/// 
/// ```
/// use debruijn::summarizer::{SummaryConfig, SampleInfo, StatTest};
/// use debruijn::BaseQuality;
/// 
/// let sample_info = SampleInfo::new(0b1100, 0b0011, vec![100, 100, 100, 100]);
/// let mut summary_config = SummaryConfig::new(sample_info);
/// 
/// // ... construct the graph
/// 
/// summary_config.set_max_p(Some(0.05));
/// summary_config.set_stat_test(StatTest::StudentsTTest);
/// 
/// // ...
/// ```
#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct SummaryConfig {
    min_kmer_obs: usize,
    significant: Option<u32>,
    group_frac: GroupFrac,
    frac_cutoff: f32,
    sample_info: SampleInfo,
    max_p: Option<f32>,
    stat_test: StatTest,
    stat_test_changed: bool,
    min_quality: BaseQuality,
    min_quality_for_edge: BaseQuality
}

impl SummaryConfig {
    /// make a new `SummaryConfig`
    /// 
    /// arguments: 
    /// * `sample_info`: a [`SampleInfo`] with information about the sample groups,
    ///   which is required for any statistical analysis
    pub fn new(sample_info: SampleInfo) -> Self {
        SummaryConfig::empty().with_sample_info(sample_info)
    }

    /// make an empty `SummaryConfig`. A proper [`SampleInfo`] is required for any
    /// statistical analysis regarding sample groups.
    pub fn empty() -> Self {
        SummaryConfig { 
            min_kmer_obs: 0, 
            significant: None, 
            group_frac: GroupFrac::None, 
            frac_cutoff: 0., 
            sample_info: SampleInfo::empty(), 
            max_p: None, 
            stat_test: StatTest::WelchsTTest, 
            stat_test_changed: false,
            min_quality: BaseQuality::NoCall,
            min_quality_for_edge: BaseQuality::NoCall,
        }
    }

    /// produce a new `SummaryConfig` which will filter k-mers by their number 
    /// of observations
    pub fn with_min_kmer_obs(&self, min_kmer_obs: usize) -> Self {
        let mut config = self.clone();
        config.min_kmer_obs = min_kmer_obs;
        config
    }

    /// modify the number of k-mers observations required for each k-mer to be 
    /// included in the graph
    pub fn set_min_kmer_obs(&mut self, min_kmer_obs: usize) {
        self.min_kmer_obs = min_kmer_obs;
    }

    /// produce a new `SummaryConfig` which will round the number of observations
    /// of the k-mer to `significant_digits`
    pub fn with_significant(&self, significant_digits: Option<u32>) -> Self {
        let mut config = self.clone();
        config.significant = significant_digits;
        config
    }

    /// modify the number signigicant digits the number of observations stored 
    /// with the k-mer will be rounded to
    pub fn set_significant(&mut self, significant_digits: Option<u32>) {
        self.significant = significant_digits;
    }

    /// produce a new `SummaryConfig` which will require the k-mer to be ovserved
    /// in at least a fraction of `frac_cutoff` of either one, both or none of the
    /// sample groups
    pub fn with_group_frac(&self, group_frac: GroupFrac, frac_cutoff: f32) -> Self {
        let mut config = self.clone();
        config.group_frac = group_frac;
        config.frac_cutoff = frac_cutoff;
        config
    }

    /// modify the group fract settings which require the k-mer to be ovserved
    /// in at least a fraction of `frac_cutoff` of either one, both or none of the
    /// sample groups
    pub fn set_group_frac(&mut self, group_frac: GroupFrac, frac_cutoff: f32) {
        self.group_frac = group_frac;
        self.frac_cutoff = frac_cutoff;
    }

    /// produce a new `SummaryConfig` which contains the given [`SampleInfo`]
    fn with_sample_info(self, sample_info: SampleInfo) -> Self {
        let mut config = self.clone();
        config.sample_info = sample_info;
        config
    }

    /// produce a new `SummaryConfig` which will filter k-mers by their p-value
    /// regarding occurrence in the sample groups
    pub fn  with_max_p(&self, max_p: Option<f32>) -> Self {
        let mut config = self.clone();
        config.max_p = max_p;
        config
    }

    /// modify the maximum p-value a k-mer is allowed to have
    pub fn  set_max_p(&mut self, max_p: Option<f32>) {
        self.max_p = max_p;
    }

    /// produce a new `SummaryConfig` which will calculate the p-values based on
    /// the given statistical test
    pub fn with_stat_test(&self, stat_test: StatTest) -> Self {
        let mut config = self.clone();
        config.stat_test = stat_test;
        config
    }

    /// modify the statistical test, which is used to calculate p-values
    /// regarding observations in the two sample groups
    pub fn set_stat_test(&mut self, stat_test: StatTest) {
        if stat_test != self.stat_test { self.stat_test_changed = true }
        self.stat_test = stat_test;
    }

    /// produce a new `SummaryConfig` which will filter k-mers by their quality
    pub fn with_min_quality(&self, min_quality: BaseQuality) -> Self {
        let mut config = self.clone();
        config.min_quality = min_quality;
        config
    }

    /// modify the minimum quality required for each k-mer to be 
    /// included in the graph
    pub fn set_min_quality(&mut self, min_quality: BaseQuality) {
        self.min_quality = min_quality;
    }

    /// produce a new `SummaryConfig` which will keep k-mers from being counted 
    /// if it's quality is too low - should this disconnect the k-mer it will 
    /// still be counted
    pub fn with_min_quality_for_edge(&self, min_quality_for_edge: BaseQuality) -> Self {
        let mut config = self.clone();
        config.min_quality_for_edge = min_quality_for_edge;
        config
    }

    /// modify the minimum quality required for each k-mer observation for its
    /// edges to be counted
    pub fn set_min_quality_for_edge(&mut self, min_quality_for_edge: BaseQuality) {
        self.min_quality_for_edge = min_quality_for_edge;
    }

    /// get the binary encoded group affiliation of tags
    pub fn get_markers(&self) -> (Marker, Marker) {
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
#[serde(rename_all = "kebab-case")]
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
#[serde(rename_all = "kebab-case")]
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
/// use debruijn::summarizer::{SampleInfo, Marker};
/// 
/// let marker0: Marker = 0b0000111; // = 7
/// let marker1: Marker = 0b1111000; // = 120
/// 
/// let sample_kmers = vec![1232, 12323, 24342, 24234, 345456, 21234, 546456];
/// assert_eq!(marker0.count_ones() + marker1.count_ones(), sample_kmers.len() as u32);
/// 
/// let sample_info = SampleInfo::new(marker0, marker1, sample_kmers);
/// 
/// ```
/// 
#[derive(Clone, Serialize, Deserialize, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct SampleInfo {
    marker0: Marker,
    marker1: Marker,
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
    /// * `sample_kmers`: a [`Vec<u64>`] containing numbers of non-unique k-mers for each sample
    ///   at the index of the sample-id
    pub fn new(marker0: Marker, marker1: Marker, sample_kmers: Vec<u64>) -> Self {
        let count0 = marker0.count_ones() as u8;
        let count1 = marker1.count_ones() as u8;
        assert_eq!(count0+count1, sample_kmers.len() as u8);

        SampleInfo { marker0, marker1, count0, count1, sample_kmers }
    }

    /// make a new, empty [`SampleInfo`]
    pub fn empty() -> Self {
        SampleInfo { marker0: 0, marker1: 0, count0: 0, count1: 0, sample_kmers: Vec::new() }
    }

    /// get the binary encoded group affiliation of tags
    pub fn get_markers(&self) -> (Marker, Marker) {
        (self.marker0, self.marker1)
    }
}

/// count the ocurrences of the tags
fn tag_counter(tag_vec: &[Tag]) -> Vec<u32> {
    let mut tag_counter = 1;
    let mut tag_counts: Vec<u32> = Vec::new();

    // count the occurences of the labels
    for i in 1..tag_vec.len() {
        if tag_vec[i] == tag_vec[i-1] {
            tag_counter += 1;
        } else {
            tag_counts.push(tag_counter);
            tag_counter = 1;
        }
    }
    tag_counts.push(tag_counter);
    tag_counts.shrink_to_fit();

    tag_counts
}

#[derive(Debug)]
struct TagSummary {
    all_exts: Exts,
    tag_vec: Vec<Tag>,
    tag_counts: Vec<u32>,
    sum: u32,
    edge_mults: EdgeMult,
    highest_quality: Option<BaseQuality>    
}

/// summarize the k-mers, exts and labels, also include an [`EdgeMult`]
fn summarize_tags<K: Kmer, F: Iterator<Item = KmerDataItem<K, Tag>>>(items: F) -> TagSummary {
    let mut all_exts = Exts::empty();

    let mut tag_vec: Vec<Tag> = Vec::with_capacity(items.size_hint().0);
    let mut edge_mults = EdgeMult::new();
    let mut highest_quality = None;

    let mut nobs = 0;
    for item in items {
        tag_vec.push(item.data); 
        all_exts = all_exts.add(item.exts);
        edge_mults.add_exts(item.exts);
        nobs += 1;

        if let Some(q) = item.quality {
            if let Some(hq) = highest_quality {
                if q > hq { highest_quality = Some(q) }
            } else {
                highest_quality = Some(q)
            }
        }
    }

    assert_eq!(all_exts, edge_mults.exts());

    tag_vec.sort();

    let tag_counts = tag_counter(&tag_vec);

    tag_vec.dedup();

    TagSummary { all_exts, tag_vec, tag_counts, sum: nobs, edge_mults, highest_quality }
}

fn summarize_tags_edge_q<K: Kmer, F: Iterator<Item = KmerDataItem<K, Tag>>>(items: F, config: &SummaryConfig) 
-> TagSummary
{    
    // filter the k-mer occurences by their quality -> only use exts and data from k-mers with good enough quality
    let items_filtered = items.filter(|item| 
        match item.quality {
            None => true,
            Some(q) => q >= config.min_quality_for_edge
        }
    );
    
    summarize_tags(items_filtered)
}

#[derive(Debug)]
struct IDTagSummary {
    all_exts: Exts,
    tag_vec: Vec<Tag>,
    tag_counts: Vec<u32>,
    sum: u32,
    id_vec: Vec<ID>,
    edge_mults: EdgeMult,
    highest_quality: Option<BaseQuality>
}

/// summarize the k-mers, exts and labels
fn summarize_tags_ids<K: Kmer, F: Iterator<Item = KmerDataItem<K, IDTag>>>(items: F) -> IDTagSummary {
    let mut all_exts = Exts::empty();

    let mut tag_vec = Vec::with_capacity(items.size_hint().0);
    let mut id_vec = Vec::new();
    let mut edge_mults = EdgeMult::new();
    let mut highest_quality = None;

    let mut nobs = 0;
    for item in items {
        tag_vec.push(item.data.tag); 
        id_vec.push(item.data.id);
        all_exts = all_exts.add(item.exts);
        edge_mults.add_exts(item.exts);
        nobs += 1;

        if let Some(q) = item.quality {
            if let Some(hq) = highest_quality {
                if q > hq { highest_quality = Some(q) }
            } else {
                highest_quality = Some(q)
            }
        }
    }

    tag_vec.sort();
    id_vec.sort();

    let tag_counts = tag_counter(&tag_vec);

    tag_vec.dedup();
    id_vec.dedup();
    id_vec.shrink_to_fit();

    IDTagSummary {all_exts, tag_vec, tag_counts, sum: nobs, id_vec, edge_mults, highest_quality}
}

fn summarize_tags_ids_edge_q<K: Kmer, F: Iterator<Item = KmerDataItem<K, IDTag>>>(items: F, config: &SummaryConfig) 
-> IDTagSummary
{
    // filter the k-mer occurences by their quality -> only use exts and data from k-mers with good enough quality
    let items_filtered = items.filter(|item|
        match item.quality {
            None => true,
            Some(q) => q >= config.min_quality_for_edge
        }
    );

    summarize_tags_ids(items_filtered)
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
fn valid_counts(tags: Tags, nobs: Option<u32>, config: &SummaryConfig) -> bool {
    let nobs_valid = match nobs {
        Some(n) => n as usize >= config.min_kmer_obs,
        None => true
    };

    match config.group_frac {
        GroupFrac::None => nobs_valid,
        GroupFrac::Both => {
            // get amount of labels in tags from each group
            let dist0= tags.bit_and_dist(config.sample_info.marker0);
            let dist1= tags.bit_and_dist(config.sample_info.marker1);
    
            assert_eq!(dist0 + dist1, tags.to_tag_vec().len());
    
            // valid if:
            // - n obs >= min obs AND
            // - observed in at least one third of samples in both groups
            nobs_valid
                && (dist0 as f32 / config.sample_info.count0 as f32 >= config.frac_cutoff) 
                && (dist1 as f32 / config.sample_info.count1 as f32 >= config.frac_cutoff)
        },
        GroupFrac::One => {
            // get amount of labels in tags from each group
            let dist0= tags.bit_and_dist(config.sample_info.marker0);
            let dist1= tags.bit_and_dist(config.sample_info.marker1);

    
            assert_eq!(dist0 + dist1, tags.to_tag_vec().len());
    
            // valid if:
            // - n obs >= min obs AND
            // - observed in at least one third of samples in one group
            nobs_valid
                && ((dist0 as f32 / config.sample_info.count0 as f32 >= config.frac_cutoff) 
                    | (dist1 as f32 / config.sample_info.count1 as f32 >= config.frac_cutoff))
        }
    }
}

enum PInfo<'a> {
    PValue { p: f32 },
    Calculate { tag_vec: &'a [Tag], tag_counts: &'a Vec<u32> }
}

fn valid_p(p_info: PInfo, config: &SummaryConfig) -> bool {
    match config.max_p {
        Some(max_p) => {
            match p_info {
                PInfo::PValue { p } => p <= max_p,
                PInfo::Calculate { tag_vec, tag_counts } => {
                    match p_value(tag_vec, tag_counts, config) {
                        Ok(p) => p <= max_p,
                        Err(_) => true
                    } 
                }
            }
        },
        None => true
    }
}

fn p_value(tag_vec: &[Tag], tag_counts: &[u32], config: &SummaryConfig) -> Result<f32, NotEnoughSamplesError> {
    match config.stat_test {
        StatTest::StudentsTTest => students_t_test(tag_vec, tag_counts, &config.sample_info),
        StatTest::WelchsTTest => welchs_t_test(tag_vec, tag_counts, &config.sample_info),
        StatTest::UTest => u_test(tag_vec, tag_counts, &config.sample_info),
    }
}

// perform a student's t-test
fn students_t_test(tag_vec: &[Tag], tag_counts: &[u32], sample_info: &SampleInfo) -> Result<f32, NotEnoughSamplesError> {
    let n0 = sample_info.count0 as f64;
    let n1 = sample_info.count1 as f64;

    if (n0 < 2.) | (n1 < 2.) { return Err(NotEnoughSamplesError {})}

    let mut counts_g0 = Vec::new();
    let mut counts_g1 = Vec::new();

    let (m0, m1) = sample_info.get_markers();

    for (label, count) in tag_vec.iter().zip(tag_counts) {
        let bin_rep = (2 as Marker).pow(*label as u32);
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
fn welchs_t_test(tag_vec: &[Tag], tag_counts: &[u32], sample_info: &SampleInfo) -> Result<f32, NotEnoughSamplesError> {
    let n0 = sample_info.count0 as f64;
    let n1 = sample_info.count1 as f64;

    if (n0 < 2.) | (n1 < 2.) { return Err(NotEnoughSamplesError {})}

    let mut counts_g0 = Vec::new();
    let mut counts_g1 = Vec::new();

    let (m0, m1) = sample_info.get_markers();

    for (label, count) in tag_vec.iter().zip(tag_counts) {
        let bin_rep = (2 as Marker).pow(*label as u32);
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
fn u_test(tag_vec: &[Tag], tag_counts: &[u32], sample_info: &SampleInfo) -> Result<f32, NotEnoughSamplesError> {

    let n0 = sample_info.count0 as f64;
    let n1 = sample_info.count1 as f64;

    if (n0 < 2.) | (n1 < 2.) { return Err(NotEnoughSamplesError {})}


    let mut counts_g0 = Vec::new();
    let mut counts_g1 = Vec::new();

    let (m0, m1) = sample_info.get_markers();

    for (label, count) in tag_vec.iter().zip(tag_counts) {
        let bin_rep = (2 as Marker).pow(*label as u32);
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
fn log2_fold_change(tags: Tags, counts: &[u32], sample_info: &SampleInfo) -> f32 {
    let mut norm_count_g0 = 0.;
    let mut norm_count_g1 = 0.;

    let (m0, m1) = sample_info.get_markers();

    for (label, count) in tags.to_tag_vec().iter().zip(counts) {
        let bin_rep = (2 as Marker).pow(*label as u32);
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
pub trait SummaryData<DI>: Clone + Debug + Send + Sync + PartialEq + Serialize + DeserializeOwned {
    /// format the noda data 
    fn print(&self, translator: &Translator, config: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String;
    /// format the noda data in for json
    fn print_ol(&self, translator: &Translator, config: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String;
    /// format the noda data in one line
    fn print_json(&self, translator: &Translator, config: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String;
    /// get `Tags` and the overall count, returns `None` if data is insufficient
    fn tags(&self) -> Option<Tags>;
    /// get the size of the structure, including contents of boxed slices
    fn mem(&self) -> usize;
    /// get the number of observations, returns `None` if data is insufficient
    fn sum(&self) -> Option<usize>;
    /// get the IDs, returns `None` if data is insufficient
    fn ids(&self) -> Option<&[ID]>;
    /// get the p-value, returns `None` if data is insufficient
    fn p_value(&self, config: &SummaryConfig) -> Option<f32>;
    /// get the log2(fold change), returns `None` if data is insufficient
    fn fold_change(&self, config: &SummaryConfig) -> Option<f32>;
    /// get the number of samples the sequence was observed in, returns `None` if data is insufficient
    fn sample_count(&self) -> Option<usize>;
    /// get the coverage of the node edges
    fn edge_mults(&self) -> Option<&EdgeMult>;
    /// get the quality of the node k-mer
    fn quality(&self) -> Option<BaseQuality>;
    /// fix the [`EdgeMult`] by removing hanging edges
    fn fix_edge_mults(&mut self, exts: Exts);
    /// set the edge mults
    fn set_edge_mults(&mut self, edge_mults: Option<EdgeMult>);
    /// get a reference to the mapped ids,  returns `None` if data is insufficient
    fn mapped_ids(&self) -> Option<&[ID]>;
    /// add mapped ids to the data
    fn set_mapped_ids(&mut self, mapped_ids: Box<[ID]>);
    /// check if the data can be joined into one
    fn join_test(&self, other: &Self) -> bool;
    /// check if node is valid according to: min kmer obs, group fraction, p-value
    fn valid(&self, config: &SummaryConfig) -> bool;
    /// summarize k-mers
    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, DI>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self);
    /// check summerizer kind
    fn summarizer() -> Summarizers;
}
// TODO: move SummaryData::print functionality to Display trait?

/// Number of observations for the k-mer
impl SummaryData<Tag> for u32 {
    fn print(&self, _: &Translator, _: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        format!("sum: {}", self)
    }

    fn print_ol(&self, _: &Translator, _: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        format!("sum: {}", self)
    }

    fn print_json(&self, _: &Translator, _: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        format!("\"sum\": {}", self)
    }

    fn tags(&self) -> Option<Tags> { None }

    fn mem(&self) -> usize {
        mem::align_of_val(self)
    }

    fn sum(&self) -> Option<usize> {
        Some(*self as usize)
    }

    fn ids(&self) -> Option<&[ID]> { None }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn sample_count(&self) -> Option<usize> { None }

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn quality(&self) -> Option<BaseQuality> { None }

    fn fix_edge_mults(&mut self, _: Exts) { }

    fn set_edge_mults(&mut self, _: Option<EdgeMult>) { }

    fn mapped_ids(&self) -> Option<&[ID]> { None }

    fn set_mapped_ids(&mut self, _: Box<[ID]>) { }

    fn join_test(&self, other: &Self) -> bool {
        self == other
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        *self >= config.min_kmer_obs as u32
    }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, Tag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let summary = summarize_tags_edge_q(items, config);

        let valid_p = valid_p(PInfo::Calculate { tag_vec: &summary.tag_vec, tag_counts: &summary.tag_counts}, config);
        let valid_q = if let Some(q) = summary.highest_quality { q >= config.min_quality } else {true };

        let tags = Tags::from_tag_vec(summary.tag_vec);

        let valid  = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;

        let sum = match config.significant {
            Some(digits) => round_digits(summary.sum, digits),
            None => summary.sum  
        };

        (valid, summary.all_exts, sum)
    }

    fn summarizer() -> Summarizers {
        Summarizers::Sum
    }
}

/// data the k-mer was observed with
impl SummaryData<Tag> for Vec<Tag> {
    fn print(&self, translator: &Translator, _: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        if let Some(tag_translator) = translator.tag_translator() {
            let samples = self
                .iter()
                .map(|sample_id| tag_translator.get_by_right(sample_id).expect("Error: sample does not exist"))
                .collect::<Vec<_>>();
            format!("samples: {:?}", samples).replace("\"", "\'")
        } else {
            format!("samples: {:?}", self).replace("\"", "\'")
        }         
    }

    fn print_ol(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        // print is only one line anyways
        self.print(translator, config, None)
    }

    fn print_json(&self, translator: &Translator, _: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        if let Some(tag_translator) = translator.tag_translator() {
            let samples = self
                .iter()
                .map(|sample_id| tag_translator.get_by_right(sample_id).expect("Error: sample does not exist"))
                .collect::<Vec<_>>();
            format!("\"samples\": {:?}", samples)
        } else {
            format!("\"samples\": {:?}", self)
        } 
    }

    fn tags(&self) -> Option<Tags> { 
        Some(Tags::from_tag_vec(self.clone()))
    }

    fn mem(&self) -> usize {
        mem::size_of_val(&**self) + mem::size_of_val(self)
    }

    fn sum(&self) -> Option<usize> { None }

    fn ids(&self) -> Option<&[ID]> { None }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn sample_count(&self) -> Option<usize> {
        Some(self.len())
    }

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn quality(&self) -> Option<BaseQuality> { None }
    
    fn fix_edge_mults(&mut self, _: Exts) { }

    fn set_edge_mults(&mut self, _: Option<EdgeMult>) { }

        fn mapped_ids(&self) -> Option<&[ID]> { None }

    fn set_mapped_ids(&mut self, _: Box<[ID]>) { }

    fn join_test(&self, other: &Self) -> bool {
        self == other
    }

    fn valid(&self, _: &SummaryConfig) -> bool {
        true
    }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, Tag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let summary = summarize_tags_edge_q(items, config);

        let valid_p = valid_p(PInfo::Calculate { tag_vec: &summary.tag_vec, tag_counts: &summary.tag_counts}, config);
        let valid_q = if let Some(q) = summary.highest_quality { q >= config.min_quality } else {true };

        let tags = Tags::from_tag_vec(summary.tag_vec.clone());

        let valid  = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;
        
        (valid, summary.all_exts, summary.tag_vec)
    }

    fn summarizer() -> Summarizers {
        Summarizers::VecTags
    }
} 

/// the IDs the k-mer was observed with and its number of observations
/// ID could be gene-, read-, or orthogroup-ID
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
// aligned would be 16 Bytes, packed would be 12 Bytes
pub struct IDData {
    ids: Box<[ID]>,
}

impl SummaryData<IDTag> for IDData {
    fn print(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {       
        format!("IDs: {}", id_format(&self.ids, translator, id_group_translator)).replace("\"", "\'") // replace " with ' to avoid conflicts in dot file
    }

    fn print_ol(&self, translator: &Translator, config: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        // print is only one line anyways
        self.print(translator, config, id_group_translator)
    }

    fn print_json(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        format!("\"ids\": {}", id_format(&self.ids, translator, id_group_translator)) // rempve " to avoid conflicts in json file
    }

    fn tags(&self) -> Option<Tags> { None }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.ids)
    }

    fn sum(&self) -> Option<usize> { None }

    fn ids(&self) -> Option<&[ID]> {
        Some(&self.ids[..])
    }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn sample_count(&self) -> Option<usize> { None }

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn quality(&self) -> Option<BaseQuality> { None }

    fn fix_edge_mults(&mut self, _: Exts) { }
    
    fn set_edge_mults(&mut self, _: Option<EdgeMult>) { }

    fn mapped_ids(&self) -> Option<&[ID]> { None }

    fn set_mapped_ids(&mut self, _: Box<[ID]>) { }

    fn join_test(&self, other: &Self) -> bool {
        self == other
    }

    fn valid(&self, _: &SummaryConfig) -> bool { true }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, IDTag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
                let summary = summarize_tags_ids_edge_q(items, config);

        // caluclate p-value with chosen test
        let valid_p = valid_p(PInfo::Calculate { tag_vec: &summary.tag_vec, tag_counts: &summary.tag_counts}, config);
        let valid_q = if let Some(q) = summary.highest_quality { q >= config.min_quality  } else { true };            

        let ids = summary.id_vec.into();
        let tags = Tags::from_tag_vec(summary.tag_vec);

        let valid = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;

        (valid, summary.all_exts, IDData { ids }) 
    }

    fn summarizer() -> Summarizers {
        Summarizers::ID
    }
}

/// the IDs the k-mer was observed with and its number of observations
/// ID could be gene-, read-, or orthogroup-ID
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
// aligned would be 16 Bytes, packed would be 12 Bytes
pub struct IDSumData {
    ids: Box<[ID]>,
    sum: u32,
}

impl SummaryData<IDTag> for IDSumData {
    fn print(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {       
        format!("IDs: {}, sum: {}", id_format(&self.ids, translator, id_group_translator), self.sum).replace("\"", "\'") // replace " with ' to avoid conflicts in dot file
    }

    fn print_ol(&self, translator: &Translator, config: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        // print is only one line anyways
        self.print(translator, config, id_group_translator)
    }

    fn print_json(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        format!("\"ids\": {}, \"sum\": {}", id_format(&self.ids, translator, id_group_translator), self.sum) // rempve " to avoid conflicts in json file
    }

    fn tags(&self) -> Option<Tags> { None }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.ids)
    }

    fn sum(&self) -> Option<usize> {
        Some(self.sum as usize)
    }

    fn ids(&self) -> Option<&[ID]> {
        Some(&self.ids[..])
    }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn sample_count(&self) -> Option<usize> { None }

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn quality(&self) -> Option<BaseQuality> { None }

    fn fix_edge_mults(&mut self, _: Exts) { }
    
    fn set_edge_mults(&mut self, _: Option<EdgeMult>) { }

    fn mapped_ids(&self) -> Option<&[ID]> { None }

    fn set_mapped_ids(&mut self, _: Box<[ID]>) { }

    fn join_test(&self, other: &Self) -> bool {
        self == other
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        self.sum >= config.min_kmer_obs as u32
    }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, IDTag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let summary = summarize_tags_ids_edge_q(items, config);

        // caluclate p-value with chosen test
        let valid_p = valid_p(PInfo::Calculate { tag_vec: &summary.tag_vec, tag_counts: &summary.tag_counts}, config);
        let valid_q = if let Some(q) = summary.highest_quality { q >= config.min_quality  } else { true };            

        let ids = summary.id_vec.into();
        let tags = Tags::from_tag_vec(summary.tag_vec);

        let valid = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;

        let sum = match config.significant {
            Some(digits) => round_digits(summary.sum, digits),
            None => summary.sum  
        };

        (valid, summary.all_exts, IDSumData { ids, sum }) 
    }

    fn summarizer() -> Summarizers {
        Summarizers::IDSum
    }
}

/// the tags the k-mer was observed with
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TagsData {
    tags: Tags,
}

impl SummaryData<Tag> for TagsData {
    fn print(&self, translator: &Translator, _: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        // replace " with ' to avoid conflicts in dot file
        format!("{}", TagsFormatter::new(self.tags, translator)).replace("\"", "\'")
    }

    fn print_ol(&self, translator: &Translator, _: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        if let Some(tag_translator) = translator.tag_translator() {
            format!("samples: {:?}", self.tags.to_string_vec(tag_translator))
        } else {
            format!("samples: {:?}", self.tags.to_tag_vec())
        }.replace("\"", "\'") // replace " with ' to avoid conflicts in dot file
    }

    fn print_json(&self, translator: &Translator, _: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let labels = if let Some(tag_translator) = translator.tag_translator() {
            format!("{:?}", self.tags.to_string_vec(tag_translator))
        } else {
            format!("{:?}", self.tags.to_tag_vec())
        }; // rempve " to avoid conflicts in json file

        format!("\"samples\": {labels}")
    }

    fn tags(&self) -> Option<Tags> { 
        Some(self.tags)
    }

    fn mem(&self) -> usize {
        mem::size_of_val(self)
    }

    fn sum(&self) -> Option<usize> { None }

    fn ids(&self) -> Option<&[ID]> { None }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn sample_count(&self) -> Option<usize> {
        Some(self.tags.len())
    }

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn quality(&self) -> Option<BaseQuality> { None }

    fn fix_edge_mults(&mut self, _: Exts) { }
    
    fn set_edge_mults(&mut self, _: Option<EdgeMult>) { }

    fn mapped_ids(&self) -> Option<&[ID]> { None }

    fn set_mapped_ids(&mut self, _: Box<[ID]>) { }

    fn join_test(&self, other: &Self) -> bool {
        self == other
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        valid_counts(self.tags, None, config)
    }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, Tag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let summary = summarize_tags_edge_q(items, config);

        let valid_p = valid_p(PInfo::Calculate { tag_vec: &summary.tag_vec, tag_counts: &summary.tag_counts}, config);
        let valid_q = if let Some(q) = summary.highest_quality { q >= config.min_quality } else {true };

        let tags = Tags::from_tag_vec(summary.tag_vec);

        let valid  = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;
        
        (valid, summary.all_exts, TagsData { tags })
    }

    fn summarizer() -> Summarizers {
        Summarizers::Tags
    }
}

/// the tags the k-mer was observed with and its number of observations
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
// aligned would be 16 Bytes, packed would be 12 Bytes
pub struct TagsSumData {
    tags: Tags,
    sum: u32,
}

impl SummaryData<Tag> for TagsSumData {
    fn print(&self, translator: &Translator, _: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        // replace " with ' to avoid conflicts in dot file
        format!("{}sum: {}", TagsFormatter::new(self.tags, translator), self.sum).replace("\"", "\'")
    }

    fn print_ol(&self, translator: &Translator, _: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        if let Some(tag_translator) = translator.tag_translator() {
            format!("samples: {:?}, sum: {}", self.tags.to_string_vec(tag_translator), self.sum)
        } else {
            format!("samples: {:?}, sum: {}", self.tags.to_tag_vec(), self.sum)
        }.replace("\"", "\'") // replace " with ' to avoid conflicts in dot file
    }

    fn print_json(&self, translator: &Translator, _: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let labels = if let Some(tag_translator) = translator.tag_translator() {
            format!("{:?}", self.tags.to_string_vec(tag_translator))
        } else {
            format!("{:?}", self.tags.to_tag_vec())
        }; // rempve " to avoid conflicts in json file

        format!("\"samples\": {labels}, \"sum\": {}", self.sum)
    }

    fn tags(&self) -> Option<Tags> {
        Some(self.tags)
    }

    fn mem(&self) -> usize {
        mem::size_of_val(self)
    }

    fn sum(&self) -> Option<usize> {
        Some(self.sum as usize)
    }

    fn ids(&self) -> Option<&[ID]> { None }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn sample_count(&self) -> Option<usize> {
        Some(self.tags.len())
    }

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn quality(&self) -> Option<BaseQuality> { None }

    fn fix_edge_mults(&mut self, _: Exts) { }
    
    fn set_edge_mults(&mut self, _: Option<EdgeMult>) { }

    fn mapped_ids(&self) -> Option<&[ID]> { None }

    fn set_mapped_ids(&mut self, _: Box<[ID]>) { }

    fn join_test(&self, other: &Self) -> bool {
        self == other
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        valid_counts(self.tags, Some(self.sum), config)
    }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, Tag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let summary = summarize_tags_edge_q(items, config);

        let valid_p = valid_p(PInfo::Calculate { tag_vec: &summary.tag_vec, tag_counts: &summary.tag_counts}, config);
        let valid_q = if let Some(q) = summary.highest_quality { q >= config.min_quality } else {true };

        let tags = Tags::from_tag_vec(summary.tag_vec);

        let valid  = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;

        let sum = match config.significant {
            Some(digits) => round_digits(summary.sum, digits),
            None => summary.sum  
        };
        
        (valid, summary.all_exts, TagsSumData { tags, sum })
    }

    fn summarizer() -> Summarizers {
        Summarizers::TagsSum
    }
}

/// Implementation of [`SummaryData<Tag>`]
/// 
/// Contains the tags the k-mer was observed with, how many times it 
/// was observed with each label, and how many times it was observed overall
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TagsCountsSumData {
    tags: Tags,
    sum: u32,
    counts: Box<[u32]>,
}

impl SummaryData<Tag> for TagsCountsSumData {
    fn print(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        format!("{}sum: {}{}{}", TagsCountsFormatter::new(self.tags, &self.counts, translator), self.sum, p, fc).replace("\"", "\'")
    }

    fn print_ol(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        if let Some(tag_translator) = translator.tag_translator() {
            format!("samples: {:?}, counts: {:?}, sum: {}{}{}", self.tags.to_string_vec(tag_translator), self.counts, self.sum, p, fc)
        } else {
            format!("samples: {:?}, counts: {:?}, sum: {}{}{}", self.tags.to_tag_vec(), self.counts, self.sum, p, fc)
        }.replace("\"", "\'")
    }

    fn print_json(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", \"p_value\": {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", \"fold_change\": {}", fc),
            None => "".to_string()
        };

        let labels = if let Some(tag_translator) = translator.tag_translator() {
            format!("{:?}", self.tags.to_string_vec(tag_translator))
        } else {
            format!("{:?}", self.tags.to_tag_vec())
        };

        format!("\"sum\": {}, \"samples\": {labels}, \"counts\": {:?}{p}{fc}", self.sum, self.counts)
    }

    fn tags(&self) -> Option<Tags> {
        Some(self.tags)
    }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.counts)
    }

    fn sum(&self) -> Option<usize> {
        Some(self.sum as usize)
    }

    fn ids(&self) -> Option<&[ID]> { None }

    fn sample_count(&self) -> Option<usize> {
        Some(self.counts.len())
    }

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn quality(&self) -> Option<BaseQuality> { None }

    fn fix_edge_mults(&mut self, _: Exts) { }
    
    fn set_edge_mults(&mut self, _: Option<EdgeMult>) { }

    fn mapped_ids(&self) -> Option<&[ID]> { None }

    fn set_mapped_ids(&mut self, _: Box<[ID]>) { }

    fn join_test(&self, other: &Self) -> bool {
        self == other
    }

    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {      
        p_value(&self.tags.to_tag_vec(), &self.counts, config).ok()
    }


    fn fold_change(&self, config: &SummaryConfig) -> Option<f32> {
        Some(log2_fold_change(self.tags, &self.counts, &config.sample_info))
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        let valid_p = match config.max_p {
            Some(p) => self.p_value(config).expect("error calculating p_value") <= p,
            None => true,
        };

        valid_counts(self.tags, Some(self.sum), config) && valid_p
    }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, Tag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let summary = summarize_tags_edge_q(items, config);

        let valid_p = valid_p(PInfo::Calculate { tag_vec: &summary.tag_vec, tag_counts: &summary.tag_counts}, config);
        let valid_q = if let Some(q) = summary.highest_quality { q >= config.min_quality } else {true };

        let counts: Box<[u32]> = summary.tag_counts.into();
        let tags = Tags::from_tag_vec(summary.tag_vec);

        let valid  = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;

        let sum = match config.significant {
            Some(digits) => round_digits(summary.sum, digits),
            None => summary.sum  
        };

        (valid, summary.all_exts, TagsCountsSumData { tags, counts, sum })
    }

    fn summarizer() -> Summarizers {
        Summarizers::TagsCountsSum
    }
}

/// Implementation of [`SummaryData<Tag>`]
/// 
/// Contains the tags the k-mer was observed with and how many times it 
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

impl SummaryData<Tag> for TagsCountsData {
    fn print(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        format!("{}sum: {}{}{}", TagsCountsFormatter::new(self.tags, &self.counts, translator), self.sum(), p, fc).replace("\"", "\'")
    }

    fn print_ol(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        if let Some(tag_translator) = translator.tag_translator() {
            format!("samples: {:?}, counts: {:?}, sum: {}{}{}", self.tags.to_string_vec(tag_translator), self.counts, self.sum(), p, fc)
        } else {
            format!("samples: {:?}, counts: {:?}, sum: {}{}{}", self.tags.to_tag_vec(), self.counts, self.sum(), p, fc)
        }.replace("\"", "\'")
    }

    fn print_json(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", \"p_value\": {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", \"fold_change\": {}", fc),
            None => "".to_string()
        };

        let labels = if let Some(tag_translator) = translator.tag_translator() {
            format!("{:?}", self.tags.to_string_vec(tag_translator))
        } else {
            format!("{:?}", self.tags.to_tag_vec())
        };

        format!("\"sum\": {}, \"samples\": {labels}, \"counts\": {:?}{p}{fc}", self.sum(), self.counts)
    }

    fn tags(&self) -> Option<Tags> {
        Some(self.tags)
    }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.counts)
    }

    fn sum(&self) -> Option<usize> {
        Some(self.counts.iter().sum::<u32>() as usize)
    }

    fn ids(&self) -> Option<&[ID]> { None }

    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {      
        p_value(&self.tags.to_tag_vec(), &self.counts, config).ok()
    }


    fn fold_change(&self, config: &SummaryConfig) -> Option<f32> {
        Some(log2_fold_change(self.tags, &self.counts, &config.sample_info))
    }

    fn sample_count(&self) -> Option<usize> {
        Some(self.counts.len())
    }

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn quality(&self) -> Option<BaseQuality> { None }

    fn fix_edge_mults(&mut self, _: Exts) { }
    
    fn set_edge_mults(&mut self, _: Option<EdgeMult>) { }

    fn mapped_ids(&self) -> Option<&[ID]> { None }

    fn set_mapped_ids(&mut self, _: Box<[ID]>) { }

    fn join_test(&self, other: &Self) -> bool {
        self == other
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        let valid_p = match config.max_p {
            Some(p) => self.p_value(config).expect("error calculating p-value") <= p,
            None => true,
        }; 

        valid_counts(self.tags, Some(self.sum()), config) && valid_p
    }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, Tag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let summary = summarize_tags_edge_q(items, config);

        let valid_p = valid_p(PInfo::Calculate { tag_vec: &summary.tag_vec, tag_counts: &summary.tag_counts}, config);
        let valid_q = if let Some(q) = summary.highest_quality { q >= config.min_quality } else {true };

        let counts: Box<[u32]> = summary.tag_counts.into();
        let tags = Tags::from_tag_vec(summary.tag_vec);

        let valid  = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;

        (valid && valid_p, summary.all_exts, TagsCountsData { tags, counts }) 
    }

    fn summarizer() -> Summarizers {
        Summarizers::TagsCounts
    }
}

/// Implementation of [`SummaryData<Tag>`]
/// 
/// Contains the tags the k-mer was observed with, how many times it 
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

impl SummaryData<Tag> for TagsCountsPData {
    fn print(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        format!("{}sum: {}{}{}", TagsCountsFormatter::new(self.tags, &self.counts, translator), self.sum(), p, fc).replace("\"", "\'")
    }

    fn print_ol(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        if let Some(tag_translator) = translator.tag_translator() {
            format!("samples: {:?}, counts: {:?}, sum: {}{}{}", self.tags.to_string_vec(tag_translator), self.counts, self.sum(), p, fc)
        } else {
            format!("samples: {:?}, counts: {:?}, sum: {}{}{}", self.tags.to_tag_vec(), self.counts, self.sum(), p, fc)
        }.replace("\"", "\'")
    }

    fn print_json(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", \"p_value\": {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", \"fold_change\": {}", fc),
            None => "".to_string()
        };

        let labels = if let Some(tag_translator) = translator.tag_translator() {
            format!("{:?}", self.tags.to_string_vec(tag_translator))
        } else {
            format!("{:?}", self.tags.to_tag_vec())
        };

        format!("\"sum\": {}, \"samples\": {labels}, \"counts\": {:?}{p}{fc}", self.sum(), self.counts)
    }

    fn tags(&self) -> Option<Tags> {
        Some(self.tags)
    }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.counts)
    }

    fn sum(&self) -> Option<usize> {
        Some(self.sum() as usize)
    }

    fn ids(&self) -> Option<&[ID]> { None }

    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {
        if config.stat_test_changed {
            Some(p_value(&self.tags.to_tag_vec(), &self.counts, config).unwrap())
        } else {
            Some(self.p_value)
        } 
    }

    fn fold_change(&self, config: &SummaryConfig) -> Option<f32> {
        Some(log2_fold_change(self.tags, &self.counts, &config.sample_info))
    }
    
    fn sample_count(&self) -> Option<usize> {
        Some(self.counts.len())
    }

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn quality(&self) -> Option<BaseQuality> { None }

    fn fix_edge_mults(&mut self, _: Exts) { }
    
    fn set_edge_mults(&mut self, _: Option<EdgeMult>) { }

    fn mapped_ids(&self) -> Option<&[ID]> { None }

    fn set_mapped_ids(&mut self, _: Box<[ID]>) { }

    fn join_test(&self, other: &Self) -> bool {
        self == other
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        valid_counts(self.tags, Some(self.sum()), config) 
            && valid_p(PInfo::PValue { p: self.p_value(config).expect("error getting p-values") }, config)
    }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, Tag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let summary = summarize_tags_edge_q(items, config);

        let p_value = p_value(&summary.tag_vec, &summary.tag_counts, config).unwrap();

        let valid_p = valid_p(PInfo::PValue { p: p_value }, config);
        let valid_q = if let Some(q) = summary.highest_quality { q >= config.min_quality } else {true };

        let counts: Box<[u32]> = summary.tag_counts.into();
        let tags = Tags::from_tag_vec(summary.tag_vec);

        let valid  = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;

        (valid, summary.all_exts, TagsCountsPData { tags, counts, p_value }) 
    }

    fn summarizer() -> Summarizers {
        Summarizers::TagsCountsP
    }
}

/// Implementation of [`SummaryData<Tag>`]
/// 
/// Contains the tags the k-mer was observed with, how many times it 
/// was observed with each label, and the edge multiplicites/coverage
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

impl SummaryData<Tag> for TagsCountsEMData {
    fn print(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        format!("{}sum: {}{}{}, edge coverage: \n{}", TagsCountsFormatter::new(self.tags, &self.counts, translator), self.sum(), p, fc, self.edge_mults).replace("\"", "\'")
    }

    fn print_ol(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        if let Some(tag_translator) = translator.tag_translator() {
            format!("samples: {:?}, counts: {:?}, sum: {}{}{}, edge coverage: {:?}", self.tags.to_string_vec(tag_translator), self.counts, self.sum(), p, fc, self.edge_mults)
        } else {
            format!("samples: {:?}, counts: {:?}, sum: {}{}{}, edge coverage: {:?}", self.tags.to_tag_vec(), self.counts, self.sum(), p, fc, self.edge_mults)
        }.replace("\"", "\'")
    }

    fn print_json(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", \"p_value\": {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", \"fold_change\": {}", fc),
            None => "".to_string()
        };

        let labels = if let Some(tag_translator) = translator.tag_translator() {
            format!("{:?}", self.tags.to_string_vec(tag_translator))
        } else {
            format!("{:?}", self.tags.to_tag_vec())
        };

        format!("\"sum\": {}, \"samples\": {labels}, \"counts\": {:?}{p}{fc}", self.sum(), self.counts)
    }

    fn tags(&self) -> Option<Tags> {
        Some(self.tags)
    }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.counts)
    }

    fn sum(&self) -> Option<usize> {
        Some(self.counts.iter().sum::<u32>() as usize)
    }

    fn ids(&self) -> Option<&[ID]> { None }

    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {      
        p_value(&self.tags.to_tag_vec(), &self.counts, config).ok()
    }

    fn fold_change(&self, config: &SummaryConfig) -> Option<f32> {
        Some(log2_fold_change(self.tags, &self.counts, &config.sample_info))
    }

    fn sample_count(&self) -> Option<usize> {
        Some(self.counts.len())
    }

    fn edge_mults(&self) -> Option<&EdgeMult> {
        Some(&self.edge_mults)
    }

    fn quality(&self) -> Option<BaseQuality> { None }

    fn fix_edge_mults(&mut self, exts: Exts) {
        self.edge_mults.clean_edges(exts);
    }

    fn set_edge_mults(&mut self, edge_mults: Option<EdgeMult>) {
        self.edge_mults = edge_mults.expect("Error: no edge mults")
    }

    fn mapped_ids(&self) -> Option<&[ID]> { None }

    fn set_mapped_ids(&mut self, _: Box<[ID]>) { }

    fn join_test(&self, other: &Self) -> bool {
        self.counts == other.counts
            && self.tags == other.tags
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        let valid_p = match config.max_p {
            Some(p) => self.p_value(config).expect("error calculating p-value") <= p,
            None => true,
        }; 

        valid_counts(self.tags, Some(self.sum()), config) && valid_p
    }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, Tag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let summary = summarize_tags_edge_q(items, config);

        let valid_p = valid_p(PInfo::Calculate { tag_vec: &summary.tag_vec, tag_counts: &summary.tag_counts}, config);
        let valid_q = if let Some(q) = summary.highest_quality { q >= config.min_quality } else {true };

        let counts: Box<[u32]> = summary.tag_counts.into();
        let tags = Tags::from_tag_vec(summary.tag_vec);

        let valid  = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;

        (valid, summary.all_exts, TagsCountsEMData { tags, counts, edge_mults: summary.edge_mults }) 
    }

    fn summarizer() -> Summarizers {
        Summarizers::TagsCountsEM
    }
}

/// Implementation of [`SummaryData<Tag>`]
/// 
/// Contains the tags the k-mer was observed with, how many times it 
/// was observed with each label, a p-value, and the edge multiplicites/coverage
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

impl SummaryData<Tag> for TagsCountsPEMData{
    fn print(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        format!("{}sum: {}{}{}, edge coverage: \n{}", TagsCountsFormatter::new(self.tags, &self.counts, translator), self.sum(), p, fc, self.edge_mults).replace("\"", "\'")
    }

    fn print_ol(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        if let Some(tag_translator) = translator.tag_translator() {
            format!("samples: {:?}, counts: {:?}, sum: {}{}{}, edge coverage: {:?}", self.tags.to_string_vec(tag_translator), self.counts, self.sum(), p, fc, self.edge_mults)
        } else {
            format!("samples: {:?}, counts: {:?}, sum: {}{}{}, edge coverage: {:?}", self.tags.to_tag_vec(), self.counts, self.sum(), p, fc, self.edge_mults)
        }.replace("\"", "\'")    
    }

    fn print_json(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", \"p_value\": {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", \"fold_change\": {}", fc),
            None => "".to_string()
        };

        let labels = if let Some(tag_translator) = translator.tag_translator() {
            format!("{:?}", self.tags.to_string_vec(tag_translator))
        } else {
            format!("{:?}", self.tags.to_tag_vec())
        };

        format!("\"sum\": {}, \"samples\": {labels}, \"counts\": {:?}{p}{fc}", self.sum(), self.counts)
    }

    fn tags(&self) -> Option<Tags> {
        Some(self.tags)
    }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.counts)
    }

    fn sum(&self) -> Option<usize> {
        Some(self.counts.iter().sum::<u32>() as usize)
    }

    fn ids(&self) -> Option<&[ID]> { None }

    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {
        if config.stat_test_changed {
            Some(p_value(&self.tags.to_tag_vec(), &self.counts, config).unwrap())
        } else {
            Some(self.p_value)
        } 
    }

    fn fold_change(&self, config: &SummaryConfig) -> Option<f32> {
        Some(log2_fold_change(self.tags, &self.counts, &config.sample_info))
    }

    fn sample_count(&self) -> Option<usize> {
        Some(self.counts.len())
    }

    fn edge_mults(&self) -> Option<&EdgeMult> {
        Some(&self.edge_mults)
    }

    fn quality(&self) -> Option<BaseQuality> { None }

    fn fix_edge_mults(&mut self, exts: Exts) {
        self.edge_mults.clean_edges(exts);
    }

    fn set_edge_mults(&mut self, edge_mults: Option<EdgeMult>) {
        self.edge_mults = edge_mults.expect("Error: no edge mults")
    }

    fn mapped_ids(&self) -> Option<&[ID]> { None }

    fn set_mapped_ids(&mut self, _: Box<[ID]>) { }

    fn join_test(&self, other: &Self) -> bool {
        self.counts == other.counts
            && self.tags == other.tags
            && self.p_value == other.p_value
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        valid_counts(self.tags, Some(self.sum()), config) 
            && valid_p(PInfo::PValue { p: self.p_value(config).expect("error getting p-values") }, config)
    }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, Tag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let summary = summarize_tags_edge_q(items, config);

        // caluclate p-value with chosen test
        let p_value = p_value(&summary.tag_vec, &summary.tag_counts, config).unwrap();

        let valid_p = valid_p(PInfo::PValue { p: p_value }, config);
        let valid_q = if let Some(q) = summary.highest_quality { q >= config.min_quality } else { true };

        let counts: Box<[u32]> = summary.tag_counts.into();
        let tags = Tags::from_tag_vec(summary.tag_vec);

        let valid  = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;  

        (valid, summary.all_exts, TagsCountsPEMData { tags, counts, p_value, edge_mults: summary.edge_mults }) 
    }

    fn summarizer() -> Summarizers {
        Summarizers::TagsCountsPEM
    }
}

// Implementation of [`SummaryData<Tag>`]
/// 
/// Contains the tags the k-mer was observed with, how many times it 
/// was observed with each label, a p-value, and the edge multiplicites/coverage
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TagsCountsPEMQualityData {
    tags: Tags,
    counts: Box<[u32]>,
    p_value: f32,
    edge_mults: EdgeMult,
    quality: BaseQuality,
}

impl TagsCountsPEMQualityData {
    #[inline]
    pub fn sum(&self) -> u32 {
        self.counts.iter().sum::<u32>()
    }
}

impl SummaryData<Tag> for TagsCountsPEMQualityData{
    fn print(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        format!("{}sum: {}{}{}, quality: {}, edge coverage: \n{}", 
            TagsCountsFormatter::new(self.tags, &self.counts, translator), 
            self.sum(), 
            p, 
            fc, 
            self.quality, 
            self.edge_mults
        ).replace("\"", "\'")
    }

    fn print_ol(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        if let Some(tag_translator) = translator.tag_translator() {
            format!("samples: {:?}, counts: {:?}, sum: {}{}{}, quality: {}, edge coverage: {:?}", 
                self.tags.to_string_vec(tag_translator), 
                self.counts, self.sum(), 
                p, 
                fc, 
                self.quality, 
                self.edge_mults
            )
        } else {
            format!("samples: {:?}, counts: {:?}, sum: {}{}{}, quality: {}, edge coverage: {:?}", 
                self.tags.to_tag_vec(), 
                self.counts, 
                self.sum(), 
                p, 
                fc, 
                self.quality, 
                self.edge_mults
            )
        }.replace("\"", "\'")    
    }

    fn print_json(&self, translator: &Translator, config: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", \"p_value\": {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", \"fold_change\": {}", fc),
            None => "".to_string()
        };

        let labels = if let Some(tag_translator) = translator.tag_translator() {
            format!("{:?}", self.tags.to_string_vec(tag_translator))
        } else {
            format!("{:?}", self.tags.to_tag_vec())
        };

        format!("\"sum\": {}, \"samples\": {labels}, \"counts\": {:?}{p}{fc}, \"quality\": {}", self.sum(), self.counts, self.quality as u8)
    }

    fn tags(&self) -> Option<Tags> {
        Some(self.tags)
    }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.counts)
    }

    fn sum(&self) -> Option<usize> {
        Some(self.counts.iter().sum::<u32>() as usize)
    }

    fn ids(&self) -> Option<&[ID]> { None }

    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {
        if config.stat_test_changed {
            Some(p_value(&self.tags.to_tag_vec(), &self.counts, config).unwrap())
        } else {
            Some(self.p_value)
        } 
    }

    fn fold_change(&self, config: &SummaryConfig) -> Option<f32> {
        Some(log2_fold_change(self.tags, &self.counts, &config.sample_info))
    }

    fn sample_count(&self) -> Option<usize> {
        Some(self.counts.len())
    }

    fn edge_mults(&self) -> Option<&EdgeMult> {
        Some(&self.edge_mults)
    }

    fn quality(&self) -> Option<BaseQuality> { 
        Some(self.quality)
    }

    fn fix_edge_mults(&mut self, exts: Exts) {
        self.edge_mults.clean_edges(exts);
    }

    fn set_edge_mults(&mut self, edge_mults: Option<EdgeMult>) {
        self.edge_mults = edge_mults.expect("Error: no edge mults")
    }

    fn mapped_ids(&self) -> Option<&[ID]> { None }

    fn set_mapped_ids(&mut self, _: Box<[ID]>) { }

    fn join_test(&self, other: &Self) -> bool {
        self.counts == other.counts
            && self.tags == other.tags
            && self.p_value == other.p_value
            && self.quality == other.quality
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        valid_counts(self.tags, Some(self.sum()), config) 
            && valid_p(PInfo::PValue { p: self.p_value(config).expect("error getting p-values") }, config)
            && self.quality >= config.min_quality
    }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, Tag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let summary = summarize_tags_edge_q(items, config);

        // caluclate p-value with chosen test
        let p_value = p_value(&summary.tag_vec, &summary.tag_counts, config).unwrap();
        let quality = summary.highest_quality.expect("missing quality score - required for summarizer");

        let valid_p = valid_p(PInfo::PValue { p: p_value }, config);
        let valid_q = quality >= config.min_quality;

        let counts: Box<[u32]> = summary.tag_counts.into();
        let tags = Tags::from_tag_vec(summary.tag_vec);

        let valid  = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;  

        (valid, summary.all_exts, TagsCountsPEMQualityData { tags, counts, p_value, edge_mults: summary.edge_mults, quality }) 
    }

    fn summarizer() -> Summarizers {
        Summarizers::TagsCountsPEMQuality
    }
}


/// Implementation of [`SummaryData<IDTag>`]
/// 
/// Contains the tags the k-mer was observed with and how many times it 
/// was observed with each label
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IDTagsCountsData {
    ids: Box<[ID]>,
    tags: Tags,
    counts: Box<[u32]>
}

impl IDTagsCountsData {
    #[inline]
    pub fn sum(&self) -> u32 {
        self.counts.iter().sum::<u32>()
    }
}

impl SummaryData<IDTag> for IDTagsCountsData {
    fn print(&self, translator: &Translator, config: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        let ids_format = id_format(&self.ids, translator, id_group_translator);
        format!("IDs: {}, {}sum: {}{}{}", ids_format, TagsCountsFormatter::new(self.tags, &self.counts, translator), self.sum(), p, fc).replace("\"", "\'") // replace " with ' to avoid conflicts in dot file
    }

    fn print_ol(&self, translator: &Translator, config: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        let tags_format = if let Some(tag_translator) = translator.tag_translator() {
            format!("{:?}", self.tags.to_string_vec(tag_translator))
        } else {
            format!("{:?}", self.tags.to_tag_vec())
        };

        let ids_format = id_format(&self.ids, translator, id_group_translator);
        format!("IDs: {}, samples: {}, counts: {:?}, sum: {}{}{}", ids_format, tags_format, self.counts, self.sum(), p, fc).replace("\"", "\'")
    }

    fn print_json(&self, translator: &Translator, config: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", \"p_value\": {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", \"fold_change\": {}", fc),
            None => "".to_string()
        };

        let labels = if let Some(tag_translator) = translator.tag_translator() {
            format!("{:?}", self.tags.to_string_vec(tag_translator))
        } else {
            format!("{:?}", self.tags.to_tag_vec())
        };

        let ids = id_format(&self.ids, translator, id_group_translator);

        format!("\"ids\": {ids}, \"sum\": {}, \"samples\": {labels}, \"counts\": {:?}{p}{fc}", self.sum(), self.counts)
    }

    fn tags(&self) -> Option<Tags> {
        Some(self.tags)
    }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.counts) + mem::size_of_val(&*self.ids)
    }

    fn sum(&self) -> Option<usize> {
        Some(self.counts.iter().sum::<u32>() as usize)
    }

    fn ids(&self) -> Option<&[ID]> {
        Some(&self.ids)
    }

    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {      
        p_value(&self.tags.to_tag_vec(), &self.counts, config).ok()
    }


    fn fold_change(&self, config: &SummaryConfig) -> Option<f32> {
        Some(log2_fold_change(self.tags, &self.counts, &config.sample_info))
    }

    fn sample_count(&self) -> Option<usize> {
        Some(self.counts.len())
    }

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn quality(&self) -> Option<BaseQuality> { None }

    fn fix_edge_mults(&mut self, _: Exts) { }
    
    fn set_edge_mults(&mut self, _: Option<EdgeMult>) { }

    fn mapped_ids(&self) -> Option<&[ID]> { None }

    fn set_mapped_ids(&mut self, _: Box<[ID]>) { }

    fn join_test(&self, other: &Self) -> bool {
        self == other
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        let valid_p = match config.max_p {
            Some(p) => self.p_value(config).expect("error calculating p-value") <= p,
            None => true,
        }; 

        valid_counts(self.tags, Some(self.sum()), config) && valid_p
    }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, IDTag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let summary = summarize_tags_ids_edge_q(items, config);
        
        // caluclate p-value with chosen test
        let p_value = p_value(&summary.tag_vec, &summary.tag_counts, config).unwrap();

        let valid_p = valid_p(PInfo::PValue { p: p_value }, config);
        let valid_q = if let Some(q) = summary.highest_quality { q >= config.min_quality } else {true };

        let counts: Box<[u32]> = summary.tag_counts.into();
        let tags = Tags::from_tag_vec(summary.tag_vec);
        let ids: Box<[ID]> = summary.id_vec.into();

        let valid  = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;     

        (valid, summary.all_exts, IDTagsCountsData { tags, counts, ids }) 
    }

    fn summarizer() -> Summarizers {
        Summarizers::IDTagsCounts
    }
}

/// Implementation of [`SummaryData<Tag>`]
/// 
/// Contains the tags the k-mer was observed with, how many times it 
/// was observed with each label, a p-value, and the edge multiplicites/coverage
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IDTagsCountsPEMData {
    tags: Tags,
    counts: Box<[u32]>,
    ids: Box<[ID]>,
    p_value: f32,
    edge_mults: EdgeMult,
}

impl IDTagsCountsPEMData {
    #[inline]
    pub fn sum(&self) -> u32 {
        self.counts.iter().sum::<u32>()
    }
}

impl SummaryData<IDTag> for IDTagsCountsPEMData{
    fn print(&self, translator: &Translator, config: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        let ids_format = id_format(&self.ids, translator, id_group_translator);

        format!("IDs: {}, {}sum: {}{}{}, edge coverage: \n{}", 
            ids_format, 
            TagsCountsFormatter::new(self.tags, &self.counts, translator), 
            self.sum(), 
            p, 
            fc, 
            self.edge_mults
        ).replace("\"", "\'") // replace " with ' to avoid conflicts in dot file
    }

    fn print_ol(&self, translator: &Translator, config: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", p-value: {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", log2(fold change): {}", fc),
            None => "".to_string()
        };

        let tags_format = if let Some(tag_translator) = translator.tag_translator() {
            format!("{:?}", self.tags.to_string_vec(tag_translator))
        } else {
            format!("{:?}", self.tags.to_tag_vec())
        };

        let ids_format = id_format(&self.ids, translator, id_group_translator);

        format!("IDs: {}, samples: {}, counts: {:?}, sum: {}{}{}, edge coverage: {:?}", 
            ids_format, 
            tags_format, 
            self.counts, 
            self.sum(), 
            p, 
            fc, 
            self.edge_mults
        ).replace("\"", "\'")
    }

    fn print_json(&self, translator: &Translator, config: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let p = match self.p_value(config) {
            Some(p) => format!(", \"p_value\": {}", p),
            None => "".to_string()
        };

        let fc = match self.fold_change(config) {
            Some(fc) => format!(", \"fold_change\": {}", fc),
            None => "".to_string()
        };

        let labels = if let Some(tag_translator) = translator.tag_translator() {
            format!("{:?}", self.tags.to_string_vec(tag_translator))
        } else {
            format!("{:?}", self.tags.to_tag_vec())
        };

        let ids = id_format(&self.ids, translator, id_group_translator);

        format!("\"ids\": {ids}, \"sum\": {}, \"samples\": {labels}, \"counts\": {:?}{p}{fc}", self.sum(), self.counts)
    }

    fn tags(&self) -> Option<Tags> {
        Some(self.tags)
    }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.counts) + mem::size_of_val(&*self.ids)
    }

    fn sum(&self) -> Option<usize> {
        Some(self.counts.iter().sum::<u32>() as usize)
    }

    fn ids(&self) -> Option<&[ID]> {
        Some(&self.ids)
    }

    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {
        if config.stat_test_changed {
            Some(p_value(&self.tags.to_tag_vec(), &self.counts, config).unwrap())
        } else {
            Some(self.p_value)
        } 
    }

    fn fold_change(&self, config: &SummaryConfig) -> Option<f32> {
        Some(log2_fold_change(self.tags, &self.counts, &config.sample_info))
    }

    fn sample_count(&self) -> Option<usize> {
        Some(self.counts.len())
    }

    fn edge_mults(&self) -> Option<&EdgeMult> {
        Some(&self.edge_mults)
    }

    fn quality(&self) -> Option<BaseQuality> { None }

    fn fix_edge_mults(&mut self, exts: Exts) {
        self.edge_mults.clean_edges(exts);
    }

    fn set_edge_mults(&mut self, edge_mults: Option<EdgeMult>) {
        self.edge_mults = edge_mults.expect("Error: no edge mults")
    }

    fn mapped_ids(&self) -> Option<&[ID]> { None }

    fn set_mapped_ids(&mut self, _: Box<[ID]>) { }

    fn join_test(&self, other: &Self) -> bool {
        self.counts == other.counts
            && self.tags == other.tags
            && self.p_value == other.p_value
            && self.ids == other.ids
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        valid_counts(self.tags, Some(self.sum()), config) 
            && valid_p(PInfo::PValue { p: self.p_value(config).expect("error getting p-values") }, config)
    }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, IDTag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let summary = summarize_tags_ids_edge_q(items, config);

        // caluclate p-value with chosen test
        let p_value = p_value(&summary.tag_vec, &summary.tag_counts, config).unwrap();

        let valid_p = valid_p(PInfo::PValue { p: p_value }, config);
        let valid_q = if let Some(q) = summary.highest_quality { q >= config.min_quality } else {true };

        let counts: Box<[u32]> = summary.tag_counts.into();
        let tags = Tags::from_tag_vec(summary.tag_vec);
        let ids: Box<[ID]> = summary.id_vec.into();

        let valid  = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;  

        (valid, summary.all_exts, IDTagsCountsPEMData { tags, counts, p_value, ids, edge_mults: summary.edge_mults }) 
    }

    fn summarizer() -> Summarizers {
        Summarizers::IDTagsCountsPEM
    }
}

/// Implementation of [`SummaryData<Tag>`]
/// 
/// Contains the IDs the k-mer was observed with and the edge multiplicites/coverage
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IDEMData {
    ids: Box<[ID]>,
    edge_mults: EdgeMult,
}

impl SummaryData<IDTag> for IDEMData{
    fn print(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let ids_format = id_format(&self.ids, translator, id_group_translator);

        format!("IDs: {}, edge coverage: \n{}", 
            ids_format, 
            self.edge_mults
        ).replace("\"", "\'") // replace " with ' to avoid conflicts in dot file
    }

    fn print_ol(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let ids_format = id_format(&self.ids, translator, id_group_translator);

        format!("IDs: {}, edge coverage: {:?}", 
            ids_format, 
            self.edge_mults
        ).replace("\"", "\'") // replace " with ' to avoid conflicts in dot file

    }

    fn print_json(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        format!("\"ids\": {}", id_format(&self.ids, translator, id_group_translator)) // rempve " to avoid conflicts in json file
    }

    fn tags(&self) -> Option<Tags> { None }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.ids)
    }

    fn sum(&self) -> Option<usize> { None }

    fn ids(&self) -> Option<&[ID]> {
        Some(&self.ids)
    }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn sample_count(&self) -> Option<usize> { None }

    fn edge_mults(&self) -> Option<&EdgeMult> {
        Some(&self.edge_mults)
    }

    fn quality(&self) -> Option<BaseQuality> { None }

    fn fix_edge_mults(&mut self, exts: Exts) {
        self.edge_mults.clean_edges(exts);
    }

    fn set_edge_mults(&mut self, edge_mults: Option<EdgeMult>) {
        self.edge_mults = edge_mults.expect("Error: no edge mults")
    }

    fn mapped_ids(&self) -> Option<&[ID]> { None }

    fn set_mapped_ids(&mut self, _: Box<[ID]>) { }

    fn join_test(&self, other: &Self) -> bool {
        self.ids == other.ids
    }

    fn valid(&self, _: &SummaryConfig) -> bool { true }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, IDTag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let summary = summarize_tags_ids_edge_q(items, config);

        // caluclate p-value with chosen test, valid if not enough samples
        let valid_p = valid_p(PInfo::Calculate { tag_vec: &summary.tag_vec, tag_counts: &summary.tag_counts}, config);
        let valid_q = if let Some(q) = summary.highest_quality { q >= config.min_quality  } else { true };

        let ids = summary.id_vec.into();
        let tags = Tags::from_tag_vec(summary.tag_vec);

        let valid = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;

        (valid, summary.all_exts, IDEMData { ids, edge_mults: summary.edge_mults }) 
    }

    fn summarizer() -> Summarizers {
        Summarizers::IDEM
    }
}

/// Implementation of [`SummaryData<Tag>`]
/// 
/// Contains the IDs the k-mer was observed with, a placeholder for mapped ids, and edge multiplicites/coverage
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IDMapEMData {
    ids: Box<[ID]>,
    map_ids: Box<[ID]>,
    edge_mults: EdgeMult,
}

impl SummaryData<IDTag> for IDMapEMData{
    fn print(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let ids_format = id_format(&self.ids, translator, id_group_translator);
        let map_ids_format = id_format(&self.map_ids, translator, id_group_translator);

        format!("IDs: {}, mapped IDs: {}, edge coverage: {}", 
            ids_format, 
            map_ids_format,
            self.edge_mults
        ).replace("\"", "\'") // replace " with ' to avoid conflicts in dot file
    }

    fn print_ol(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let ids_format = id_format(&self.ids, translator, id_group_translator);
        let map_ids_format = id_format(&self.map_ids, translator, id_group_translator);

        format!("IDs: {}, mapped IDs: {}, edge coverage: {:?}", 
            ids_format, 
            map_ids_format,
            self.edge_mults
        ).replace("\"", "\'") // replace " with ' to avoid conflicts in dot file
    }

    fn print_json(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let ids_format = id_format(&self.ids, translator, id_group_translator);
        let map_ids_format = id_format(&self.map_ids, translator, id_group_translator);

        let has_mapped = !self.map_ids.is_empty() as usize;

        format!("\"ids\": {ids_format}, \"mapped_ids\": {map_ids_format}, \"has_mapped_ids\": {has_mapped}", ) // rempve " to avoid conflicts in json file
    }

    fn tags(&self) -> Option<Tags> { None }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.ids) + mem::size_of_val(&*self.map_ids)
    }

    fn sum(&self) -> Option<usize> { None }

    fn ids(&self) -> Option<&[ID]> {
        Some(&self.ids)
    }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn sample_count(&self) -> Option<usize> { None }

    fn edge_mults(&self) -> Option<&EdgeMult> {
        Some(&self.edge_mults)
    }

    fn quality(&self) -> Option<BaseQuality> { None }

    fn fix_edge_mults(&mut self, exts: Exts) {
        self.edge_mults.clean_edges(exts);
    }

    fn set_edge_mults(&mut self, edge_mults: Option<EdgeMult>) {
        self.edge_mults = edge_mults.expect("Error: no edge mults")
    }

    fn mapped_ids(&self) -> Option<&[ID]> {
        Some(&self.map_ids)
    }

    fn set_mapped_ids(&mut self, mapped_ids: Box<[ID]>) {
        self.map_ids = mapped_ids
    }

    fn join_test(&self, other: &Self) -> bool {
        self.ids == other.ids 
        && self.map_ids == other.map_ids
    }

    fn valid(&self, _: &SummaryConfig) -> bool { true }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, IDTag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let summary = summarize_tags_ids_edge_q(items, config);

        // caluclate p-value with chosen test
        let valid_p = valid_p(PInfo::Calculate { tag_vec: &summary.tag_vec, tag_counts: &summary.tag_counts}, config);
        let valid_q = if let Some(q) = summary.highest_quality { q >= config.min_quality  } else { true };            

        let ids = summary.id_vec.into();
        let tags = Tags::from_tag_vec(summary.tag_vec);

        let valid = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;

        (valid, summary.all_exts, IDMapEMData { ids, map_ids: Vec::new().into(), edge_mults: summary.edge_mults }) 
    }

    fn summarizer() -> Summarizers {
        Summarizers::IDMapEM
    }
}

/// Implementation of [`SummaryData<Tag>`]
/// 
/// Contains the IDs the k-mer was observed with, a placeholder for mapped ids, and edge multiplicites/coverage
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IDMapEMQualityData {
    ids: Box<[ID]>,
    map_ids: Box<[ID]>,
    edge_mults: EdgeMult,
    quality: BaseQuality
}

impl SummaryData<IDTag> for IDMapEMQualityData{
    fn print(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let ids_format = id_format(&self.ids, translator, id_group_translator);
        let map_ids_format = id_format(&self.map_ids, translator, id_group_translator);

        format!("IDs: {}, mapped IDs: {}, quality: {}, edge coverage: {}", 
            ids_format, 
            map_ids_format,
            self.quality,
            self.edge_mults,
        ).replace("\"", "\'") // replace " with ' to avoid conflicts in dot file
    }

    fn print_ol(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let ids_format = id_format(&self.ids, translator, id_group_translator);
        let map_ids_format = id_format(&self.map_ids, translator, id_group_translator);

        format!("IDs: {}, mapped IDs: {}, quality: {}, edge coverage: {:?}", 
            ids_format, 
            map_ids_format,
            self.quality,
            self.edge_mults
        ).replace("\"", "\'") // replace " with ' to avoid conflicts in dot file
    }

    fn print_json(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let ids_format = id_format(&self.ids, translator, id_group_translator);
        let map_ids_format = id_format(&self.map_ids, translator, id_group_translator);

        let has_mapped = !self.map_ids.is_empty() as usize;

        format!("\"ids\": {ids_format}, \"mapped_ids\": {map_ids_format}, \"has_mapped_ids\": {has_mapped}, \"quality\": {}", self.quality as u8) // rempve " to avoid conflicts in json file
    }

    fn tags(&self) -> Option<Tags> { None }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.ids) + mem::size_of_val(&*self.map_ids)
    }

    fn sum(&self) -> Option<usize> { None }

    fn ids(&self) -> Option<&[ID]> {
        Some(&self.ids)
    }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn sample_count(&self) -> Option<usize> { None }

    fn edge_mults(&self) -> Option<&EdgeMult> {
        Some(&self.edge_mults)
    }

    fn quality(&self) -> Option<BaseQuality> {
        Some(self.quality)
    }

    fn fix_edge_mults(&mut self, exts: Exts) {
        self.edge_mults.clean_edges(exts);
    }

    fn set_edge_mults(&mut self, edge_mults: Option<EdgeMult>) {
        self.edge_mults = edge_mults.expect("Error: no edge mults")
    }

    fn mapped_ids(&self) -> Option<&[ID]> {
        Some(&self.map_ids)
    }

    fn set_mapped_ids(&mut self, mapped_ids: Box<[ID]>) {
        self.map_ids = mapped_ids
    }

    fn join_test(&self, other: &Self) -> bool {
        self.ids == other.ids 
        && self.map_ids == other.map_ids
        && self.quality == other.quality
    }

    fn valid(&self, config: &SummaryConfig) -> bool { 
        self.quality >= config.min_quality
    }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, IDTag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let summary = summarize_tags_ids_edge_q(items, config);
        
        let quality = summary.highest_quality.expect("missing quality score - required for summarizer");
        
        // caluclate p-value with chosen test
        let valid_p = valid_p(PInfo::Calculate { tag_vec: &summary.tag_vec, tag_counts: &summary.tag_counts}, config);
        let valid_q = quality >= config.min_quality;

        let ids = summary.id_vec.into();
        let tags = Tags::from_tag_vec(summary.tag_vec);

        let valid = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;

        (valid, summary.all_exts, IDMapEMQualityData { ids, map_ids: Vec::new().into(), edge_mults: summary.edge_mults, quality }) 
    }

    fn summarizer() -> Summarizers {
        Summarizers::IDMapEMQuality
    }
}

/// Implementation of [`SummaryData<Tag>`]
/// 
/// Contains the IDs the k-mer was observed with, a placeholder for mapped ids, and edge multiplicites/coverage
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SumMapEMQualityData {
    sum: u32,
    map_ids: Box<[ID]>,
    edge_mults: EdgeMult,
    quality: BaseQuality
}

impl SummaryData<Tag> for SumMapEMQualityData{
    fn print(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let map_ids_format = id_format(&self.map_ids, translator, id_group_translator);

        format!("sum: {}, mapped IDs: {}, quality: {}, edge coverage: {}", 
            self.sum, 
            map_ids_format,
            self.quality,
            self.edge_mults,
        ).replace("\"", "\'") // replace " with ' to avoid conflicts in dot file
    }

    fn print_ol(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let map_ids_format = id_format(&self.map_ids, translator, id_group_translator);

        format!("sum: {}, mapped IDs: {}, quality: {}, edge coverage: {:?}", 
            self.sum, 
            map_ids_format,
            self.quality,
            self.edge_mults
        ).replace("\"", "\'") // replace " with ' to avoid conflicts in dot file
    }

    fn print_json(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let map_ids_format = id_format(&self.map_ids, translator, id_group_translator);

        let has_mapped = !self.map_ids.is_empty() as usize;

        format!("\"sum\": {}, \"mapped_ids\": {map_ids_format}, \"has_mapped_ids\": {has_mapped}, \"quality\": {}", self.sum, self.quality as u8) // rempve " to avoid conflicts in json file
    }

    fn tags(&self) -> Option<Tags> { None }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.map_ids)
    }

    fn sum(&self) -> Option<usize> { 
        Some(self.sum as usize)
    }

    fn ids(&self) -> Option<&[ID]> { None }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn sample_count(&self) -> Option<usize> { None }

    fn edge_mults(&self) -> Option<&EdgeMult> {
        Some(&self.edge_mults)
    }

    fn quality(&self) -> Option<BaseQuality> {
        Some(self.quality)
    }

    fn fix_edge_mults(&mut self, exts: Exts) {
        self.edge_mults.clean_edges(exts);
    }

    fn set_edge_mults(&mut self, edge_mults: Option<EdgeMult>) {
        self.edge_mults = edge_mults.expect("Error: no edge mults")
    }

    fn mapped_ids(&self) -> Option<&[ID]> {
        Some(&self.map_ids)
    }

    fn set_mapped_ids(&mut self, mapped_ids: Box<[ID]>) {
        self.map_ids = mapped_ids
    }

    fn join_test(&self, other: &Self) -> bool {
        self.sum == other.sum 
        && self.map_ids == other.map_ids
        && self.quality == other.quality
    }

    fn valid(&self, config: &SummaryConfig) -> bool { 
        self.quality >= config.min_quality
        && self.sum as usize >= config.min_kmer_obs
    }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, Tag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let summary = summarize_tags_edge_q(items, config);
        
        let quality = summary.highest_quality.expect("missing quality score - required for summarizer");
        
        // caluclate p-value with chosen test
        let valid_p = valid_p(PInfo::Calculate { tag_vec: &summary.tag_vec, tag_counts: &summary.tag_counts}, config);
        let valid_q = quality >= config.min_quality;

        let tags = Tags::from_tag_vec(summary.tag_vec);

        let valid = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;

        let sum = match config.significant {
            Some(digits) => round_digits(summary.sum, digits),
            None => summary.sum  
        };

        (valid, summary.all_exts, SumMapEMQualityData { sum, map_ids: Vec::new().into(), edge_mults: summary.edge_mults, quality }) 
    }

    fn summarizer() -> Summarizers {
        Summarizers::SumMapEMQuality
    }
}

/// Implementation of [`SummaryData<Tag>`]
/// 
/// Contains the IDs the k-mer was observed with, a placeholder for mapped ids, and edge multiplicites/coverage
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MapEMQualityData {
    map_ids: Box<[ID]>,
    edge_mults: EdgeMult,
    quality: BaseQuality
}

impl SummaryData<Tag> for MapEMQualityData{
    fn print(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let map_ids_format = id_format(&self.map_ids, translator, id_group_translator);

        format!("mapped IDs: {}, quality: {}, edge coverage: {}", 
            map_ids_format,
            self.quality,
            self.edge_mults,
        ).replace("\"", "\'") // replace " with ' to avoid conflicts in dot file
    }

    fn print_ol(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let map_ids_format = id_format(&self.map_ids, translator, id_group_translator);

        format!("mapped IDs: {}, quality: {}, edge coverage: {:?}", 
            map_ids_format,
            self.quality,
            self.edge_mults
        ).replace("\"", "\'") // replace " with ' to avoid conflicts in dot file
    }

    fn print_json(&self, translator: &Translator, _: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
        let map_ids_format = id_format(&self.map_ids, translator, id_group_translator);

        let has_mapped = !self.map_ids.is_empty() as usize;

        format!("\"mapped_ids\": {map_ids_format}, \"has_mapped_ids\": {has_mapped}, \"quality\": {}", self.quality as u8) // rempve " to avoid conflicts in json file
    }

    fn tags(&self) -> Option<Tags> { None }

    fn mem(&self) -> usize {
        mem::size_of_val(self) + mem::size_of_val(&*self.map_ids)
    }

    fn sum(&self) -> Option<usize> { None }

    fn ids(&self) -> Option<&[ID]> { None }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn sample_count(&self) -> Option<usize> { None }

    fn edge_mults(&self) -> Option<&EdgeMult> {
        Some(&self.edge_mults)
    }

    fn quality(&self) -> Option<BaseQuality> {
        Some(self.quality)
    }

    fn fix_edge_mults(&mut self, exts: Exts) {
        self.edge_mults.clean_edges(exts);
    }

    fn set_edge_mults(&mut self, edge_mults: Option<EdgeMult>) {
        self.edge_mults = edge_mults.expect("Error: no edge mults")
    }

    fn mapped_ids(&self) -> Option<&[ID]> {
        Some(&self.map_ids)
    }

    fn set_mapped_ids(&mut self, mapped_ids: Box<[ID]>) {
        self.map_ids = mapped_ids
    }

    fn join_test(&self, other: &Self) -> bool {
        self.map_ids == other.map_ids
        && self.quality == other.quality
    }

    fn valid(&self, config: &SummaryConfig) -> bool { 
        self.quality >= config.min_quality
    }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, Tag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let summary = summarize_tags_edge_q(items, config);
        
        let quality = summary.highest_quality.expect("missing quality score - required for summarizer");
        
        // caluclate p-value with chosen test
        let valid_p = valid_p(PInfo::Calculate { tag_vec: &summary.tag_vec, tag_counts: &summary.tag_counts}, config);
        let valid_q = quality >= config.min_quality;

        let tags = Tags::from_tag_vec(summary.tag_vec);

        let valid = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;

        let sum = match config.significant {
            Some(digits) => round_digits(summary.sum, digits),
            None => summary.sum  
        };

        (valid, summary.all_exts, MapEMQualityData { map_ids: Vec::new().into(), edge_mults: summary.edge_mults, quality }) 
    }

    fn summarizer() -> Summarizers {
        Summarizers::MapEMQuality
    }
}


/// Implementation of [`SummaryData<Tag>`]
/// 
/// Contains how many times the k-mer was observed in each group, only validates count
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

impl SummaryData<Tag> for GroupCountData {
    fn print(&self, _: &Translator, _: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        format!("count 1: {}\ncount 2: {}", self.group1, self.group2)
    }

    fn print_ol(&self, _: &Translator, _: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        format!("count 1: {}, count 2: {}", self.group1, self.group2)
    }

    fn print_json(&self, _: &Translator, _: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        format!("\"count1\": {}, \"count2\": {}", self.group1, self.group2)
    }

    fn tags(&self) -> Option<Tags> { None }

    fn mem(&self) -> usize {
        mem::size_of_val(self)
    }

    fn sum(&self) -> Option<usize> {
        Some((self.group1 + self.group2) as usize)
    }

    fn ids(&self) -> Option<&[ID]> { None }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn sample_count(&self) -> Option<usize> { None }

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn quality(&self) -> Option<BaseQuality> { None }

    fn fix_edge_mults(&mut self, _: Exts) { }
    
    fn set_edge_mults(&mut self, _: Option<EdgeMult>) { }

    fn mapped_ids(&self) -> Option<&[ID]> { None }

    fn set_mapped_ids(&mut self, _: Box<[ID]>) { }

    fn join_test(&self, other: &Self) -> bool {
        self == other
    }

    fn valid(&self, config: &SummaryConfig) -> bool {
        self.sum() >= config.min_kmer_obs as u32
    }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, Tag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let mut all_exts = Exts::empty();
        let mut count1 = 0;
        let mut count2 = 0;

        let mut nobs = 0u32;
        for item in items {
            let tag = (2 as Marker).pow(item.data as u32);
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
            all_exts = all_exts.add(item.exts);
        }

        let (group1, group2) = match config.significant {
            Some(digits) => (round_digits(count1, digits), round_digits(count2, digits)),
            None => (count1, count2)
        };

        assert_eq!((count1 + count2),nobs);
        (nobs as usize >= config.min_kmer_obs, all_exts, GroupCountData { group1, group2 })
    }

    fn summarizer() -> Summarizers {
        Summarizers::GroupCount
    }
}

/// Implementation of [`SummaryData<Tag>`]
/// 
/// Contains the relative number of observations for the k-mer (in percent) 
/// in group 1 and the absolute overall count, only validates count
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RelCountData {
    percent: u32,
    count: u32
}

impl SummaryData<Tag> for RelCountData {    
    fn print(&self, _: &Translator, _: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        format!("relative amount group 1: {}\ncount both: {}", self.percent, self.count)
    }

    fn print_ol(&self, _: &Translator, _: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        format!("relative amount group 1: {}, count both: {}", self.percent, self.count)
    }

    fn print_json(&self, _: &Translator, _: &SummaryConfig, _: Option<&HashMap<ID, ID>>) -> String {
        format!("\"rel_count_1\": {}, \"sum\": {}", self.percent, self.count)
    }

    fn tags(&self) -> Option<Tags> { None }

    fn mem(&self) -> usize {
        mem::size_of_val(self)
    }

    fn sum(&self) -> Option<usize> {
        Some(self.count as usize)
    }

    fn ids(&self) -> Option<&[ID]> { None }

    fn p_value(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn fold_change(&self, _: &SummaryConfig) -> Option<f32> { None }

    fn sample_count(&self) -> Option<usize> { None }

    fn edge_mults(&self) -> Option<&EdgeMult> { None }

    fn quality(&self) -> Option<BaseQuality> { None }

    fn fix_edge_mults(&mut self, _: Exts) { }
    
    fn set_edge_mults(&mut self, _: Option<EdgeMult>) { }

    fn mapped_ids(&self) -> Option<&[ID]> { None }

    fn set_mapped_ids(&mut self, _: Box<[ID]>) { }

    fn join_test(&self, other: &Self) -> bool {
        self == other
    }
    
    fn valid(&self, config: &SummaryConfig) -> bool {
        self.count >= config.min_kmer_obs as u32
    }

    fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, Tag>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
        let mut all_exts = Exts::empty();
        let mut count1 = 0;
        let mut count2 = 0;

        let mut nobs = 0u32;
        for item in items {
            let tag = (2 as Marker).pow(item.data as u32);
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
            all_exts = all_exts.add(item.exts);
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
    
    fn summarizer() -> Summarizers {
        Summarizers::RelCount
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
pub enum Summarizers {
    Sum,
    VecTags,
    ID,
    IDSum,
    Tags,
    TagsSum,
    TagsCounts,
    TagsCountsSum,
    TagsCountsP,
    TagsCountsEM,
    TagsCountsPEM,
    TagsCountsPEMQuality,
    IDTagsCounts,
    IDTagsCountsPEM,
    IDEM,
    IDMapEM,
    IDMapEMQuality,
    SumMapEMQuality,
    MapEMQuality,
    GroupCount,
    RelCount
}

#[cfg(test)]
mod test {
    

    use bimap::BiMap;
    use crate::{Tags, clean_graph::CleanGraph, compression::{ CheckCompress, ScmapCompress, compress_graph, compress_kmers_with_hash}, dna_string::DnaString, filter::filter_kmers, graph::{BaseGraph, Node}, kmer::{Kmer8, Kmer16}, reads::{Reads, ReadsPaired}, summarizer::{self, ID, IDData, IDTag, NotEnoughSamplesError, SampleInfo, SummaryData, TagsData, Translator, id_format, p_value, students_t_test, u_test, valid_p, welchs_t_test}};
    use crate::Exts;

    use super::{log2_fold_change, round_digits, SummaryConfig, TagsCountsSumData};

    #[test]
    fn test_p_value() -> Result<(), NotEnoughSamplesError> {

        /*
        group 1: 111111100000 = 4064
        group 2: 000000011111 = 31
         */

        let sample_kmers = vec![1; 12];
        let sample_info = SampleInfo::new(31, 4064, sample_kmers);

        let summary_config_w = SummaryConfig::new(sample_info.clone()).with_stat_test(summarizer::StatTest::WelchsTTest);
        let summary_config_t = SummaryConfig::new(sample_info.clone()).with_stat_test(summarizer::StatTest::StudentsTTest);
        let summary_config_u = SummaryConfig::new(sample_info.clone()).with_stat_test(summarizer::StatTest::UTest);
        
        let tag_vec = [0, 1, 2, 3, 4, 8];
        let tag_counts = vec![1; 6];
        
        let p = welchs_t_test(&tag_vec, &tag_counts, &sample_info)?;
        assert_eq!((p * 1000.).round(), 1.);
        let p = students_t_test(&tag_vec, &tag_counts, &sample_info)?;
        assert_eq!((p * 10000.).round(), 5.);
        let p = u_test(&tag_vec, &tag_counts, &sample_info)?;
        assert_eq!((p * 1000.).round(), 5.);

        let p: f32 = p_value(&tag_vec, &tag_counts, &summary_config_w)?;
        assert_eq!((p * 1000.).round(), 1.);
        let p = p_value(&tag_vec, &tag_counts, &summary_config_t)?;
        assert_eq!((p * 10000.).round(), 5.);
        let p = p_value(&tag_vec, &tag_counts, &summary_config_u)?;
        assert_eq!((p * 1000.).round(), 5.);


        // test with different kmer counts
        let sample_kmers = vec![12, 3345, 3478, 87, 1, 2, 666, 98111, 23982938, 555, 122, 7238];

        let sample_info = SampleInfo::new(31, 4064, sample_kmers);
        let summary_config_w = SummaryConfig::new(sample_info.clone()).with_stat_test(summarizer::StatTest::WelchsTTest);
        let summary_config_t = SummaryConfig::new(sample_info.clone()).with_stat_test(summarizer::StatTest::StudentsTTest);
        let summary_config_u = SummaryConfig::new(sample_info.clone()).with_stat_test(summarizer::StatTest::UTest);

        let tag_vec = [0, 1, 7, 8, 9, 10];
        let tag_counts = vec![1; 6];

        let p = welchs_t_test(&tag_vec, &tag_counts, &sample_info)?;
        assert_eq!((p * 10000.).round(), 4113.);
        let p = students_t_test(&tag_vec, &tag_counts, &sample_info)?;
        assert_eq!((p * 10000.).round(), 2955.);
        let p = u_test(&tag_vec, &tag_counts, &sample_info)?;
        assert_eq!((p * 1000.).round(), 862.);

        let p: f32 = p_value(&tag_vec, &tag_counts, &summary_config_w)?;
        assert_eq!((p * 10000.).round(), 4113.);
        let p: f32 = p_value(&tag_vec, &tag_counts, &summary_config_t)?;
        assert_eq!((p * 10000.).round(), 2955.);
        let p: f32 = p_value(&tag_vec, &tag_counts, &summary_config_u)?;
        assert_eq!((p * 1000.).round(), 862.);


        // test: not enough samples for statistical test

        /*
        group 1: 000000100 = 4
        group 2: 000000011 = 3
         */

        let sample_kmers = vec![3, 3, 3];
        let sample_info = SampleInfo::new(0b100, 0b11, sample_kmers);
        let summary_config_w = SummaryConfig::new(sample_info.clone()).with_stat_test(summarizer::StatTest::WelchsTTest);
        let summary_config_t = SummaryConfig::new(sample_info.clone()).with_stat_test(summarizer::StatTest::StudentsTTest);
        let summary_config_u = SummaryConfig::new(sample_info.clone()).with_stat_test(summarizer::StatTest::UTest);

        let tag_vec = [0, 1, 2];
        let tag_counts = vec![2, 3, 1];

        if welchs_t_test(&tag_vec, &tag_counts, &sample_info).is_ok() { panic!("should throw err") }
        if students_t_test(&tag_vec, &tag_counts, &sample_info).is_ok() { panic!("should throw err") }
        if u_test(&tag_vec, &tag_counts, &sample_info).is_ok() { panic!("should throw err") }

        if p_value(&tag_vec, &tag_counts, &summary_config_w).is_ok() { panic!("should throw err") }
        if p_value(&tag_vec, &tag_counts, &summary_config_t).is_ok() { panic!("should throw err") }
        if p_value(&tag_vec, &tag_counts, &summary_config_u).is_ok() { panic!("should throw err") }

        Ok(())
    }

    #[test]
    fn test_valid_p() {
        let sample_kmers = vec![1; 12];
        let sample_info = SampleInfo::new(31, 4064, sample_kmers);
        let summary_config_m = SummaryConfig::new(sample_info.clone()).with_max_p(None);
        let summary_config_p = SummaryConfig::new(sample_info.clone()).with_max_p(Some(0.05));

        let tag_vec = [0, 1, 2, 3, 4, 8];  // p should be 0.001
        let tag_counts = vec![1; 6];
        let vp = valid_p(summarizer::PInfo::Calculate { tag_vec: &tag_vec, tag_counts: &tag_counts }, &summary_config_m);
        assert!(vp);
        let vp = valid_p(summarizer::PInfo::Calculate { tag_vec: &tag_vec, tag_counts: &tag_counts }, &summary_config_p);
        assert!(vp);


        let tag_vec = [0, 7, 8, 9, 10]; // p should be 0.2238
        let tag_counts = vec![1; 5];

        let vp = valid_p(summarizer::PInfo::Calculate { tag_vec: &tag_vec, tag_counts: &tag_counts }, &summary_config_m);
        assert!(vp);
        let vp = valid_p(summarizer::PInfo::Calculate { tag_vec: &tag_vec, tag_counts: &tag_counts }, &summary_config_p);
        assert!(!vp);
    }

    #[test]
    #[cfg(not(feature = "sample128"))]
    fn test_data_valid() {
        use crate::{kmer::Kmer16, test::build_test_graph};

        let mut graph: BaseGraph<Kmer8, TagsCountsSumData> = BaseGraph::new(false);

        let tags = Tags::from_tag_vec(vec![0, 2, 6]);
        let counts: Box<[u32]> = [1, 3, 5].into();
        let sum = counts.iter().sum::<u32>();
        graph.add(&DnaString::from_acgt_bytes("AAAAAAAA".as_bytes()), Exts::empty(), TagsCountsSumData { tags, counts, sum });

        let tags = Tags::from_tag_vec(vec![0]);
        let counts: Box<[u32]> = [1].into();
        let sum = counts.iter().sum::<u32>();
        graph.add(&DnaString::from_acgt_bytes("CCCCCCCC".as_bytes()), Exts::empty(), TagsCountsSumData { tags, counts, sum });
                
        
        let graph = graph.finish();

        graph.print();

        let sample_kmers = vec![123, 234, 12334, 34, 1232, 123, 123, 34];
        let sample_info = SampleInfo::new(0b00100101, 0b11011010, sample_kmers);
        let config = SummaryConfig::new(sample_info).with_stat_test(summarizer::StatTest::StudentsTTest);

        let censor_nodes = CleanGraph::new(|node: &Node<'_, Kmer8, TagsCountsSumData>| !node.data().valid(&config))
                    .find_bad_nodes(&graph);
        println!("censor nodes: {:?}", censor_nodes);
        let filter_graph = compress_graph(false, &ScmapCompress::new(), graph, Some(censor_nodes));

        filter_graph.print();

        // larger test

        let (_, _, ser_graph) = build_test_graph::<Kmer16, TagsCountsSumData, _>();
        let (graph, _translator, mut config) = ser_graph.dissolve();

        config.set_min_kmer_obs(3);

        let node3 = graph.get_node(3);
        assert!(!node3.data().valid(&config));

        let bad_nodes = graph.find_bad_nodes(|node| node.data().valid(&config));
        let bad_node_correct = vec![3, 15, 20, 21, 23, 35, 36, 37, 40, 41, 46, 58, 61, 66, 67, 76, 77, 80, 89, 98, 99];
        assert_eq!(bad_nodes, bad_node_correct);
        let _filtered_graph = compress_graph(false, &ScmapCompress::new(), graph, Some(bad_nodes));

    }

    #[test]
    fn test_fold_change() {    
        let marker0 = 0b0000000001111;
        let marker1 = 0b1111111110000;
        let sample_kmers = vec![2834, 2343, 12, 1234, 345345, 122, 234, 23455, 231, 2, 3564, 12344, 34555];
        let sample_info = SampleInfo::new(marker0, marker1, sample_kmers);
        //let summary_config = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), None, summarizer::StatTest::WelchsTTest);

        let labels = vec![0, 1, 2, 3, 7, 8];
        let tags = Tags::from_tag_vec(labels);
        let counts = vec![1, 6, 9, 3, 6, 10];
        let fold_change = log2_fold_change(tags, &counts, &sample_info);
        assert_eq!(fold_change, 5.286_453_2);

        let labels = vec![0, 6, 7, 8, 10, 11];
        let tags = Tags::from_tag_vec(labels);
        let counts = vec![12, 3, 7, 1, 22, 6];
        let fold_change = log2_fold_change(tags, &counts, &sample_info);
        assert_eq!(fold_change, -1.339_324_5);

        // x/0 = inf -> log(inf) = inf
        let labels = vec![0, 1];
        let tags = Tags::from_tag_vec(labels);
        let counts = vec![12, 3];
        let fold_change = log2_fold_change(tags, &counts, &sample_info);
        assert_eq!(fold_change, f32::INFINITY);

        // 0/x = 0 -> log2(0) = -inf
        let labels = vec![7, 8];
        let tags = Tags::from_tag_vec(labels);
        let counts = vec![12, 3];
        let fold_change = log2_fold_change(tags, &counts, &sample_info);
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


    #[test]
    fn test_id_format() {
        let id_translator = [("A", 0), ("B", 1), ("C", 2), ("D", 3), ("E", 4), ("F", 5), ("G", 6)].into_iter().map(|(a, b)| (a.to_string(), b as ID)).collect();
        let id_gr_tr = [(0 as ID, 0 as ID), (1, 0), (2, 0), (3, 0), (4, 1), (5, 1), (6, 1)].into_iter().collect();

        let translator = Translator::new_id_translator(id_translator);
        let e_tr = Translator::new_tag_translator(BiMap::new());

        assert_eq!("[\"A\", \"B\", \"C\"]", id_format(&[0, 1, 2], &translator, None));
        assert_eq!("[0, 1]", id_format(&[0, 1, 5], &translator, Some(&id_gr_tr)));
        assert_eq!("[0, 1, 2]", id_format(&[0, 1, 2], &e_tr, None));
    }

    #[test]
    fn test_summary_config() {
        let sample_info = SampleInfo::new(0b11, 0b1100, vec![12, 12, 12, 12]);
        let config1 = SummaryConfig::new(sample_info.clone())
            .with_group_frac(summarizer::GroupFrac::One, 0.3)
            .with_max_p(Some(0.3))
            .with_min_kmer_obs(2)
            .with_min_quality(crate::BaseQuality::Marginal)
            .with_min_quality_for_edge(crate::BaseQuality::Medium)
            .with_significant(Some(4))
            .with_stat_test(summarizer::StatTest::StudentsTTest);

        let config2 = SummaryConfig {
            min_kmer_obs: 2,
            significant: Some(4),
            group_frac: summarizer::GroupFrac::One,
            frac_cutoff: 0.3,
            sample_info,
            max_p: Some(0.3),
            stat_test: summarizer::StatTest::StudentsTTest,
            stat_test_changed: false,
            min_quality: crate::BaseQuality::Marginal,
            min_quality_for_edge: crate::BaseQuality::Medium
        };

        assert_eq!(config1, config2)
    }

    #[test]
    fn test_summarize_edge_q() {
        let read1 = "AGCTAGCGATGCTAGCTAGCATCGTAGCTAGCAAGCTGATCAAGTCGATGCTGACTGATGCTAGCTGACTGATCGATGCTAGCTGATC";
        let qual1 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
        let read2 = "GAAATCGTAGCTGAAAAATCGTAGCTGATCGTAGCTGATCTATTCTAGCTGATCAAGTCGATCGGAGGGGTTTCGGAGTTTCGGGATTCGTAT";
        let qual2 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";

        /*
        should produce three sequences, not connected to each other:
        -> the first read as a whole and the second read, with the false k-mer missing and producing a gap
        seq: AGCTAGCGATGCTAGCTAGCATCGTAGCTAGCAAGCTGATCAAGTCGATGCTGACTGATGCTAGCTGACTGATCGATGCTAGCTGATC, node: Node { id:0, Exts: |, L:[] R:[], Seq: 88, Data: TagsData { tags: [0] } }
        seq: GAAATCGTAGCTGAAAAATCGTAGCTGATCGTAGCTGATCTATTCTAGCTGATCAAGTCGA, node: Node { id:1, Exts: |, L:[] R:[], Seq: 61, Data: TagsData { tags: [1] } }
        seq: GCTGATCAAGTCGATCGGAGGGGTTTCGGAGTTTCGGGATTCGTAT, node: Node { id:2, Exts: |, L:[] R:[], Seq: 46, Data: TagsData { tags: [1] } }
        [0]
        [1]
        [2]
        */

        // without IDs

        let mut reads = Reads::new_with_quality(crate::reads::Strandedness::Forward);
        reads.add_read(DnaString::from_acgt_bytes(read1.as_bytes()), None, 0, Some(qual1.as_bytes()));
        reads.add_read(DnaString::from_acgt_bytes(read2.as_bytes()), None, 1, Some(qual2.as_bytes()));

        let reads_paired = ReadsPaired::Unpaired { reads };

        let sample_info = SampleInfo::new(0b10, 0b01, vec![73, 78]);
        let summary_config = SummaryConfig::new(sample_info)
            .with_min_quality_for_edge(crate::BaseQuality::Medium);

        let (kmers, _) = filter_kmers::<TagsData, Kmer16, _>(&reads_paired, &summary_config, false, 1., false);

        let comp_spec = CheckCompress::new(|a: TagsData, _b| a, |a, b| a.join_test(b));
        let mut graph = compress_kmers_with_hash(true, &comp_spec, kmers, false, false).finish();
        graph.fix_exts(None);

        for node_id in 0..graph.len() {
            let node = graph.get_node(node_id);
            println!{"seq: {:?}, node: {:?}", node.sequence(), node}
        }

        let expected_components = [vec![0], vec![1], vec![2]];

        for (i, component) in graph.iter_components().enumerate() {
            println!("{:?}", component);
            assert_eq!(component, expected_components[i]);
        }

        // with IDs

        let mut reads = Reads::new_with_quality(crate::reads::Strandedness::Forward);
        reads.add_read(DnaString::from_acgt_bytes(read1.as_bytes()), None, IDTag::new(0, 0), Some(qual1.as_bytes()));
        reads.add_read(DnaString::from_acgt_bytes(read2.as_bytes()), None, IDTag::new(1, 1), Some(qual2.as_bytes()));

        let reads_paired = ReadsPaired::Unpaired { reads };

        let (kmers, _) = filter_kmers::<IDData, Kmer16, _>(&reads_paired, &summary_config, false, 1., false);

        let comp_spec = CheckCompress::new(|a: IDData, _b| a, |a, b| a.join_test(b));
        let mut graph = compress_kmers_with_hash(true, &comp_spec, kmers, false, false).finish();
        graph.fix_exts(None);

        for node_id in 0..graph.len() {
            let node = graph.get_node(node_id);
            println!{"seq: {:?}, node: {:?}", node.sequence(), node}
        }

        let expected_components = [vec![0], vec![1], vec![2]];

        for (i, component) in graph.iter_components().enumerate() {
            println!("{:?}", component);
            assert_eq!(component, expected_components[i]);
        }

    }
}