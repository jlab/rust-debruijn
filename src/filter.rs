// Copyright 2017 10x Genomics

//! Methods for converting sequences into kmers, filtering observed kmers before De Bruijn graph construction, and summarizing 'color' annotations.

use core::f32;
use std::fmt::Debug;
use std::marker::PhantomData;
use std::mem;
use std::ops::Deref;
use std::sync::Arc;
use std::sync::Mutex;
use std::time::Instant;

use bimap::BiMap;
use boomphf::hashmap::BoomHashMap2;
use itertools::Itertools;
use log::debug;
use num_traits::abs;
use rayon::current_num_threads;
use rayon::prelude::*;

use crate::Dir;
use crate::Exts;
use crate::Kmer;
use crate::Tags;
use crate::Vmer;

fn bucket<K: Kmer>(kmer: K) -> usize {
    (kmer.get(0) as usize) << 6
        | (kmer.get(1) as usize) << 4
        | (kmer.get(2) as usize) << 2
        | (kmer.get(3) as usize)
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
}

impl<> SummaryData<u16> for u16 {
    fn new(data: u16) -> Self {
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

}

#[derive(Debug, Clone, PartialEq)]
pub struct TagsSumData {
    tags: Tags,
    sum: i32,
}

impl SummaryData<(Tags, i32)> for TagsSumData {
    fn new(data: (Tags, i32)) -> Self {
        TagsSumData { tags: data.0, sum: data.1 }
    }

    fn print(&self, tag_translator: &BiMap<&str, u8>) -> String {
        format!("tags: {:?}, sum: {}", self.tags.to_string_vec(tag_translator), self.sum).replace("\"", "\'")
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
}

#[derive(Debug, Clone, PartialEq)]
pub struct TagsCountData {
    tags: Tags,
    counts: Box<[u32]>,
    sum: i32,
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
    
    fn new(min_kmer_obs: usize) -> Self;
    fn summarize<K: Kmer, F: Iterator<Item = (K, Exts, DI)>>(&self, items: F) -> (bool, Exts, DO);
}

/// A simple KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times. The metadata returned about a Kmer
/// is the number of times it was observed, capped at 2^16.
pub struct CountFilter<D> {
    min_kmer_obs: usize,
    phantom: PhantomData<D>
}

impl<D> KmerSummarizer<D, u16, u16> for CountFilter<D> {
    fn new(min_kmer_obs: usize) -> Self {
        CountFilter {
            min_kmer_obs,
            phantom: PhantomData,
        }
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, D)>>(&self, items: F) -> (bool, Exts, u16) {
        let mut all_exts = Exts::empty();
        let mut count = 0u16;
        for (_, exts, _) in items {
            count = count.saturating_add(1);
            all_exts = all_exts.add(exts);
        }

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
    fn new(min_kmer_obs: usize) -> Self {
        CountFilterSet {
            min_kmer_obs,
            phantom: PhantomData,
        }
    }

    fn summarize<K, F: Iterator<Item = (K, Exts, D)>>(&self, items: F) -> (bool, Exts, Vec<D>) {
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
    pub sum_datas: u64,
    phantom: PhantomData<u8>,
}

impl KmerSummarizer<u8, TagsSumData, (Tags, i32)> for CountFilterComb {
    fn new(min_kmer_obs: usize) -> Self {
        CountFilterComb {
            min_kmer_obs,
            sum_datas: 0,
            phantom: PhantomData,
        }
    }

    fn summarize<K: Kmer, F: Iterator<Item = (K, Exts, u8)>>(&self, items: F) -> (bool, Exts, TagsSumData) {
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
    fn new(min_kmer_obs: usize) -> Self {
        CountFilterStats {
            min_kmer_obs,
            phantom: PhantomData,
        }
    }

    fn summarize<K: Kmer, F: Iterator<Item = (K, Exts, u8)>>(&self, items: F) -> (bool, Exts, TagsCountData) {
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

pub enum Summarizer {
    Count,
    Tags,
    CountNTags,
    CountNCountTags,
}



/// Process DNA sequences into kmers and determine the set of valid kmers,
/// their extensions, and summarize associated label/'color' data. The input
/// sequences are converted to kmers of type `K`, and like kmers are grouped together.
/// All instances of each kmer, along with their label data are passed to
/// `summarizer`, an implementation of the `KmerSummarizer` which decides if
/// the kmer is 'valid' by an arbitrary predicate of the kmer data, and
/// summarizes the the individual label into a single label data structure
/// for the kmer. Care is taken to keep the memory consumption small.
/// Less than 4G of temporary memory should be allocated to hold intermediate kmers.
///
///
/// # Arguments
///
/// * `seqs` a slice of (sequence, extensions, data) tuples. Each tuple
///   represents an input sequence. The input sequence must implement `Vmer<K`> The data slot is an arbitrary data
///   structure labeling the input sequence.
///   If complete sequences are passed in, the extensions entry should be
///   set to `Exts::empty()`.
///   In sharded DBG construction (for example when minimizer-based partitioning
///   of the input strings), the input sequence is a sub-string of the original input string.
///   In this case the extensions of the sub-string in the original string
///   should be passed in the extensions.
/// * `summarizer` is an implementation of `KmerSummarizer<D1,DS>` that decides
///   whether a kmer is valid (e.g. based on the number of observation of the kmer),
///   and summarizes the data about the individual kmer observations. See `CountFilter`
///   and `CountFilterSet` for examples.
/// * `stranded`: if true, preserve the strandedness of the input sequences, effectively
///   assuming they are all in the positive strand. If false, the kmers will be canonicalized
///   to the lexicographic minimum of the kmer and it's reverse complement.
/// * `report_all_kmers`: if true returns the vector of all the observed kmers and performs the
///   kmer based filtering
/// * `memory_size`: gives the size bound on the memory in GB to use and automatically determines
///   the number of passes needed.
/// # Returns
/// BoomHashMap2 Object, check rust-boomphf for details
#[inline(never)]
//pub fn filter_kmers_parallel<K: Kmer + Sync + Send, V: Vmer + Sync, D1: Clone + Debug + Sync, DS: Clone + Sync + Send, S: KmerSummarizer<D1, DS, (usize, usize)> +  Send>(
pub fn filter_kmers_parallel<K: Kmer + Sync + Send, V: Vmer + Sync, DO, DS: Clone + std::fmt::Debug + Send + SummaryData<DO>, S: KmerSummarizer<u8, DS, DO>>(
    seqs: &[(V, Exts, u8)],
    //summarizer: &dyn Deref<Target = S>,
    // summarizer without wrapper, why wrapper???
    _summarizer: Box<S>,
    //summarizer: Summarizer,
    min_kmer_obs: usize,
    stranded: bool,
    report_all_kmers: bool,
    memory_size: usize,
    time: bool,
    progress: bool,
) -> (BoomHashMap2<K, Exts, DS>, Vec<K>)
{
    // take timestamp before all processes
    let before_all = Instant::now();

    let rc_norm = !stranded;
    const BUCKETS: usize = 256;

    // Estimate 6 consumed by Kmer vectors, and set iteration count appropriately
    let input_kmers: usize = seqs
        .iter()
        .map(|&(ref vmer, _, _)| vmer.len().saturating_sub(K::k() - 1))
        .sum();

    if time { println!("time counting kmers (s): {}", before_all.elapsed().as_secs_f32()) }

    let kmer_mem = input_kmers * mem::size_of::<(K, Exts, u8)>();
    let max_mem = memory_size * 10_usize.pow(9);
    let slices_seq = kmer_mem / max_mem + 1;
    let slices = slices_seq;
    //let sz = buckets / slices + 1;

    debug!("kmers: {}, mem per kmer: {}, kmer_mem: {} Bytes, slices: {}", input_kmers, mem::size_of::<(K, Exts, u8)>(), kmer_mem, slices);
    
    // exponential probabilty distribution: p(x) = lambda * exp(-lambda * x)
    let lambda = 0.008;

    let mut bucket_ranges = Vec::with_capacity(if slices < BUCKETS {BUCKETS} else {slices});

    for i in 1..=slices {
        // calculate lower and upper bound with Quantile function of exponential probability distribution
        // I(i) = [1/lambda * |ln(1 - (i-1)/slices)|, 1/lambda * |ln(1 - i/slices)|]
        let lbound = ((1./lambda) * abs((1.-(i as f32 - 1.)/slices as f32).ln())) as usize;
        let ubound = ((1./lambda) * abs((1.-i as f32/slices as f32).ln())) as usize;

        // if upper bound is above no of buckets (256), reduce to no of buckets
        let ubound = if ubound > BUCKETS { BUCKETS } else { ubound };
        if ubound > lbound && lbound < BUCKETS { bucket_ranges.push(lbound..ubound) };
    }

    debug!("bucket_ranges: {:?}, len br: {}", bucket_ranges, bucket_ranges.len());
    assert!(bucket_ranges[bucket_ranges.len() - 1].end >= BUCKETS);
    let n_buckets = bucket_ranges.len();

    if bucket_ranges.len() > 1 {
        debug!(
            "{} sequences, {} kmers, {} passes",
            seqs.len(),
            input_kmers,
            bucket_ranges.len()
        );
    }

    debug!("n of seqs: {}", seqs.len());

    let mut time_picking_par = 0.;
    let mut time_picking = 0.;
    let mut time_summarizing = 0.;

    let shared_target_vecs = Arc::new(Mutex::new((Vec::new(), Vec::new(), Vec::new(), Vec::new())));

    if time { println!("time all prepariations before sliced in filter_kmers (s): {}", before_all.elapsed().as_secs_f32()) }

    for (i, bucket_range) in bucket_ranges.into_iter().enumerate() {

        debug!("Processing bucket {} of {}", i+1, n_buckets);

        let before_kmer_picking = Instant::now();
        // first step: picking kmers with their exts & data from the reads
        // go through all kmers and sort into bucket according to first four bases
        // all kmers starting with "AAAA" go in kmer_buckets[0], all starting with AAAC go in kmer_buckets[1] and so on
        // when using the first four bases, this needs 256 buckets
        // the buckets are split in to the bucket_ranges to save memory

        // split all reads into ranges to be processed in parallel for counting capacities and picking kmers

        let n_threads = current_num_threads();
        let n_reads = seqs.len();
        let sz = n_reads / n_threads + 1;

        debug!("n_reads: {}", n_reads);
        debug!("sz: {}", sz);

        let mut parallel_ranges = Vec::with_capacity(slices);
        let mut start = 0;
        while start < n_reads {
            parallel_ranges.push(start..start + sz);
            start += sz;
        }

        let last_start = parallel_ranges.pop().expect("no kmers in parallel ranges").start;
        parallel_ranges.push(last_start..n_reads);
        debug!("parallel ranges: {:?}", parallel_ranges);


        let kmer_buckets = Arc::new(Mutex::new(vec![Vec::new(); n_threads]));

        let before_picking_parallel = Instant::now();

        parallel_ranges.clone().into_par_iter().enumerate().for_each(|(i, range)| {

            // first go trough all kmers to find the length of all buckets (to reserve capacity)
            let mut capacities = [0usize; BUCKETS];
            for &(ref seq, _, _) in &seqs[range.clone()] { 
                // iterate through all kmers in seq
                for kmer in seq.iter_kmers::<K>() {
                    // if not stranded choose lexiographically lesser of kmer and rc of kmer
                    let min_kmer = if rc_norm {
                        let (min_kmer, _) = kmer.min_rc_flip();
                        min_kmer
                    } else {
                        kmer
                    };

                    // calculate which bucket this kmer belongs to
                    let bucket = bucket(min_kmer);
                    // check if bucket is in current range and if so, add one to needed capacity
                    let in_range = bucket >= bucket_range.start && bucket < bucket_range.end;
                    if in_range { 
                        capacities[bucket] += 1;
                    }
                }
            }

            let mut kmer_buckets1d = Vec::with_capacity(BUCKETS); 
            // reserve capacities needed for current range in each bucket
            for capacity in capacities.into_iter() {
                kmer_buckets1d.push(Vec::with_capacity(capacity));
            }
            for &(ref seq, seq_exts, ref d) in &seqs[range] {
                for (kmer, exts) in seq.iter_kmer_exts::<K>(seq_exts) {
                    let (min_kmer, flip_exts) = if rc_norm {
                        let (min_kmer, flip) = kmer.min_rc_flip();
                        let flip_exts = if flip { exts.rc() } else { exts };
                        (min_kmer, flip_exts)
                    } else {
                        (kmer, exts)
                    };

                    // calculate which bucket this kmer belongs to
                    let bucket = bucket(min_kmer);

                    // check if bucket is in current range and if so, push kmer to bucket
                    let in_range = bucket >= bucket_range.start && bucket < bucket_range.end;
                    if in_range {
                        kmer_buckets1d[bucket].push((min_kmer, flip_exts, d.clone()));
                    }
                }
            }

            // clone and lock kmer_buckets to safely share across threads
            let kb_clone = Arc::clone(&kmer_buckets);
            let mut kb2d = kb_clone.lock().expect("lock kmer buckets 2d");
            // replace empty vec too keep order
            kb2d[i] = kmer_buckets1d;

        });

        time_picking_par += before_picking_parallel.elapsed().as_secs_f32();

        // unlock kmer buckets and move out of guard so they can be turned into iterator
        let mut kmer_buckets = kmer_buckets.lock().expect("unlock kmer_buckets final");
        let kmer_buckets = mem::take(&mut *kmer_buckets);

        // all combined buckets go into this vector
        let mut new_buckets = vec![Vec::new(); BUCKETS];
        // flatten kmer buckets
        for thread_vec in kmer_buckets.into_iter() {
            for (i, mut bucket) in thread_vec.into_iter().enumerate() {
                new_buckets[i].reserve_exact(bucket.len());
                new_buckets[i].append(&mut bucket);
            }
        }

        time_picking += before_kmer_picking.elapsed().as_secs_f32();
        
        let before_parallel = Instant::now();

        if progress {
            println!("Processing bucket {} of {}", i+1, n_buckets);
            for _i in 0..127 {
                print!("-");
            }
            print!("|");
            print!("\n");
        }
        
        // parallel start
        new_buckets.into_par_iter().enumerate().for_each(|(j, mut kmer_vec)| {
            //debug!("kmers in bucket #{}: {}", j, kmer_vec.len());
            if progress & (j % 2 == 0) { print!("|") };
            debug!("starting bucket {} with {} kmers, capacity of {}", j, kmer_vec.len(), kmer_vec.capacity());
            kmer_vec.sort_by_key(|elt| elt.0);

            let size = kmer_vec.iter().group_by(|elt| elt.0).into_iter().count();

            let mut all_kmers = Vec::with_capacity(size);
            let mut valid_kmers = Vec::with_capacity(size);
            let mut valid_exts = Vec::with_capacity(size);
            let mut valid_data = Vec::with_capacity(size);


            for (kmer, kmer_obs_iter) in kmer_vec.into_iter().group_by(|elt| elt.0).into_iter() {
                let summarizer = S::new(min_kmer_obs);
                let (is_valid, exts, summary_data) = summarizer.summarize(kmer_obs_iter);
                if report_all_kmers {
                    all_kmers.push(kmer);
                }
                if is_valid {
                    valid_kmers.push(kmer);
                    valid_exts.push(exts);
                    valid_data.push(summary_data);
                }
            }

            // if there are valid k-mers in this bucket, append them to the shared target vectors
            // important that this is done in one step so each kmer has the same index with its exts and data
            if valid_kmers.len() > 0 {

                let stv_clone = Arc::clone(&shared_target_vecs);
                let mut stv = stv_clone.lock().expect("lock target vectors");
                // valid kmers
                stv.0.reserve_exact(valid_kmers.len());
                stv.0.append(&mut valid_kmers); 
                // valid exts
                stv.1.reserve_exact(valid_exts.len());
                stv.1.append(&mut valid_exts); 
                // valid data
                stv.2.reserve_exact(valid_data.len());
                stv.2.append(&mut valid_data); 
            }

            // if kmers were collected into all_kmers, append them to shared target vector
            if all_kmers.len() > 0 {
                let stv_clone = Arc::clone(&shared_target_vecs);
                let mut stv = stv_clone.lock().expect("lock target vectors");
                // all kmers
                stv.3.reserve_exact(all_kmers.len());
                stv.3.append(&mut all_kmers);
            }            
            
            /* let mut data_out = shared_data.lock().expect("unlock shared filter data");
            debug!("bucket {} processed, mem of valid_kmers: {} Bytes, mem of valid_exts: {} Bytes, mem of valid_data: {} Bytes", 
                j, mem::size_of_val(&*valid_kmers), mem::size_of_val(&*valid_exts), mem::size_of_val(&*valid_data));
            // i is the slice/bucket_range
            data_out[i].push((all_kmers, valid_kmers, valid_exts, valid_data)); */
        });
        // parallel end
        if progress { print!("\n") }

        time_summarizing += before_parallel.elapsed().as_secs_f32();

        debug!("processed bucket {i}");
    }

    if time { 
        println!("time counting + collecting par (s): {}", time_picking_par);
        println!("time counting + collecting (s): {}", time_picking);
        println!("time summarizing (s): {}", time_summarizing);
    }

    /* let data_out = shared_data.lock().expect("final unlock shared filter data");
    debug!("data out capacity:: {}, size: {}", data_out[0].capacity(), data_out[0].len()); */

    /* let mut all_kmers = Vec::new();
    let mut valid_kmers = Vec::new();
    let mut valid_exts = Vec::new();
    let mut valid_data = Vec::new();

    for slice in data_out.iter() {
        for bucket in slice {
            all_kmers.reserve_exact(bucket.0.len());
            valid_kmers.reserve_exact(bucket.1.len());
            valid_exts.reserve_exact(bucket.2.len());
            valid_data.reserve_exact(bucket.3.len());

            all_kmers.append(&mut bucket.0.clone());
            valid_kmers.append(&mut bucket.1.clone());
            valid_exts.append(&mut bucket.2.clone());
            valid_data.append(&mut bucket.3.clone());
        }
    }  */

    let stv = shared_target_vecs.lock().expect("final lock target vectors");

    debug!("valid kmers - capacity: {}, size: {}, mem: {}", stv.0.capacity(), stv.0.len(), mem::size_of_val(&*stv.0));
    debug!("valid exts - capacity: {}, size: {}, mem: {}", stv.1.capacity(), stv.1.len(), mem::size_of_val(&*stv.1));
    debug!("valid data - capacity: {}, size: {}, mem: {}", stv.2.capacity(), stv.2.len(), mem::size_of_val(&*stv.2));
    debug!("all kmers - capacity: {}, size: {}, mem: {}", stv.3.capacity(), stv.3.len(), mem::size_of_val(&*stv.3));

    let before_hash = Instant::now();
    let hm = BoomHashMap2::new(stv.0.to_vec(), stv.1.to_vec(), stv.2.to_vec());
    let after_hash = before_hash.elapsed().as_secs_f32();
    let all_kmers = stv.3.to_vec();

    let filter_kmers_inner = before_all.elapsed().as_secs_f32();
    if time { 
        println!("time filter_kmers inner (s): {}", filter_kmers_inner);
        println!("time build filtered hash map (s): {}", after_hash);
    }

    (
        hm,
        all_kmers,
    )
}

/// Process DNA sequences into kmers and determine the set of valid kmers,
/// their extensions, and summarize associated label/'color' data. The input
/// sequences are converted to kmers of type `K`, and like kmers are grouped together.
/// All instances of each kmer, along with their label data are passed to
/// `summarizer`, an implementation of the `KmerSummarizer` which decides if
/// the kmer is 'valid' by an arbitrary predicate of the kmer data, and
/// summarizes the the individual label into a single label data structure
/// for the kmer. Care is taken to keep the memory consumption small.
/// Less than 4G of temporary memory should be allocated to hold intermediate kmers.
///
///
/// # Arguments
///
/// * `seqs` a slice of (sequence, extensions, data) tuples. Each tuple
///   represents an input sequence. The input sequence must implement `Vmer<K`> The data slot is an arbitrary data
///   structure labeling the input sequence.
///   If complete sequences are passed in, the extensions entry should be
///   set to `Exts::empty()`.
///   In sharded DBG construction (for example when minimizer-based partitioning
///   of the input strings), the input sequence is a sub-string of the original input string.
///   In this case the extensions of the sub-string in the original string
///   should be passed in the extensions.
/// * `summarizer` is an implementation of `KmerSummarizer<D1,DS>` that decides
///   whether a kmer is valid (e.g. based on the number of observation of the kmer),
///   and summarizes the data about the individual kmer observations. See `CountFilter`
///   and `CountFilterSet` for examples.
/// * `stranded`: if true, preserve the strandedness of the input sequences, effectively
///   assuming they are all in the positive strand. If false, the kmers will be canonicalized
///   to the lexicographic minimum of the kmer and it's reverse complement.
/// * `report_all_kmers`: if true returns the vector of all the observed kmers and performs the
///   kmer based filtering
/// * `memory_size`: gives the size bound on the memory in GB to use and automatically determines
///   the number of passes needed.
/// # Returns
/// BoomHashMap2 Object, check rust-boomphf for details
#[inline(never)]
pub fn filter_kmers<K: Kmer, V: Vmer, D1: Clone + Debug, DO, DS: SummaryData<DO>, S: KmerSummarizer<D1, DS, DO>>(
    seqs: &[(V, Exts, D1)],
    summarizer: &dyn Deref<Target = S>,
    stranded: bool,
    report_all_kmers: bool,
    memory_size: usize,
    time: bool,
    progress: bool,
) -> (BoomHashMap2<K, Exts, DS>, Vec<K>)
where
    DS: Debug,
{
    let before_all = Instant::now();
    let rc_norm = !stranded;
    const BUCKETS: usize = 256;

    // Estimate memory consumed by Kmer vectors, and set iteration count appropriately
    let input_kmers: usize = seqs
        .iter()
        .map(|&(ref vmer, _, _)| vmer.len().saturating_sub(K::k() - 1))
        .sum();

    if time { println!("time counting kmers (s): {}", before_all.elapsed().as_secs_f32()) }

    let kmer_mem = input_kmers * mem::size_of::<(K, D1)>();
    debug!("size used for calculation: {}B", mem::size_of::<(K, D1)>());
    debug!("size of kmer, E, D: {} B", mem::size_of::<(K, Exts, D1)>());
    debug!("size of K: {} B, size of Exts: {} B, size of D1: {}", mem::size_of::<K>(), mem::size_of::<Exts>(), mem::size_of::<D1>());
    debug!("type D1: {}", std::any::type_name::<D1>());

    let max_mem: usize = memory_size * 10_usize.pow(9);
    let slices: usize = kmer_mem / max_mem + 1;
  
    // exponential probabilty distribution: p(x) = lambda * exp(-lambda * x)
    // for the probability distribution a higher LAMBDA (~ 0.011) would be more accurate 
    // but since vectors with summarized kmers start to grow later, we pretend we have a lower value for LAMBDA 
    const LAMBDA: f32 = 0.008;

    let mut bucket_ranges: Vec<std::ops::Range<usize>> = Vec::with_capacity(if slices < BUCKETS {slices} else {BUCKETS});

    for i in 1..=slices {
        // calculate lower and upper bound with Quantile function of exponential probability distribution
        // I(i) = [1/LAMBDA * |ln(1 - (i-1)/slices)|, 1/LAMBDA * |ln(1 - i/slices)|]
        let lbound: usize = ((1./LAMBDA) * abs((1.-(i as f32 - 1.)/slices as f32).ln())) as usize;
        let ubound: usize = ((1./LAMBDA) * abs((1.-i as f32/slices as f32).ln())) as usize;

        // if upper bound is above no of buckets (256), reduce to no of buckets
        let ubound: usize = if ubound > BUCKETS { BUCKETS } else { ubound };
        if ubound > lbound && lbound < BUCKETS { bucket_ranges.push(lbound..ubound) };
    }

    debug!("bucket ranges: {:?}", bucket_ranges);

    debug!("kmer_mem: {} B, max_mem: {}B, slices: {}", kmer_mem, max_mem, slices);

    debug!("bucket_ranges: {:?}, len br: {}", bucket_ranges, bucket_ranges.len());
    assert!(bucket_ranges[bucket_ranges.len() - 1].end >= BUCKETS);
    let n_buckets = bucket_ranges.len();

    if bucket_ranges.len() > 1 {
        debug!(
            "{} sequences, {} kmers, {} passes",
            seqs.len(),
            input_kmers,
            bucket_ranges.len()
        );
    }

    debug!("n of seqs: {}", seqs.len());


    let mut all_kmers = Vec::new();
    let mut valid_kmers = Vec::new();
    let mut valid_exts = Vec::new();
    let mut valid_data = Vec::new();

    let mut time_picking = 0.;
    let mut time_summarizing = 0.;
    let mut time_picking_par = 0.;

    if time { println!("time all prepariations before sliced in filter_kmers (s): {}", before_all.elapsed().as_secs_f32()) }

    for (i, bucket_range) in bucket_ranges.into_iter().enumerate() {
        debug!("Processing slice {} of {}", i+1, n_buckets);

        let before_kmer_picking = Instant::now();
        // first step: picking kmers with their exts & data from the reads
        // go through all kmers and sort into bucket according to first four bases
        // all kmers starting with "AAAA" go in kmer_buckets[0], all starting with AAAC go in kmer_buckets[1] and so on
        // when using the first four bases, this needs 256 buckets
        // the buckets are split in to the bucket_ranges to save memory

        // first go trough all kmers to find the length of all buckets (to reserve capacity)
        let mut capacities: [usize; 256] = [0; BUCKETS];

        for &(ref seq, _, _) in seqs {
            // iterate through all kmers in seq
            for kmer in seq.iter_kmers::<K>() {
                // if not stranded choose lexiographically lesser of kmer and rc of kmer
                let min_kmer = if rc_norm {
                    let (min_kmer, _) = kmer.min_rc_flip();
                    min_kmer
                } else {
                    kmer
                };

                // calculate which bucket this kmer belongs to
                let bucket = bucket(min_kmer);
                // check if bucket is in current range and if so, add one to needed capacity
                let in_range = bucket >= bucket_range.start && bucket < bucket_range.end;
                if in_range { capacities[bucket] += 1 }
            }
        }

        debug!("kmer capacities: {:?}, times {}", capacities, mem::size_of::<(K, Exts, D1)>());
        
        let mut kmer_buckets = Vec::new();
        // reserve needed capacity in each bucket
        for capacity in capacities {
            kmer_buckets.push(Vec::with_capacity(capacity));
        }

        // then go through all kmers and add to bucket according to first four bases and current bucket_range
        for &(ref seq, seq_exts, ref d) in seqs {
            // iterate trough all kmers in seq
            for (kmer, exts) in seq.iter_kmer_exts::<K>(seq_exts) {
                // // if not stranded choose lexiographically lesser of kmer and rc of kmer, flip exts if needed
                let (min_kmer, flip_exts) = if rc_norm {
                    let (min_kmer, flip) = kmer.min_rc_flip();
                    let flip_exts = if flip { exts.rc() } else { exts };
                    (min_kmer, flip_exts)
                } else {
                    (kmer, exts)
                };

                // calculate which bucket this kmer belongs to
                let bucket = bucket(min_kmer);

                // check if bucket is in current range and if so, push kmer to bucket
                let in_range = bucket >= bucket_range.start && bucket < bucket_range.end;
                if in_range {
                    kmer_buckets[bucket].push((min_kmer, flip_exts, d.clone()));
                }
            }
        }

        time_picking_par += before_kmer_picking.elapsed().as_secs_f32();

        debug!("size of the slice: {} B", mem::size_of_val(&*kmer_buckets));
        let mut slice_elements: usize = 0;
        for bucket in kmer_buckets.iter() {
            slice_elements += bucket.len();
        }
        debug!("overall elements in this bucket: {slice_elements}");
        debug!("slice size guess (advanced version):
            {} (len) * 24B (ref vec) + {} B (elements size) * {} (elements)
            = {} B", kmer_buckets.len(), mem::size_of::<(K, Exts, D1)>(), slice_elements, kmer_buckets.len() * 24 + (mem::size_of::<(K, Exts, D1)>() * slice_elements));
        
        debug!("no of kmer buckets: {}", kmer_buckets.len());

        time_picking += before_kmer_picking.elapsed().as_secs_f32();

        let before_summarizing = Instant::now();

        if progress {
            println!("Processing bucket {} of {}", i+1, n_buckets);
            for _i in 0..127 {
                print!("-");
            }
            print!("|");
            print!("\n")
        }

        let mut progress_counter = 0;

        for mut kmer_vec in kmer_buckets {
            debug!("bucket {} with {} kmers, capacity of {}", progress_counter, kmer_vec.len(), kmer_vec.capacity());
            progress_counter += 1;
            //debug!("kmers in this bucket: {}", kmer_vec.len());
            if progress & (progress_counter % 2 == 0) { print!("|") };
            kmer_vec.sort_by_key(|elt| elt.0);

            
            // predict amount of unique k-mers found in this bucket
            let size = kmer_vec.iter().group_by(|elt| elt.0).into_iter().count();

            // only works perfectly if min k-mer count is 1, else this might reserve too much 
            // still better than doubling the vector
            // also reserve_exact considers pre-existing free capacity -> this might mostly add to runtime
            valid_kmers.reserve_exact(size);
            valid_exts.reserve_exact(size);
            valid_data.reserve_exact(size);

            // if all k-mers should be reported, also grow all_kmers by exact amount -> this should always be the perfect capacity
            if report_all_kmers {
                all_kmers.reserve_exact(size);
            }

            // group the tuples by the k-mers and iterate over the groups
            for (kmer, kmer_obs_iter) in kmer_vec.into_iter().group_by(|elt| elt.0).into_iter() {
                // summarize group with chosen summarizer and add result to vectors
                let (is_valid, exts, summary_data) = summarizer.summarize(kmer_obs_iter);
                if report_all_kmers {
                    all_kmers.push(kmer);
                }
                if is_valid {
                    valid_kmers.push(kmer);
                    valid_exts.push(exts);
                    valid_data.push(summary_data); 
                }
            }
            debug!("finished bucket {}, current mems: valid_kmers {} Bytes, valid_exts {} Bytes, valid_data {} Bytes", 
                progress_counter, mem::size_of_val(&*valid_kmers), mem::size_of_val(&*valid_exts), mem::size_of_val(&*valid_data))
        }
        if progress { print!("\n") };

        time_summarizing += before_summarizing.elapsed().as_secs_f32();

        debug!("valid kmers - capacity: {}, size: {}, mem: {} Bytes", valid_kmers.capacity(), valid_kmers.len(), mem::size_of_val(&*valid_kmers));
        debug!("valid exts - capacity: {}, size: {}, mem: {} Bytes", valid_exts.capacity(), valid_exts.len(), mem::size_of_val(&*valid_exts));
        debug!("valid data - capacity: {}, size: {}, mem: {} Bytes", valid_data.capacity(), valid_data.len(), mem::size_of_val(&*valid_data));
    }

    if time { 
        println!("time picking par (s): {}", time_picking_par);
        println!("time picking (s): {}", time_picking);
        println!("time summarizing (s): {}", time_summarizing);
    }


    debug!(
        "Unique kmers: {}, All kmers (if returned): {}",
        valid_kmers.len(),
        all_kmers.len(),
    );

    debug!("size of valid kmers: {} Bytes
        size of valid exts: {} Bytes
        size of valid data: {} Bytes", mem::size_of_val(&*valid_kmers), mem::size_of_val(&*valid_exts), mem::size_of_val(&*valid_data));

    let before_hash = Instant::now();
    let hm = BoomHashMap2::new(valid_kmers, valid_exts, valid_data);
    let after_hash = before_hash.elapsed().as_secs_f32();

    let filter_kmers_inner = before_all.elapsed().as_secs_f32();
    if time { 
        println!("time filter_kmers inner (s): {}", filter_kmers_inner);
        println!("time build filter hash map (s): {}", after_hash);
    }
    (
        hm,
        all_kmers,
    )
}

/// Remove extensions in valid_kmers that point to censored kmers. A censored kmer
/// exists in all_kmers but not valid_kmers. Since the kmer exists in this partition,
/// but was censored, we know that we can delete extensions to it.
/// In sharded kmer processing, we will have extensions to kmers in other shards. We don't
/// know whether these are censored until later, so we retain these extension.
pub fn remove_censored_exts_sharded<K: Kmer, D>(
    stranded: bool,
    valid_kmers: &mut [(K, (Exts, D))],
    all_kmers: &[K],
) {
    for idx in 0..valid_kmers.len() {
        let mut new_exts = Exts::empty();
        let kmer = valid_kmers[idx].0;
        let exts = (valid_kmers[idx].1).0;

        for dir in [Dir::Left, Dir::Right].iter() {
            for i in 0..4 {
                if exts.has_ext(*dir, i) {
                    let _ext_kmer = kmer.extend(i, *dir);

                    let ext_kmer = if stranded {
                        _ext_kmer
                    } else {
                        _ext_kmer.min_rc()
                    };

                    let censored = if valid_kmers.binary_search_by_key(&ext_kmer, |d| d.0).is_ok() {
                        // ext_kmer is valid. not censored
                        false
                    } else {
                        // ext_kmer is not valid. if it was in this shard, then we censor it
                        all_kmers.binary_search(&ext_kmer).is_ok()
                    };

                    if !censored {
                        new_exts = new_exts.set(*dir, i);
                    }
                }
            }
        }

        (valid_kmers[idx].1).0 = new_exts;
    }
}

/// Remove extensions in valid_kmers that point to censored kmers. Use this method in a non-partitioned
/// context when valid_kmers includes _all_ kmers that will ultimately be included in the graph.
pub fn remove_censored_exts<K: Kmer, D>(stranded: bool, valid_kmers: &mut [(K, (Exts, D))]) {
    for idx in 0..valid_kmers.len() {
        let mut new_exts = Exts::empty();
        let kmer = valid_kmers[idx].0;
        let exts = (valid_kmers[idx].1).0;

        for dir in [Dir::Left, Dir::Right].iter() {
            for i in 0..4 {
                if exts.has_ext(*dir, i) {
                    let ext_kmer = if stranded {
                        kmer.extend(i, *dir)
                    } else {
                        kmer.extend(i, *dir).min_rc()
                    };

                    let kmer_valid = valid_kmers.binary_search_by_key(&ext_kmer, |d| d.0).is_ok();

                    if kmer_valid {
                        new_exts = new_exts.set(*dir, i);
                    }
                }
            }
        }

        (valid_kmers[idx].1).0 = new_exts;
    }
}
