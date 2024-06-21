// Copyright 2017 10x Genomics

//! Methods for converting sequences into kmers, filtering observed kmers before De Bruijn graph construction, and summarizing 'color' annotations.
use core::slice;
use std::any::Any;
use std::any::TypeId;
use std::collections::btree_map::Range;
use std::fmt::Debug;
use std::marker::PhantomData;
use std::mem;
use std::ops::Deref;
use std::sync::Mutex;

use boomphf::hashmap::BoomHashMap2;
use itertools::Itertools;
use log::debug;
use serde_json::de;
use rayon::prelude::*;

use crate::Dir;
use crate::Exts;
use crate::Kmer;
use crate::Vmer;

fn bucket<K: Kmer>(kmer: K) -> usize {
    (kmer.get(0) as usize) << 6
        | (kmer.get(1) as usize) << 4
        | (kmer.get(2) as usize) << 2
        | (kmer.get(3) as usize)
}

/// Implement this trait to control how multiple observations of a kmer
/// are carried forward into a DeBruijn graph.
pub trait KmerSummarizer<DI, DO> {
    /// The input `items` is an iterator over kmer observations. Input observation
    /// is a tuple of (kmer, extensions, data). The summarize function inspects the
    /// data and returns a tuple indicating:
    /// * whether this kmer passes the filtering criteria (e.g. is there a sufficient number of observation)
    /// * the accumulated Exts of the kmer
    /// * a summary data object of type `DO` that will be used as a color annotation in the DeBruijn graph.
    fn summarize<K: Kmer, F: Iterator<Item = (K, Exts, DI)>>(&self, items: F) -> (bool, Exts, DO, usize);
}

/// A simple KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times. The metadata returned about a Kmer
/// is the number of times it was observed, capped at 2^16.
pub struct CountFilter {
    min_kmer_obs: usize,
}

impl CountFilter {
    /// Construct a `CountFilter` KmerSummarizer only accepts kmers that are observed
    /// at least `min_kmer_obs` times.
    pub fn new(min_kmer_obs: usize) -> CountFilter {
        CountFilter { min_kmer_obs }
    }
}

impl<D> KmerSummarizer<D, u16> for CountFilter {
    fn summarize<K, F: Iterator<Item = (K, Exts, D)>>(&self, items: F) -> (bool, Exts, u16, usize) {
        let mut all_exts = Exts::empty();
        let mut count = 0u16;
        for (_, exts, _) in items {
            count = count.saturating_add(1);
            all_exts = all_exts.add(exts);
        }

        (count as usize >= self.min_kmer_obs, all_exts, count, 0)
    }
}

/// A simple KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times. The metadata returned about a Kmer
/// is a vector of the unique data values observed for that kmer.
pub struct CountFilterSet<D> {
    min_kmer_obs: usize,
    phantom: PhantomData<D>,
}

impl<D> CountFilterSet<D> {
    /// Construct a `CountFilterSet` KmerSummarizer only accepts kmers that are observed
    /// at least `min_kmer_obs` times.
    pub fn new(min_kmer_obs: usize) -> CountFilterSet<D> {
        CountFilterSet {
            min_kmer_obs,
            phantom: PhantomData,
        }
    }
}

impl<D: Ord + Debug> KmerSummarizer<D, Vec<D>> for CountFilterSet<D> {
    fn summarize<K, F: Iterator<Item = (K, Exts, D)>>(&self, items: F) -> (bool, Exts, Vec<D>, usize) {
        let mut all_exts = Exts::empty();

        let mut out_data: Vec<D> = Vec::with_capacity(items.size_hint().0);

        let mut nobs = 0i32;
        for (_, exts, d) in items {
            out_data.push(d);
            all_exts = all_exts.add(exts);
            nobs += 1;
        }

        let max = out_data.len();

        out_data.sort();
        out_data.dedup();

        
        //debug!("count filter set out data len: {}", out_data.len());
        (nobs as usize >= self.min_kmer_obs, all_exts, out_data, max)
        
    }
}

/// A simple KmerSummarizer that only accepts kmers that are observed
/// at least a given number of times. The metadata returned about a Kmer
/// is a vector of the unique data values observed for that kmer.
pub struct CountFilterComb<D> {
    min_kmer_obs: usize,
    pub sum_datas: u64,
    phantom: PhantomData<D>,
}

impl<D> CountFilterComb<D> {
    /// Construct a `CountFilterSet` KmerSummarizer only accepts kmers that are observed
    /// at least `min_kmer_obs` times.
    /// data is ([vec with tags], count)
    pub fn new(min_kmer_obs: usize) -> CountFilterComb<D> {
        CountFilterComb {
            min_kmer_obs,
            sum_datas: 0,
            phantom: PhantomData,
        }
    }

}

impl<D: Ord + Debug> KmerSummarizer<D, (Vec<D>, i32)> for CountFilterComb<D> {
    fn summarize<K: Kmer, F: Iterator<Item = (K, Exts, D)>>(&self, items: F) -> (bool, Exts, (Vec<D>, i32), usize) {
        let mut all_exts = Exts::empty();

        let mut out_data: Vec<D> = Vec::with_capacity(items.size_hint().0);
        if out_data.len() > 9999 {
            println!("od size hint: {:?}", items.size_hint());
        }
        let mut kmer: K = Kmer::empty();

        let mut nobs = 0i32;
        for (k, exts, d) in items {
            kmer = k;
            out_data.push(d); // uses a shit ton of heap memory
            all_exts = all_exts.add(exts);
            nobs += 1;
        }
        

        if out_data.len() > 9999 {debug!(
            "odl {:?}
            kmer: {:?}
            size: {:?}B", out_data.len(), kmer, mem::size_of_val(&out_data)

        )}

        let max = out_data.len();

        out_data.sort();
        out_data.dedup();

        /* debug!("out_data post sorting: {:?}", out_data);
        debug!("nobs: {}", nobs);
        debug!("result: {:?}", (nobs as usize >= self.min_kmer_obs, all_exts, &out_data)); */
        //debug!("count filter set out data len: {}", out_data.len());
        (nobs as usize >= self.min_kmer_obs, all_exts, (out_data, nobs), max)
        
    }
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
pub fn filter_kmers_parallel<K: Kmer + Sync + Send, V: Vmer + Sync, D1: Clone + Debug + Sync, DS: Clone + Sync + Send, S: KmerSummarizer<D1, DS> +  Send>(
    seqs: &[(V, Exts, D1)],
    //summarizer: &dyn Deref<Target = S>,
    // summarizer without wrapper, why wrapper???
    summarizer: Box<S>,
    stranded: bool,
    report_all_kmers: bool,
    memory_size: usize,
) -> (BoomHashMap2<K, Exts, DS>, Vec<K>)
where
    DS: Debug,
{
    let rc_norm = !stranded;

    let shared_summarizer = Mutex::new(summarizer);

    //println!("sequence input: {:?}", seqs);

    // Estimate memory consumed by Kmer vectors, and set iteration count appropriately
    let input_kmers: usize = seqs
        .iter()
        .map(|&(ref vmer, _, _)| vmer.len().saturating_sub(K::k() - 1))
        .sum();
    let kmer_mem = input_kmers * mem::size_of::<(K, D1)>();
    let max_mem = memory_size * 10_usize.pow(9);
    let slices_seq = kmer_mem / max_mem + 1;
    let slices = slices_seq * rayon::current_num_threads();
    let sz = 256 / slices + 1;

    debug!("kmer_mem: {}, max_mem: {}, slices: {}, sz: {}", kmer_mem, max_mem, slices, sz);

    let mut bucket_ranges = Vec::with_capacity(slices);
    let mut start = 0;
    while start < 256 {
        bucket_ranges.push(start..start + sz);
        start += sz;
    }
    debug!("start: {}, bucket_ranges: {:?}, len br: {}", start, bucket_ranges, bucket_ranges.len());
    assert!(bucket_ranges[bucket_ranges.len() - 1].end >= 256);
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

    let mut shared_data: Mutex<Vec<Vec<(Vec<K>, Vec<K>, Vec<Exts>, Vec<DS>)>>> = Mutex::new(vec![vec![]; n_buckets]);
    debug!("data_out empty: {:?}", shared_data.lock());

    /* let mut all_kmers = Vec::new();
    let mut valid_kmers = Vec::new();
    let mut valid_exts = Vec::new();
    let mut valid_data = Vec::new(); */



    // parallel start

    bucket_ranges.into_par_iter().enumerate().for_each(&|(i, bucket_range): (usize, std::ops::Range<usize>)| {

        debug!("Processing bucket {} of {}", i, n_buckets);

        let mut all_kmers = Vec::new();
        let mut valid_kmers = Vec::new();
        let mut valid_exts = Vec::new();
        let mut valid_data = Vec::new();

        let mut kmer_buckets = vec![Vec::new(); 256];

        for &(ref seq, seq_exts, ref d) in seqs {
            for (kmer, exts) in seq.iter_kmer_exts::<K>(seq_exts) {
                //println!("kmer: {:?}, exts: {:?}", kmer, exts);
                let (min_kmer, flip_exts) = if rc_norm {
                    let (min_kmer, flip) = kmer.min_rc_flip();
                    let flip_exts = if flip { exts.rc() } else { exts };
                    (min_kmer, flip_exts)
                } else {
                    (kmer, exts)
                };
                let bucket = bucket(min_kmer);

                if bucket >= bucket_range.start && bucket < bucket_range.end {
                    kmer_buckets[bucket].push((min_kmer, flip_exts, d.clone()));
                }
            }
        }
        
        debug!("no of kmer buckets: {}", kmer_buckets.len());

        for mut kmer_vec in kmer_buckets {
            debug!("kmers in this bucket: {}", kmer_vec.len());
            //println!("kmers in this bucket: {:?}", kmer_vec);
            kmer_vec.sort_by_key(|elt| elt.0);

            for (kmer, kmer_obs_iter) in kmer_vec.into_iter().group_by(|elt| elt.0).into_iter() {
                let summarizer = shared_summarizer.lock().expect("unlock shared filter summarizer");
                let (is_valid, exts, summary_data, _) = summarizer.summarize(kmer_obs_iter);
                if report_all_kmers {
                    all_kmers.push(kmer);
                }
                if is_valid {
                    valid_kmers.push(kmer);
                    valid_exts.push(exts);
                    valid_data.push(summary_data);
                }
            }
        }

        debug!("processed bucket {i}");

        //println!("all k: {:?}\n v k: {:?}\n v e: {:?}\n v d: {:?}", all_kmers, valid_kmers, valid_exts, valid_data);

        let mut data_out = shared_data.lock().expect("unlock shared filter data");

        data_out[i].push((all_kmers, valid_kmers, valid_exts, valid_data));
    });
    // parallel end

    /* debug!(
        "Unique kmers: {}, All kmers (if returned): {}",
        valid_kmers.len(),
        all_kmers.len(),
    ); */

    let data_out = shared_data.lock().expect("final unlock shared filter data");

    let mut all_kmers = Vec::new();
    let mut valid_kmers = Vec::new();
    let mut valid_exts = Vec::new();
    let mut valid_data = Vec::new();

    debug!("data_out: {:?}", data_out);

    for bucket in data_out.iter() {
        all_kmers.append(&mut bucket[0].0.clone());
        valid_kmers.append(&mut bucket[0].1.clone());
        valid_exts.append(&mut bucket[0].2.clone());
        valid_data.append(&mut bucket[0].3.clone());
    } 
    
    debug!("data_out2: {:?}", data_out);

    (
        BoomHashMap2::new(valid_kmers, valid_exts, valid_data),
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
pub fn filter_kmers<K: Kmer, V: Vmer, D1: Clone + Debug, DS, S: KmerSummarizer<D1, DS>>(
    seqs: &[(V, Exts, D1)],
    summarizer: &dyn Deref<Target = S>,
    stranded: bool,
    report_all_kmers: bool,
    memory_size: usize,
) -> (BoomHashMap2<K, Exts, DS>, Vec<K>)
where
    DS: Debug,
{
    let rc_norm = !stranded;

    //println!("sequence input: {:?}", seqs);

    // Estimate memory consumed by Kmer vectors, and set iteration count appropriately
    let input_kmers: usize = seqs
        .iter()
        .map(|&(ref vmer, _, _)| vmer.len().saturating_sub(K::k() - 1))
        .sum();
    let kmer_mem = input_kmers * mem::size_of::<(K, D1)>();
    debug!("size used for calculation: {}B", mem::size_of::<(K, D1)>());
    debug!("size of kmer, E, D: {}B", mem::size_of::<(K, Exts, D1)>());
    debug!("size of K: {}B, size of Exts: {}B, size of D1: {}", mem::size_of::<K>(), mem::size_of::<Exts>(), mem::size_of::<D1>());
    debug!("type D1: {}", std::any::type_name::<D1>());

    let max_mem = memory_size * 10_usize.pow(9);
    let slices = kmer_mem / max_mem + 1;
    //let slices = 257;
    let sz = 256 / slices + 1;

    debug!("kmer_mem: {}B, max_mem: {}B, slices: {}, sz: {}", kmer_mem, max_mem, slices, sz);

    let mut bucket_ranges = Vec::with_capacity(slices);
    let mut start = 0;
    while start < 256 {
        bucket_ranges.push(start..start + sz);
        start += sz;
    }
    debug!("start: {}, bucket_ranges: {:?}, len br: {}", start, bucket_ranges, bucket_ranges.len());
    assert!(bucket_ranges[bucket_ranges.len() - 1].end >= 256);
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

    let mut data_lengths: usize = 0;

    let mut all_kmers = Vec::new();
    let mut valid_kmers = Vec::new();
    let mut valid_exts = Vec::new();
    let mut valid_data = Vec::new();
    for (i, bucket_range) in bucket_ranges.into_iter().enumerate() {
        debug!("Processing bucket {} of {}", i, n_buckets);

        let mut kmer_buckets = vec![Vec::new(); 256];

        for &(ref seq, seq_exts, ref d) in seqs {
            for (kmer, exts) in seq.iter_kmer_exts::<K>(seq_exts) {
                //println!("kmer: {:?}, exts: {:?}", kmer, exts);
                let (min_kmer, flip_exts) = if rc_norm {
                    let (min_kmer, flip) = kmer.min_rc_flip();
                    let flip_exts = if flip { exts.rc() } else { exts };
                    (min_kmer, flip_exts)
                } else {
                    (kmer, exts)
                };
                let bucket = bucket(min_kmer);

                if bucket >= bucket_range.start && bucket < bucket_range.end {
                    kmer_buckets[bucket].push((min_kmer, flip_exts, d.clone())); // also lots of heap memory
                }
            }
        }
        debug!("size of the bucket: {}B
            len of kmer bucket: {}", mem::size_of_val(&*kmer_buckets), kmer_buckets.len());
        
        debug!("no of kmer buckets: {}", kmer_buckets.len());

        for mut kmer_vec in kmer_buckets {
            debug!("kmers in this bucket: {}", kmer_vec.len());
            //println!("kmers in this bucket: {:?}", kmer_vec);
            kmer_vec.sort_by_key(|elt| elt.0);

            for (kmer, kmer_obs_iter) in kmer_vec.into_iter().group_by(|elt| elt.0).into_iter() {
                let (is_valid, exts, summary_data, max) = summarizer.summarize(kmer_obs_iter);
                data_lengths += max;
                if report_all_kmers {
                    all_kmers.push(kmer);
                }
                if is_valid {
                    valid_kmers.push(kmer);
                    valid_exts.push(exts);
                    valid_data.push(summary_data); // also lots of memory // but below limit/in the end amount around limit
                }
            }
        }
        debug!("sum data lengths after this bucket: {}", data_lengths);
    }

    debug!(
        "Unique kmers: {}, All kmers (if returned): {}",
        valid_kmers.len(),
        all_kmers.len(),
    );
    debug!("size of valid kmers: {}B
        size of valid exts: {}B
        size of valid data: {}B", mem::size_of_val(&*valid_kmers), mem::size_of_val(&*valid_exts), mem::size_of_val(&*valid_data));
    (
        BoomHashMap2::new(valid_kmers, valid_exts, valid_data),
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
