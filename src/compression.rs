// Copyright 2017 10x Genomics

//! Create compressed DeBruijn graphs from uncompressed DeBruijn graphs, or a collection of disjoint DeBruijn graphs.
use bit_set::BitSet;
use indicatif::{ProgressBar, ProgressIterator, ProgressStyle};
use log::debug;
use rayon::current_num_threads;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::collections::{HashMap, VecDeque};
use std::fmt::Debug;
use std::marker::PhantomData;
use std::mem;
use std::sync::{Arc, Mutex};
use std::time::Instant;

use crate::dna_string::{DnaString, PackedDnaStringSet};
use crate::graph::{BaseGraph, DebruijnGraph};
use crate::Dir;
use crate::Exts;
use crate::Kmer;
use crate::Vmer;
use boomphf::hashmap::BoomHashMap2;

#[derive(Copy, Clone)]
enum ExtMode<K: Kmer> {
    Unique(K, Dir, Exts),
    Terminal(Exts),
}

#[derive(Copy, Clone)]
enum ExtModeNode {
    Unique(usize, Dir, Exts),
    Terminal(Exts),
}

/// Customize the path-compression process. Implementing this trait lets the user
/// control how the per-kmer data (of type `D`) is summarized into a per-path
/// summary data (of type `DS`). It also let's the user introduce new breaks into
/// paths by inspecting in the per-kmer data of a proposed with `join_test_kmer`
/// function.
pub trait CompressionSpec<D> {
    /// combine the data of two nodes
    fn reduce(&self, path_object: D, kmer_object: &D) -> D;
    /// check if the data of two nodes can be combined
    fn join_test(&self, d1: &D, d2: &D) -> bool;
}

/// Simple implementation of `CompressionSpec` that lets you provide that data reduction function as a closure
pub struct SimpleCompress<D, F> {
    func: F,
    d: PhantomData<D>,
}

impl<D, F> SimpleCompress<D, F> {
    pub fn new(func: F) -> SimpleCompress<D, F> {
        SimpleCompress {
            func,
            d: PhantomData,
        }
    }
}

impl<D, F> CompressionSpec<D> for SimpleCompress<D, F>
where
    for<'r> F: Fn(D, &'r D) -> D,
{
    fn reduce(&self, d: D, other: &D) -> D {
        (self.func)(d, other)
    }

    fn join_test(&self, _: &D, _: &D) -> bool {
        true
    }
}

/// Extending trait CompressionSpec for compression
pub struct ScmapCompress<D> {
    d: PhantomData<D>,
}

impl<D> ScmapCompress<D> {
    pub fn new() -> ScmapCompress<D> {
        ScmapCompress { d: PhantomData }
    }
}

impl<D> Default for ScmapCompress<D> {
    fn default() -> Self {
        Self::new()
    }
}

impl<D: PartialEq> CompressionSpec<D> for ScmapCompress<D>
where
    D: Debug,
{
    fn reduce(&self, d: D, other: &D) -> D {
        if d != *other {
            panic!("{:?} != {:?}, Should not happen", d, *other);
        }
        d
    }

    fn join_test(&self, d1: &D, d2: &D) -> bool {
        d1 == d2
    }
}

/// CompressionSpec with check and function
pub struct CheckCompress<D, F1, F2> {
    reduce_func: F1,
    join_func: F2,
    d: PhantomData<D>
} 

impl<D, F1, F2> CheckCompress<D, F1, F2> {
    pub fn new(reduce_func: F1, join_func: F2) -> Self {
        CheckCompress {
            reduce_func,
            join_func,
            d: PhantomData,
        }
    }
}

impl<D, F1, F2> CompressionSpec<D> for CheckCompress<D, F1, F2>
where
    for<'r> F1: Fn(D, &'r D) -> D,
    for<'r> F2: Fn(&'r D, &'r D) -> bool
{
    fn reduce(&self, d: D, other: &D) -> D {
        (self.reduce_func)(d, other)
    }

    fn join_test(&self, d: &D, other: &D) -> bool {
        (self.join_func)(d, other)
    }
}


struct CompressFromGraph<'a, 'b, K: 'a + Kmer, D: 'a + PartialEq, S: CompressionSpec<D>> {
    stranded: bool,
    d: PhantomData<D>,
    spec: &'b S,
    available_nodes: BitSet,
    graph: &'a DebruijnGraph<K, D>,
}

impl<K, D, S> CompressFromGraph<'_, '_, K, D, S>
where
    K: Kmer + Send + Sync,
    D: Debug + Clone + PartialEq,
    S: CompressionSpec<D>,
{
    #[inline(never)]
    fn try_extend_node(&mut self, node: usize, dir: Dir) -> ExtModeNode {
        let node = self.graph.get_node(node);
        let bases = node.sequence();
        let exts = node.exts();

        if exts.num_ext_dir(dir) != 1
            || (!self.stranded && node.len() == K::k() && bases.get_kmer::<K>(0).is_palindrome())
        {
            ExtModeNode::Terminal(exts.single_dir(dir))
        } else {
            // Get the next kmer
            let ext_base = exts.get_unique_extension(dir).expect("should be unique");
            let end_kmer: K = bases.term_kmer(dir);

            let next_kmer = end_kmer.extend(ext_base, dir);
            let (next_node_id, next_side_incoming, rc) = match self.graph.find_link(next_kmer, dir)
            {
                Some(e) => e,
                None => {
                    println!("dir: {:?}, Lmer: {:?}, exts: {:?}", dir, bases, exts);
                    println!("end kmer: {:?}", end_kmer);
                    println!("No kmer: {:?}", next_kmer);
                    println!("rc: {:?}", next_kmer.min_rc());
                    panic!("No kmer: {:?}", next_kmer)
                }
            };

            let next_node = self.graph.get_node(next_node_id);
            let next_exts = next_node.exts();

            let consistent = (next_node.len() == K::k())
                || match (dir, next_side_incoming, rc) {
                    (Dir::Left, Dir::Right, false) => true,
                    (Dir::Left, Dir::Left, true) => true,
                    (Dir::Right, Dir::Left, false) => true,
                    (Dir::Right, Dir::Right, true) => true,
                    _ => {
                        println!("dir: {:?}, Lmer: {:?}, exts: {:?}", dir, bases, exts);
                        println!("end kmer: {:?}", end_kmer);
                        println!("next kmer: {:?}", next_kmer);
                        println!("rc: {:?}", next_kmer.min_rc());
                        println!(
                            "next bases: {:?}, next_side_incoming: {:?}, rc: {:?}",
                            next_node.sequence(),
                            next_side_incoming,
                            rc
                        );
                        false
                    }
                };
            assert!(consistent);

            // We can include this kmer in the line if:
            // a) it exists in the partition, and is still unused
            // b) the kmer we go to has a unique extension back in our direction
            // c) the new edge is not of length K and a palindrome
            // d) the color of the current and next node is same

            if !self.available_nodes.contains(next_node_id)
                || (!self.stranded && next_kmer.is_palindrome())
                || !self.spec.join_test(node.data(), next_node.data())
            {
                // Next kmer isn't in this partition,
                // or we've already used it,
                // or it's palindrom and we are not stranded
                // or the colors were not same
                return ExtModeNode::Terminal(exts.single_dir(dir));
            }

            // orientation of next edge
            let next_side_outgoing = next_side_incoming.flip();

            let incoming_count = next_exts.num_ext_dir(next_side_incoming);
            let outgoing_exts = next_exts.single_dir(next_side_outgoing);

            if incoming_count == 0 {
                println!("dir: {:?}, Lmer: {:?}, exts: {:?}", dir, bases, exts);
                println!("end kmer: {:?}", end_kmer);
                println!("next_node: {:?}", next_node);
                println!("next_node data: {:?}", next_node.sequence());
                panic!("unreachable");
            } else if incoming_count == 1 {
                // We have a unique path to next_kmer -- include it
                ExtModeNode::Unique(next_node_id, next_side_outgoing, outgoing_exts)
            } else {
                // there's more than one path
                // into the target kmer - don't include it
                ExtModeNode::Terminal(exts.single_dir(dir))
            }
        }
    }

    /// Generate complete unbranched edges
    fn extend_node(&mut self, start_node: usize, start_dir: Dir) -> (Vec<(usize, Dir)>, Exts) {
        let mut current_dir = start_dir;
        let mut current_node = start_node;
        let mut path = Vec::new();
        let final_exts: Exts; // must get set below

        self.available_nodes.remove(start_node);

        loop {
            let ext_result = self.try_extend_node(current_node, current_dir);

            match ext_result {
                ExtModeNode::Unique(next_node, next_dir_outgoing, _) => {
                    let next_dir_incoming = next_dir_outgoing.flip();
                    path.push((next_node, next_dir_incoming));
                    self.available_nodes.remove(next_node);
                    current_node = next_node;
                    current_dir = next_dir_outgoing;
                }
                ExtModeNode::Terminal(ext) => {
                    final_exts = ext;
                    break;
                }
            }
        }

        (path, final_exts)
    }

    // Determine the sequence and extensions of the maximal unbranched
    // edge, centered around the given edge number
    #[inline(never)]
    fn build_node(&mut self, seed_node: usize) -> (DnaString, Exts, VecDeque<(usize, Dir)>, D) {
        let (l_path, l_ext) = self.extend_node(seed_node, Dir::Left);
        let (r_path, r_ext) = self.extend_node(seed_node, Dir::Right);

        // Stick together edge chunks to get full edge sequence
        let mut node_path = VecDeque::new();

        let mut node_data: D = self.graph.get_node(seed_node).data().clone();
        node_path.push_back((seed_node, Dir::Left));

        // Add on the left path
        for &(next_node, incoming_dir) in l_path.iter() {
            node_path.push_front((next_node, incoming_dir.flip()));
            node_data = self
                .spec
                .reduce(node_data, self.graph.get_node(next_node).data());
        }

        // Add on the right path
        for &(next_node, incoming_dir) in r_path.iter() {
            node_path.push_back((next_node, incoming_dir));
            node_data = self
                .spec
                .reduce(node_data, self.graph.get_node(next_node).data());
        }

        let left_extend = match l_path.last() {
            None => l_ext,
            Some(&(_, Dir::Left)) => l_ext.complement(),
            Some(&(_, Dir::Right)) => l_ext,
        };

        let right_extend = match r_path.last() {
            None => r_ext,
            Some(&(_, Dir::Left)) => r_ext,
            Some(&(_, Dir::Right)) => r_ext.complement(),
        };

        let path_seq = self.graph.sequence_of_path(node_path.iter());

        // return sequence and extensions
        (
            path_seq,
            Exts::from_single_dirs(left_extend, right_extend),
            node_path,
            node_data,
        )
    }

    /// Simplify a compressed Debruijn graph by merging adjacent unbranched nodes, and optionally
    /// censoring some nodes
    fn compress_graph(
        stranded: bool,
        compression: &S,
        mut old_graph: DebruijnGraph<K, D>,
        censor_nodes: Option<Vec<usize>>,
    ) -> DebruijnGraph<K, D> {
        let n_nodes = old_graph.len();
        let mut available_nodes = BitSet::with_capacity(n_nodes);
        for i in 0..n_nodes {
            available_nodes.insert(i);
        }

        if let Some(c) = censor_nodes {
            for censor in c {
                available_nodes.remove(censor);
            }
        }

        old_graph.fix_exts(Some(&available_nodes));

        let mut comp = CompressFromGraph {
            spec: compression,
            stranded,
            graph: &old_graph,
            available_nodes,
            d: PhantomData,
        };

        // FIXME -- clarify requirements around state of extensions
        let mut graph = BaseGraph::new(stranded);

        for node_counter in 0..n_nodes {
            if comp.available_nodes.contains(node_counter) {
                let (seq, exts, _, data) = comp.build_node(node_counter);
                graph.add(&seq, exts, data);
            }
        }

        // We will have some hanging exts due to
        let mut dbg = graph.finish();
        dbg.fix_exts(None); // ????
        debug_assert!(dbg.is_compressed(compression).is_none());
        dbg
    }
}

/// Perform path-compression on a (possibly partially compressed) DeBruijn graph
pub fn compress_graph<
    K: Kmer + Send + Sync,
    D: Clone + Debug + PartialEq,
    S: CompressionSpec<D>,
>(
    stranded: bool,
    spec: &S,
    old_graph: DebruijnGraph<K, D>,
    censor_nodes: Option<Vec<usize>>,
) -> DebruijnGraph<K, D> {
    CompressFromGraph::<K, D, S>::compress_graph(stranded, spec, old_graph, censor_nodes)
}

//////////////////////////////
// Compress from Hash a new Struct
//////////////////////////////
/// Generate a compressed DeBruijn graph from hash_index
struct CompressFromHash<'a, 'b, K: 'a + Kmer, D: 'a, S: CompressionSpec<D>> {
    stranded: bool,
    k: PhantomData<K>,
    d: PhantomData<D>,
    spec: &'b S,
    available_kmers: BitSet,
    index: &'a BoomHashMap2<K, Exts, D>,
}

/// Compression of paths in Debruijn graph
impl<K: Kmer +  Send + Sync, D: Clone + Debug + Send + Sync, S: CompressionSpec<D> + Sync> CompressFromHash<'_, '_, K, D, S> {
    fn get_kmer_data(&self, kmer: &K) -> (&Exts, &D) {
        match self.index.get(kmer) {
            Some(data) => data,
            None => panic!("couldn't find kmer {:?}", kmer),
        }
    }

    fn get_kmer_id(&self, kmer: &K) -> Option<usize> {
        self.index.get_key_id(kmer)
    }

    /// Attempt to extend kmer v in direction dir. Return:
    ///  - Unique(nextKmer, nextDir) if a single unique extension
    ///    is possible.  nextDir indicates the direction to extend nextMker
    ///    to preserve the direction of the extension.
    /// - Term(ext) no unique extension possible, indicating the extensions at this end of the line
    fn try_extend_kmer(&self, kmer: K, dir: Dir) -> ExtMode<K> {
        // metadata of start kmer
        let (exts, ref kmer_data) = self.get_kmer_data(&kmer);

        // kmer is marked terminal if it has not one extension in one direction (if clear path always 1) 
        // or if the graph is not stranded and the kmer is a palindrome
        if exts.num_ext_dir(dir) != 1 || (!self.stranded && kmer.is_palindrome()) {
            ExtMode::Terminal(exts.single_dir(dir))
        } else {
            // Get the next kmer
            let ext_base = exts.get_unique_extension(dir).expect("should be unique");

            let mut next_kmer = kmer.extend(ext_base, dir);

            let mut do_flip = false;
            
            // decide if direction needs to be changed (how?????????) turn kmer into rc
            if !self.stranded {
                let flip_rc = next_kmer.min_rc_flip();
                do_flip = flip_rc.1;
                next_kmer = flip_rc.0;
            }

            let next_dir = dir.cond_flip(do_flip);
            let is_palindrome = !self.stranded && next_kmer.is_palindrome();

            // We can include this kmer in the line if:
            // a) it exists in the partition and is still unused
            // b) the kmer we go to has a unique extension back in our direction

            // Check condition a)
            match self.get_kmer_id(&next_kmer) {
                Some(id) if self.available_kmers.contains(id) => (),

                // This kmer isn't in this partition, or we've already used it
                _ => return ExtMode::Terminal(exts.single_dir(dir)),
            }

            // Check condition b)
            // Direction we're approaching the new kmer from
            let new_incoming_dir = dir.flip().cond_flip(do_flip);
            let next_kmer_r = self.get_kmer_data(&next_kmer);
            let (next_kmer_exts, ref next_kmer_data) = next_kmer_r;
            let incoming_count = next_kmer_exts.num_ext_dir(new_incoming_dir);
            let outgoing_exts = next_kmer_exts.single_dir(new_incoming_dir.flip());

            // Test if the spec let's us combine these into the same path
            let can_join = self.spec.join_test(kmer_data, next_kmer_data);

            if incoming_count == 0 && !is_palindrome {
                println!("{:?}, {:?}, {:?}", kmer, exts, kmer_data);
                println!(
                    "{:?}, {:?}, {:?}",
                    next_kmer, next_kmer_exts, next_kmer_data
                );
                panic!("unreachable");
            } else if can_join && incoming_count == 1 && !is_palindrome {
                // We have a unique path to next_kmer -- include it
                ExtMode::Unique(next_kmer, next_dir, outgoing_exts)
            } else {
                // there's more than one path
                // into the target kmer - don't include it
                ExtMode::Terminal(exts.single_dir(dir))
            }
        }
    }

    /// Attempt to extend kmer v in direction dir. Return:
    ///  - Unique(nextKmer, nextDir) if a single unique extension
    ///    is possible.  nextDir indicates the direction to extend nextMker
    ///    to preserve the direction of the extension.
    /// - Term(ext) no unique extension possible, indicating the extensions at this end of the line
    fn try_extend_kmer_par(&self, kmer: K, dir: Dir, path: &mut [(K, Dir)]) -> ExtMode<K> {
        // metadata of start kmer
        let (exts, ref kmer_data) = self.get_kmer_data(&kmer);

        // kmer is marked terminal if it has not one extension in one direction (if clear path always 1) 
        // or if the graph is not stranded and the kmer is a palindrome
        if exts.num_ext_dir(dir) != 1 || (!self.stranded && kmer.is_palindrome()) {
            ExtMode::Terminal(exts.single_dir(dir))
        } else {
            // Get the next kmer
            let ext_base = exts.get_unique_extension(dir).expect("should be unique");

            let mut next_kmer = kmer.extend(ext_base, dir);

            let mut do_flip = false;
            
            // decide if direction needs to be changed (how?????????) turn kmer into rc
            if !self.stranded {
                let flip_rc = next_kmer.min_rc_flip();
                do_flip = flip_rc.1;
                next_kmer = flip_rc.0;
            }

            let next_dir = dir.cond_flip(do_flip);
            let is_palindrome = !self.stranded && next_kmer.is_palindrome();

            // We can include this kmer in the line if:
            // a) it exists in the partition and is still unused
            // b) the kmer we go to has a unique extension back in our direction

            // condition a) is skipped cause of parallel process
            // instead check if already in path
            // might be time consuming, can try to replace with hashmap later
            // to avoid circluar graphs leading infinite loops
            // also check if just checking for next_dir would suffice

            if path.contains(&(next_kmer, Dir::Left)) {
                // The next kmer is already in the path, break loop
                return ExtMode::Terminal(exts.single_dir(dir))
            }

            if path.contains(&(next_kmer, Dir::Right)) {
                // The next kmer is already in the path, break loop
                return ExtMode::Terminal(exts.single_dir(dir))
            }

            // Check condition b)
            // Direction we're approaching the new kmer from
            let new_incoming_dir = dir.flip().cond_flip(do_flip);
            let next_kmer_r = self.get_kmer_data(&next_kmer);
            let (next_kmer_exts, ref next_kmer_data) = next_kmer_r;
            let incoming_count = next_kmer_exts.num_ext_dir(new_incoming_dir);
            let outgoing_exts = next_kmer_exts.single_dir(new_incoming_dir.flip());

            // Test if the spec let's us combine these into the same path
            let can_join = self.spec.join_test(kmer_data, next_kmer_data);

            if incoming_count == 0 && !is_palindrome {
                println!("{:?}, {:?}, {:?}", kmer, exts, kmer_data);
                println!(
                    "{:?}, {:?}, {:?}",
                    next_kmer, next_kmer_exts, next_kmer_data
                );
                panic!("unreachable");
            } else if can_join && incoming_count == 1 && !is_palindrome {
                // We have a unique path to next_kmer -- include it
                ExtMode::Unique(next_kmer, next_dir, outgoing_exts)
            } else {
                // there's more than one path
                // into the target kmer - don't include it
                ExtMode::Terminal(exts.single_dir(dir))
            }
        }
    }


    /// Build the maximal line starting at kmer in direction dir, at most max_dist long.
    /// Also return the extensions at the end of this line.
    /// Sub-lines break if their extensions are not available in this shard
    #[inline(never)]
    fn extend_kmer(&mut self, kmer: K, start_dir: Dir, path: &mut Vec<(K, Dir)>) -> Exts {
        let mut current_dir = start_dir;
        let mut current_kmer = kmer;
        path.clear();

        let final_exts: Exts; // must get set below

        // get id of kmer and remove from available kmers
        let id = self.get_kmer_id(&kmer).expect("should have this kmer");
        let _ = self.available_kmers.remove(id);

        loop {
            let ext_result = self.try_extend_kmer(current_kmer, current_dir);

            match ext_result {
                ExtMode::Unique(next_kmer, next_dir, _) => {
                    path.push((next_kmer, next_dir));
                    let next_id = self.get_kmer_id(&next_kmer).expect("should have this kmer");
                    self.available_kmers.remove(next_id);
                    current_kmer = next_kmer;
                    current_dir = next_dir;
                }
                ExtMode::Terminal(ext) => {
                    final_exts = ext;
                    break;
                }
            }
        }

        final_exts
    }

    fn extend_kmer_par(&mut self, kmer: K, start_dir: Dir, path: &mut Vec<(K, Dir)>) -> Exts {
        let mut current_dir = start_dir;
        let mut current_kmer = kmer;
        path.clear();

        let final_exts: Exts; // must get set below

        loop {
            let ext_result = self.try_extend_kmer_par(current_kmer, current_dir, path);

            match ext_result {
                ExtMode::Unique(next_kmer, next_dir, _) => {
                    path.push((next_kmer, next_dir));
                    //let next_id = self.get_kmer_id(&next_kmer).expect("should have this kmer");
                    //self.available_kmers.remove(next_id);
                    current_kmer = next_kmer;
                    current_dir = next_dir;
                }
                ExtMode::Terminal(ext) => {
                    final_exts = ext;
                    break;
                }
            }
        }

        final_exts
    }

    /// Build the edge surrounding a kmer
    #[inline(never)]
    fn  build_node(
        &mut self,
        seed_id: usize,
        path: &mut Vec<(K, Dir)>,
        edge_seq: &mut VecDeque<u8>,
    ) -> (Exts, D) {
        let seed: K = *self.index.get_key(seed_id).expect("Index out of bound");
        edge_seq.clear();
        for i in 0..K::k() {
            edge_seq.push_back(seed.get(i));
        }

        let mut node_data = self.get_kmer_data(&seed).1.clone();

        // Unique path from seed kmer with Dir Left is built
        let l_ext = self.extend_kmer(seed, Dir::Left, path);


        // Add on the left path
        for &(next_kmer, dir) in path.iter() {
            let kmer = match dir {
                Dir::Left => next_kmer,
                Dir::Right => next_kmer.rc(),
            };

            edge_seq.push_front(kmer.get(0));

            // Reduce the data object
            let (_, kmer_data) = self.get_kmer_data(&next_kmer);
            node_data = self.spec.reduce(node_data, kmer_data)
        }

        let left_extend = match path.last() {
            None => l_ext,
            Some(&(_, Dir::Left)) => l_ext,
            Some(&(_, Dir::Right)) => l_ext.complement(),
        };


        // Unique path from seed kmer with Dir Right is built
        let r_ext = self.extend_kmer(seed, Dir::Right, path);

        // Add on the right path
        for &(next_kmer, dir) in path.iter() {
            let kmer = match dir {
                Dir::Left => next_kmer.rc(),
                Dir::Right => next_kmer,
            };

            edge_seq.push_back(kmer.get(K::k() - 1));

            let (_, kmer_data) = self.get_kmer_data(&next_kmer);
            node_data = self.spec.reduce(node_data, kmer_data)
        }

        let right_extend = match path.last() {
            None => r_ext,
            Some(&(_, Dir::Left)) => r_ext.complement(),
            Some(&(_, Dir::Right)) => r_ext,
        };
        
        (Exts::from_single_dirs(left_extend, right_extend), node_data)
    }


    /// Build the edge surrounding a kmer
    #[inline(never)]
    fn  build_node_start(
        &mut self,
        seed_id: usize,
        path: &mut Vec<(K, Dir)>,
        edge_seq: &mut VecDeque<u8>,
    ) -> (K, K) {
        let seed: K = *self.index.get_key(seed_id).expect("Index out of bound");
        edge_seq.clear();
        for i in 0..K::k() {
            edge_seq.push_back(seed.get(i));
        }

        let mut node_data = self.get_kmer_data(&seed).1.clone();

        // Unique path from seed kmer with Dir Left is built
        let _ = self.extend_kmer_par(seed, Dir::Left, path);

        // Add on the left path
        for &(next_kmer, dir) in path.iter() {
            let kmer = match dir {
                Dir::Left => next_kmer,
                Dir::Right => next_kmer.rc(),
            };

            edge_seq.push_front(kmer.get(0));

            // Reduce the data object
            let (_, kmer_data) = self.get_kmer_data(&next_kmer);
            node_data = self.spec.reduce(node_data, kmer_data)
        }

        // Unique path from seed kmer with Dir Right is built
        let _ = self.extend_kmer_par(seed, Dir::Right, path);

        // Add on the right path
        for &(next_kmer, dir) in path.iter() {
            let kmer = match dir {
                Dir::Left => next_kmer.rc(),
                Dir::Right => next_kmer,
            };

            edge_seq.push_back(kmer.get(K::k() - 1));

            let (_, kmer_data) = self.get_kmer_data(&next_kmer);
            node_data = self.spec.reduce(node_data, kmer_data)
        }

        let mut path_seq = PackedDnaStringSet::new();
        path_seq.add(edge_seq);
        let first = path_seq.sequence.get_kmer::<K>(0);
        let min_first = first.min_rc();
        let last = path_seq.sequence.get_kmer::<K>(path_seq.sequence.len()-K::k());
        let min_last = last.min_rc();

        if min_first < min_last {
            (min_first, min_last)
        } else {
            (min_last, min_first)
        }
    }

    #[inline(never)]
    fn  build_node_par(
        &mut self,
        seed_id: usize,
        path: &mut Vec<(K, Dir)>,
        edge_seq: &mut VecDeque<u8>,
    ) -> (Exts, D) {
        let seed: K = *self.index.get_key(seed_id).expect("Index out of bound");
        edge_seq.clear();
        for i in 0..K::k() {
            edge_seq.push_back(seed.get(i));
        }

        let mut node_data = self.get_kmer_data(&seed).1.clone();

        // Unique path from seed kmer with Dir Left is built
        let l_ext = self.extend_kmer_par(seed, Dir::Left, path);


        // Add on the left path
        for &(next_kmer, dir) in path.iter() {
            let kmer = match dir {
                Dir::Left => next_kmer,
                Dir::Right => next_kmer.rc(),
            };

            edge_seq.push_front(kmer.get(0));

            // Reduce the data object
            let (_, kmer_data) = self.get_kmer_data(&next_kmer);
            node_data = self.spec.reduce(node_data, kmer_data)
        }

        let left_extend = match path.last() {
            None => l_ext,
            Some(&(_, Dir::Left)) => l_ext,
            Some(&(_, Dir::Right)) => l_ext.complement(),
        };


        // Unique path from seed kmer with Dir Right is built
        let r_ext = self.extend_kmer_par(seed, Dir::Right, path);

        // Add on the right path
        for &(next_kmer, dir) in path.iter() {
            let kmer = match dir {
                Dir::Left => next_kmer.rc(),
                Dir::Right => next_kmer,
            };

            edge_seq.push_back(kmer.get(K::k() - 1));

            let (_, kmer_data) = self.get_kmer_data(&next_kmer);
            node_data = self.spec.reduce(node_data, kmer_data)
        }

        let right_extend = match path.last() {
            None => r_ext,
            Some(&(_, Dir::Left)) => r_ext.complement(),
            Some(&(_, Dir::Right)) => r_ext,
        };
        
        (Exts::from_single_dirs(left_extend, right_extend), node_data)
    }

    /// Compress a set of kmers and their extensions and metadata into a base DeBruijn graph.
    #[inline(never)]
    pub fn compress_kmers(
        stranded: bool,
        spec: &S,
        index: &BoomHashMap2<K, Exts, D>,
        progress: bool,
    ) -> BaseGraph<K, D> {
        
        let n_kmers = index.len();
        let mut available_kmers = BitSet::with_capacity(n_kmers);
        let progress = if n_kmers < 128 { false } else { progress };

        for i in 0..n_kmers {
            available_kmers.insert(i);
        }

        let mut comp = CompressFromHash {
            stranded,
            spec,
            k: PhantomData,
            d: PhantomData,
            available_kmers,
            index,
        };

        // Path-compressed De Bruijn graph will be created here
        let mut graph = BaseGraph::new(stranded);

        // Paths will be get assembled here
        let mut path_buf = Vec::new();

        // Node sequences will get assembled here
        let mut edge_seq_buf = VecDeque::new();

        debug!("n of kmers: {}", n_kmers);

        let steps = n_kmers as f32 / 128.;

        if progress {
            println!("Compressing kmers");
            for _i in 0..127 {
                print!("-");
            }
            println!("|");
        }

        let pb = ProgressBar::new(n_kmers as u64);
        pb.set_style(ProgressStyle::with_template("{msg} [{elapsed_precise}] {bar:60.cyan/blue} ({pos}/{len})").unwrap().progress_chars("#/-"));
        pb.set_message(format!("{:<32}", "compressing graph"));


        for kmer_counter in (0..n_kmers).progress_with(pb) {
            if progress && (kmer_counter as f32 % steps >= 0.) & (kmer_counter as f32 % steps < 1.) { print!("|")}

            if (kmer_counter as f32 % steps >= 0.) & (kmer_counter as f32 % steps < 1.) {
                debug!("another 1/128 done: {}, data graph size: {}", (kmer_counter as f32 / steps) as i32, mem::size_of_val(&*graph.data));
            }

            if comp.available_kmers.contains(kmer_counter) {
                let (node_exts, node_data) =
                    comp.build_node(kmer_counter, &mut path_buf, &mut edge_seq_buf);
                graph.add(&edge_seq_buf, node_exts, node_data);
            }
        }

        if progress { println!() };

        graph
    }

    /// Compress a set of kmers and their extensions and metadata into a base DeBruijn graph, utilizing multithreading
    #[inline(never)]
    pub fn compress_kmers_parallel(
        stranded: bool,
        spec: &S,
        index: &BoomHashMap2<K, Exts, D>,
        progress: bool,
    ) -> BaseGraph<K, D> {
        
        let n_kmers = index.len();
        let progress = if n_kmers < 128 { false } else { progress };

        // calculate the slices of the workload according to the available threads
        debug!("current num threads: {}", current_num_threads());
        debug!("n of kmers: {}", n_kmers);

        let slices = current_num_threads();
        let sz = n_kmers / slices + 1;

        debug!("n_kmers: {}", n_kmers);
        debug!("sz: {}", sz);

        let mut parallel_ranges = Vec::with_capacity(slices);
        let mut start = 0;
        while start < n_kmers {
            parallel_ranges.push(start..start + sz);
            start += sz;
        }

        let last_start = parallel_ranges.pop().expect("no kmers in parallel ranges").start;
        parallel_ranges.push(last_start..n_kmers);
        debug!("parallel ranges: {:?}", parallel_ranges);

        let all_start_end_kmers: Arc<Mutex<HashMap<K, K>>> = Arc::new(Mutex::new(HashMap::new()));

        let steps = n_kmers as f32 / 128.;
        if progress {
            println!("Compressing kmers\nStep 1 of 2");
            for _i in 0..127 {
                print!("-");
            }
            println!("|");
        }
        
        // go through all kmers and find the start and end kmer of each compressable sequence node
        parallel_ranges.into_par_iter().for_each(|range| {

            let mut path_buf = Vec::new();
            let mut edge_seq_buf = VecDeque::new();
            //let mut start_end_kmers: Vec<(K, K)> = Vec::new();

            let range2 = range.clone();

            for kmer_counter in range {

                if progress && (kmer_counter as f32 % steps >= 0.) & (kmer_counter as f32 % steps < 1.) { print!("|")}

                let mut comp = CompressFromHash {
                    stranded,
                    spec,
                    k: PhantomData,
                    d: PhantomData,
                    available_kmers: BitSet::new(),
                    index,
                };               

                let kmers = comp.build_node_start(kmer_counter, &mut path_buf, &mut edge_seq_buf);

                let all_clone = Arc::clone(&all_start_end_kmers);
                let mut all_lock = all_clone.lock().expect("lock all_start_end_kmers");
                all_lock.entry(kmers.0).or_insert(kmers.1);
            }

            debug!("finished range: {:?}", range2);
        });

        if progress { println!() }

        // all the kmers which occurr at the beginning or end of a node are soerted and deduped
        // resulting in one (K, K) per node
        let all_start_end_kmers: Vec<(K, K)> = all_start_end_kmers.lock().expect("final lock all_start_end_kmers").clone().into_iter().collect();

        //all_start_end_kmers.sort();
        //all_start_end_kmers.dedup();

        //debug!("all start and end kmers: {:?}", all_start_end_kmers);

        // all graphs constructed by the respective threads are pushed into this and later combined
        let graphs: Arc<Mutex<Vec<BaseGraph<K, D>>>> = Arc::new(Mutex::new(Vec::with_capacity(current_num_threads())));
        
        // calculate the workload has to be sliced depending on the threads
        let n_starts = all_start_end_kmers.len();
        let slices = current_num_threads();
        let sz = n_starts / slices + 1;

        debug!("n start kmers: {}", all_start_end_kmers.len());
        debug!("sz: {}", sz);

        let mut parallel_ranges = Vec::with_capacity(slices);
        let mut start = 0;
        while start < n_starts {
            parallel_ranges.push(start..start + sz);
            start += sz;
        }

        let last_start = parallel_ranges.pop().expect("no kmers in parallel ranges").start;
        parallel_ranges.push(last_start..n_starts);

        debug!("parallel ranges 2: {:?}", parallel_ranges);

        let steps = n_starts as f32 / 128.;
        if progress {
            println!("Step 2 of 2");
            for _i in 0..127 {
                print!("-");
            }
            println!("|");
        }
        let short_progress = n_starts < 100;

        // go trough all start kmers and find the corresponding sequence, add them to the partial graph
        parallel_ranges.into_par_iter().for_each(|range| {

            let mut graph: BaseGraph<K, D> = BaseGraph::new(stranded);

            let range2 = range.clone();

            for (i, (start, _)) in all_start_end_kmers[range].iter().enumerate() {   

                if progress & !short_progress && (i as f32 % steps >= 0.) & (i as f32 % steps < 1.) { print!("|")}

                let mut path_buf = Vec::new();
                let mut edge_seq_buf = VecDeque::new();
                let mut comp = CompressFromHash {
                    stranded,
                    spec,
                    k: PhantomData,
                    d: PhantomData,
                    available_kmers: BitSet::new(),
                    index,
                };      
                let (node_exts, node_data) = comp.build_node_par(comp.index.get_key_id(start).expect("get kmer id from index, should exist"), &mut path_buf, &mut edge_seq_buf);
                graph.add(&edge_seq_buf, node_exts, node_data);
            }

            debug!("finished range: {:?}, subgraph data size: {}", range2, mem::size_of_val(&*graph.data));

            let graphs_clone = Arc::clone(&graphs);
            let mut graphs_lock = graphs_clone.lock().expect("lock graphs to push new graph");
            graphs_lock.push(graph);

        });

        if progress & short_progress {
            for _i in 0..128 {
                print!("|");
            }
        }
        
        if progress { println!() }

        // combine the graphs and return the resulting graph
        let graph = BaseGraph::combine(graphs.lock().expect("final graph lock").clone().into_iter());
        graph
    }
}


/// Take a BoomHash Object and build a compressed DeBruijn graph.
#[inline(never)]
pub fn compress_kmers_with_hash<K: Kmer + Send + Sync, D: Clone + Debug + Send + Sync, S: CompressionSpec<D> + Send + Sync>(
    stranded: bool,
    spec: &S,
    index: &BoomHashMap2<K, Exts, D>,
    time: bool,
    parallel: bool,
    progress: bool,
) -> BaseGraph<K, D> {
    let before_compression = Instant::now();
    let graph = if !parallel {CompressFromHash::<K, D, S>::compress_kmers(stranded, spec, index, progress) } else {
        CompressFromHash::<K, D, S>::compress_kmers_parallel(stranded, spec, index, progress)};
    if time { println!("time compression (s): {}", before_compression.elapsed().as_secs_f32()) }
    graph
}

/// Take (make) a BoomHash Object and build a compressed DeBruijn graph.
#[inline(never)]
pub fn compress_kmers<K: Kmer + Send + Sync, D: Clone + Debug  + Send + Sync, S: CompressionSpec<D> + Send + Sync>(
    stranded: bool,
    spec: &S,
    kmer_exts: &[(K, (Exts, D))],
) -> BaseGraph<K, D> {
    let mut keys = Vec::with_capacity(kmer_exts.len());
    let mut exts = Vec::with_capacity(kmer_exts.len());
    let mut data = Vec::with_capacity(kmer_exts.len());

    for (k, (e, d)) in kmer_exts {
        keys.push(*k);
        data.push(d.clone());
        exts.push(*e);
    }

    let index = BoomHashMap2::new(keys, exts, data);
    CompressFromHash::<K, D, S>::compress_kmers(stranded, spec, &index, false)
}

/// Build graph from a set of kmers with unknown extensions by finding the extensions on the fly.
#[inline(never)]
pub fn compress_kmers_no_exts<K: Kmer + Send + Sync, D: Clone + Debug + Send + Sync, S: CompressionSpec<D> + Send + Sync>(
    stranded: bool,
    spec: &S,
    kmer_exts: &[(K, D)],
) -> BaseGraph<K, D> {
    let kmer_set: std::collections::HashSet<_> = kmer_exts.iter().map(|(k, _)| k).collect();

    let can = |k: K| k.min_rc();

    let mut keys = Vec::with_capacity(kmer_exts.len());
    let mut exts = Vec::with_capacity(kmer_exts.len());
    let mut data = Vec::with_capacity(kmer_exts.len());
    for (k, d) in kmer_exts {
        let mut e = Exts::empty();

        for l in 0..4 {
            let new = can(k.extend_left(l));

            if kmer_set.contains(&new) {
                e = e.set(Dir::Left, l);
            }
        }

        for r in 0..4 {
            let new = can(k.extend_right(r));

            if kmer_set.contains(&new) {
                e = e.set(Dir::Right, r);
            }
        }

        keys.push(*k);
        data.push(d.clone());
        exts.push(e);
    }

    assert_eq!(kmer_set.len(), keys.len());

    let index = BoomHashMap2::new(keys, exts, data);
    CompressFromHash::<K, D, S>::compress_kmers(stranded, spec, &index,false)
}

/// assumes stranded = false
pub fn uncompressed_graph<K: Kmer, D: Clone + Debug>(
    index: &BoomHashMap2<K, Exts, D>,
) -> BaseGraph<K, D> {

    let mut graph: BaseGraph<K, D> = BaseGraph::new(false);
    let mut kmer_seq: VecDeque<u8> = VecDeque::with_capacity(K::k());

    //let n_kmers = index.len();

    for (kmer, exts, data) in index.into_iter() {
        kmer_seq.clear();
        for i in 0..K::k() {
            kmer_seq.push_back(kmer.get(i));
        }
        let data2 = data.clone();
        graph.add(&kmer_seq, *exts, data2);


    }
    graph
}
