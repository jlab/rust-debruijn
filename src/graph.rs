// Copyright 2017 10x Genomics

//! Containers for path-compressed De Bruijn graphs

use bit_set::BitSet;
use indicatif::ProgressBar;
use indicatif::ProgressIterator;
use indicatif::ProgressStyle;
use itertools::enumerate;
use log::{debug, trace};
use rayon::prelude::*;
use rayon::current_num_threads;
use serde_derive::{Deserialize, Serialize};
use smallvec::SmallVec;
use std::borrow::Borrow;

use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;
use std::f32;
use std::fmt::{self, Debug, Display};
use std::fs::{remove_file, File};
use std::hash::Hash;
use std::io::BufWriter;
use std::io::{BufReader, Error, Read};
use std::io::Write;
use std::iter::FromIterator;
use std::marker::PhantomData;
use std::path::Path;

use boomphf::hashmap::BoomHashMap;

use serde_json;
use serde_json::Value;

type SmallVec4<T> = SmallVec<[T; 4]>;
type SmallVec8<T> = SmallVec<[T; 8]>;

use crate::bits_to_base;
use crate::colors::ColorMode;
use crate::colors::Colors;
use crate::compression::CompressionSpec;
use crate::dna_string::{DnaString, DnaStringSlice, PackedDnaStringSet};
use crate::summarizer::SummaryConfig;
use crate::summarizer::SummaryData;
use crate::summarizer::ID;
use crate::BUF;
use crate::PROGRESS_STYLE;
use crate::{Dir, Exts, Kmer, Mer, Vmer};

/// A compressed DeBruijn graph carrying auxiliary data on each node of type `D`.
/// This type does not carry the sorted index arrays the allow the graph
/// to be walked efficiently. The `DeBruijnGraph` type wraps this type and add those
/// vectors.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct BaseGraph<K, D> {
    pub sequences: PackedDnaStringSet,
    pub exts: Vec<Exts>,
    pub data: Vec<D>,
    pub stranded: bool,
    phantom: PhantomData<K>,
}

impl<K, D> BaseGraph<K, D> {
    pub fn new(stranded: bool) -> Self {
        BaseGraph {
            sequences: PackedDnaStringSet::new(),
            exts: Vec::new(),
            data: Vec::new(),
            phantom: PhantomData,
            stranded,
        }
    }

    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }

    pub fn combine<I: Iterator<Item = BaseGraph<K, D>>>(graphs: I) -> Self {
        let mut sequences = PackedDnaStringSet::new();
        let mut exts = Vec::new();
        let mut data = Vec::new();
        let mut stranded = Vec::new();

        for g in graphs {
            for s in 0..g.sequences.len() {
                sequences.add(&g.sequences.get(s));
            }

            exts.extend(g.exts);
            data.extend(g.data);
            stranded.push(g.stranded);
        }

        let out_stranded = stranded.iter().all(|x| *x);

        if !out_stranded && !stranded.iter().all(|x| !*x) {
            panic!("attempted to combine stranded and unstranded graphs");
        }

        BaseGraph {
            sequences,
            stranded: out_stranded,
            exts,
            data,
            phantom: PhantomData,
        }
    }
}

impl<K: Kmer, D> BaseGraph<K, D> {
    pub fn add<R: Borrow<u8>, S: IntoIterator<Item = R>>(
        &mut self,
        sequence: S,
        exts: Exts,
        data: D,
    ) {
        self.sequences.add(sequence);
        self.exts.push(exts);
        self.data.push(data);
    }
}

impl<K: Kmer + Send + Sync, D> BaseGraph<K, D> {
    pub fn finish(self) -> DebruijnGraph<K, D> {
        let indices: Vec<u32> = (0..self.len() as u32).collect();

        let left_order = {
            let mut kmers: Vec<K> = Vec::with_capacity(self.len());
            for idx in &indices {
                kmers.push(self.sequences.get(*idx as usize).first_kmer());
                
            }
            
            BoomHashMap::new_parallel(kmers, indices.clone())
        };

        let right_order = {
            let mut kmers: Vec<K> = Vec::with_capacity(self.len());
            for idx in &indices {
                kmers.push(self.sequences.get(*idx as usize).last_kmer());
            }
            
            BoomHashMap::new_parallel(kmers, indices)
        };
        debug!("finish graph loops: 2x {}", self.len());

        DebruijnGraph {
            base: self,
            left_order,
            right_order,
        }
    }
}

impl<K: Kmer, D> BaseGraph<K, D> {
    pub fn finish_serial(self) -> DebruijnGraph<K, D> {
        let indices: Vec<u32> = (0..self.len() as u32).collect();

        let left_order = {
            let mut kmers: Vec<K> = Vec::with_capacity(self.len());
            let mut sequences: Vec<String> = Vec::new();
            for idx in &indices {
                kmers.push(self.sequences.get(*idx as usize).first_kmer());
                sequences.push(self.sequences.get(*idx as usize).to_dna_string());
            }
            println!("left kmers: {:?}", kmers);
            println!("left seqs: {:?}", sequences);
            BoomHashMap::new(kmers, indices.clone())
        };

        let right_order = {
            let mut kmers: Vec<K> = Vec::with_capacity(self.len());
            let mut sequences: Vec<String> = Vec::new();
            for idx in &indices {
                kmers.push(self.sequences.get(*idx as usize).last_kmer());
                sequences.push(self.sequences.get(*idx as usize).to_dna_string());
            }
            println!("right kmers: {:?}", kmers);
            println!("right seqs: {:?}", sequences);
            BoomHashMap::new(kmers, indices)
        };

        DebruijnGraph {
            base: self,
            left_order,
            right_order,
        }
    }
}

/// A compressed DeBruijn graph carrying auxiliary data on each node of type `D`.
/// The struct carries sorted index arrays the allow the graph
/// to be walked efficiently.
#[derive(Serialize, Deserialize, Debug)]
pub struct DebruijnGraph<K: Hash, D> {
    pub base: BaseGraph<K, D>,
    left_order: BoomHashMap<K, u32>,
    right_order: BoomHashMap<K, u32>,
}

impl<K: Kmer, D: Debug> DebruijnGraph<K, D> {
    /// Total number of nodes in the DeBruijn graph
    pub fn len(&self) -> usize {
        self.base.len()
    }

    pub fn is_empty(&self) -> bool {
        self.base.is_empty()
    }

    /// Get a node given it's `node_id`
    pub fn get_node(&self, node_id: usize) -> Node<K, D> {
        Node {
            node_id,
            graph: self,
        }
    }

    /// Get a node given it's `node_id`
    pub fn get_node_kmer(&self, node_id: usize) -> NodeKmer<K, D> {
        let node = self.get_node(node_id);
        let node_seq = node.sequence();

        NodeKmer {
            node_id,
            node_seq_slice: node_seq,
            phantom_d: PhantomData,
            phantom_k: PhantomData,
        }
    }

    /// Return an iterator over all nodes in the graph
    pub fn iter_nodes(&self) -> NodeIter<K, D> {
        NodeIter {
            graph: self,
            node_id: 0,
        }
    }

    /// Find the edges leaving node `node_id` in direction `Dir`. Should generally be
    /// accessed via a Node wrapper object
    fn find_edges(&self, node_id: usize, dir: Dir) -> SmallVec4<(u8, usize, Dir, bool)> {
        let exts = self.base.exts[node_id];
        let sequence = self.base.sequences.get(node_id);
        let kmer: K = sequence.term_kmer(dir);
        let mut edges = SmallVec4::new();

        for i in 0..4 {
            if exts.has_ext(dir, i) {
                let link = self.find_link(kmer.extend(i, dir), dir).expect("missing link");
                edges.push((i, link.0, link.1, link.2));
            }
        }

        edges
    }

    /// Find the edges leaving node `node_id` in direction `Dir`. Should generally be
    /// accessed via a Node wrapper object
    /// 
    /// allows missing links
    fn _find_edges_sharded(&self, node_id: usize, dir: Dir) -> SmallVec4<(u8, usize, Dir, bool)> {
        let exts = self.base.exts[node_id];
        let sequence = self.base.sequences.get(node_id);
        let kmer: K = sequence.term_kmer(dir);
        let mut edges = SmallVec4::new();

        for i in 0..4 {
            if exts.has_ext(dir, i) {
                let link = self.find_link(kmer.extend(i, dir), dir); //.expect("missing link");
                if let Some(l) = link {
                    edges.push((i, l.0, l.1, l.2));
                }
                // Otherwise, this edge doesn't exist within this shard, so ignore it.
                // NOTE: this should be allowed in a 'complete' DBG
            }
        }

        edges
    }

    /// Seach for the kmer `kmer`, appearing at the given `side` of a node sequence.
    fn search_kmer(&self, kmer: K, side: Dir) -> Option<usize> {
        match side {
            Dir::Left => self.left_order.get(&kmer).map(|pos| *pos as usize),
            Dir::Right => self.right_order.get(&kmer).map(|pos| *pos as usize),
        }
    }

    /// Find a link in the graph, possibly handling a RC switch.
    pub fn find_link(&self, kmer: K, dir: Dir) -> Option<(usize, Dir, bool)> {
        // Only test self-consistent paths through
        // the edges
        // Avoids problems due to single kmer edges
        // (dir, next_side_incoming, rc)
        // (Dir::Left, Dir::Right, false) => true,
        // (Dir::Left, Dir::Left,  true) => true,
        // (Dir::Right, Dir::Left, false) => true,
        // (Dir::Right, Dir::Right, true) => true,

        let rc = kmer.rc();

        match dir {
            Dir::Left => {
                if let Some(idx) = self.search_kmer(kmer, Dir::Right) {
                    return Some((idx, Dir::Right, false));
                }

                if !self.base.stranded {
                    if let Some(idx) = self.search_kmer(rc, Dir::Left) {
                        return Some((idx, Dir::Left, true));
                    }
                }
            }

            Dir::Right => {
                if let Some(idx) = self.search_kmer(kmer, Dir::Left) {
                    return Some((idx, Dir::Left, false));
                }

                if !self.base.stranded {
                    if let Some(idx) = self.search_kmer(rc, Dir::Right) {
                        return Some((idx, Dir::Right, true));
                    }
                }
            }
        }

        None
    }

    /// Check whether the graph is fully compressed. Return `None` if it's compressed,
    /// otherwise return `Some(node1, node2)` representing a pair of node that could
    /// be collapsed. Probably only useful for testing.
    pub fn is_compressed<S: CompressionSpec<D>>(&self, spec: &S) -> Option<(usize, usize)> {
        for i in 0..self.len() {
            let n = self.get_node(i);

            for dir in &[Dir::Left, Dir::Right] {
                let dir_edges = n.edges(*dir);
                if dir_edges.len() == 1 {
                    let (_, next_id, return_dir, _) = dir_edges[0];
                    let next = self.get_node(next_id);

                    let ret_edges = next.edges(return_dir);
                    if ret_edges.len() == 1 {
                        // Test for us being a palindrome: this makes it OK
                        if n.len() == K::k() && n.sequence().first_kmer::<K>().is_palindrome() {
                            continue;
                        }

                        // Test for a neighbor being a palindrome: this makes it OK
                        if next.len() == K::k() && next.sequence().first_kmer::<K>().is_palindrome()
                        {
                            continue;
                        }

                        // Test for this edge representing a smooth circle (biting it's own tail):
                        if n.node_id == next_id {
                            continue;
                        }

                        if spec.join_test(n.data(), next.data()) {
                            // Found a unbranched edge that should have been eliminated
                            return Some((i, next_id));
                        }
                    }
                }
            }
        }

        None
    }

    /// Remove non-existent extensions that may be created due to filtered kmers
    /// 
    /// if `valid_nodes` if `None`, all nodes are valid
    pub fn fix_exts(&mut self, valid_nodes: Option<&BitSet>) {
        for i in 0..self.len() {
            let valid_exts = self.get_valid_exts(i, valid_nodes);
            self.base.exts[i] = valid_exts;
        }
    }

    pub fn get_valid_exts(&self, node_id: usize, valid_nodes: Option<&BitSet>) -> Exts {
        let mut new_exts = Exts::empty();
        let node = self.get_node(node_id);
        let exts = node.exts();
        let l_kmer: K = node.sequence().first_kmer();
        let r_kmer: K = node.sequence().last_kmer();

        let check_node = |id| match valid_nodes {
            Some(bs) => bs.contains(id),
            None => true,
        };

        for i in 0..4 {
            if exts.has_ext(Dir::Left, i) {
                match self.find_link(l_kmer.extend_left(i), Dir::Left) {
                    Some((target, _, _)) if check_node(target) => {
                        new_exts = new_exts.set(Dir::Left, i)
                    }
                    _ => (),
                }
            }

            if exts.has_ext(Dir::Right, i) {
                match self.find_link(r_kmer.extend_right(i), Dir::Right) {
                    Some((target, _, _)) if check_node(target) => {
                        new_exts = new_exts.set(Dir::Right, i)
                    }
                    _ => (),
                }
            }
        }

        new_exts
    }

    /// Find the highest-scoring, unambiguous path in the graph. Each node get a score
    /// given by `score`. Any node where `solid_path(node) == True` are valid paths -
    /// paths will be terminated if there are multiple valid paths emanating from a node.
    pub fn max_path<F, F2>(&self, score: F, solid_path: F2) -> Vec<(usize, Dir)>
    where
        F: Fn(&D) -> f32,
        F2: Fn(&D) -> bool,
    {
        if self.is_empty() {
            return Vec::default();
        }

        let mut best_node = 0;
        let mut best_score = f32::MIN;
        for i in 0..self.len() {
            let node = self.get_node(i);
            let node_score = score(node.data());

            if node_score > best_score {
                best_node = i;
                best_score = node_score;
            }
        }

        let oscore = |state| match state {
            None => 0.0,
            Some((id, _)) => score(self.get_node(id).data()),
        };

        let osolid_path = |state| match state {
            None => false,
            Some((id, _)) => solid_path(self.get_node(id).data()),
        };

        // Now expand in each direction, greedily taking the best path. Stop if we hit a node we've
        // already put into the path
        let mut used_nodes = HashSet::new();
        let mut path = VecDeque::new();

        // Start w/ initial state
        used_nodes.insert(best_node);
        path.push_front((best_node, Dir::Left));

        for init in [(best_node, Dir::Left, false), (best_node, Dir::Right, true)].iter() {
            let &(start_node, dir, do_flip) = init;
            let mut current = (start_node, dir);

            loop {
                let mut next = None;
                let (cur_id, incoming_dir) = current;
                let node = self.get_node(cur_id);
                let edges = node.edges(incoming_dir.flip());

                let mut solid_paths = 0;
                for (_, id, dir, _) in edges {
                    let cand = Some((id, dir));
                    if osolid_path(cand) {
                        solid_paths += 1;
                    }

                    if oscore(cand) > oscore(next) {
                        next = cand;
                    }
                }

                if solid_paths > 1 {
                    break;
                }

                match next {
                    Some((next_id, next_incoming)) if !used_nodes.contains(&next_id) => {
                        if do_flip {
                            path.push_front((next_id, next_incoming.flip()));
                        } else {
                            path.push_back((next_id, next_incoming));
                        }

                        used_nodes.insert(next_id);
                        current = (next_id, next_incoming);
                    }
                    _ => break,
                }
            }
        }

        Vec::from_iter(path)
    }


    /// Find the highest-scoring, unambiguous path in the graph. Each node get a score
    /// given by `score`. Any node where `solid_path(node) == True` are valid paths -
    /// paths will be terminated if there are multiple valid paths emanating from a node.
    /// Returns vec with path for each component
    pub fn max_path_comp<F, F2>(&self, score: F, solid_path: F2) -> Vec<VecDeque<(usize, Dir)>>
    where
        F: Fn(&D) -> f32,
        F2: Fn(&D) -> bool,
    {
        if self.is_empty() {
            let vec: Vec<VecDeque<(usize, Dir)>> = Vec::new();
            return vec;
        }

        let components = self.iter_components();
        let mut paths: Vec<VecDeque<(usize, Dir)>> = Vec::new();

        for component in components {

            let current_comp = &component;
            

            let mut best_node = current_comp[0];
            let mut best_score = f32::MIN;
            for c in current_comp.iter() {
                let node = self.get_node(*c);
                let node_score = score(node.data());

                if node_score > best_score {
                    best_node = *c;
                    best_score = node_score;
                }
            }

            let oscore = |state| match state {
                None => 0.0,
                Some((id, _)) => score(self.get_node(id).data()),
            };

            let osolid_path = |state| match state {
                None => false,
                Some((id, _)) => solid_path(self.get_node(id).data()),
            };

            // Now expand in each direction, greedily taking the best path. Stop if we hit a node we've
            // already put into the path
            let mut used_nodes = HashSet::new();
            let mut path = VecDeque::new();

            // Start w/ initial state
            used_nodes.insert(best_node);
            path.push_front((best_node, Dir::Left));

            for init in [(best_node, Dir::Left, false), (best_node, Dir::Right, true)].iter() {
                let &(start_node, dir, do_flip) = init;
                let mut current = (start_node, dir);

                loop {
                    let mut next = None;
                    let (cur_id, incoming_dir) = current;
                    let node = self.get_node(cur_id);
                    let edges = node.edges(incoming_dir.flip());

                    let mut solid_paths = 0;
                    for (_, id, dir, _) in edges {
                        let cand = Some((id, dir));
                        if osolid_path(cand) {
                            solid_paths += 1;
                        }

                        if oscore(cand) > oscore(next) {
                            next = cand;
                        }
                    }

                    if solid_paths > 1 {
                        break;
                    }

                    match next {
                        Some((next_id, next_incoming)) if !used_nodes.contains(&next_id) => {
                            if do_flip {
                                path.push_front((next_id, next_incoming.flip()));
                            } else {
                                path.push_back((next_id, next_incoming));
                            }

                            used_nodes.insert(next_id);
                            current = (next_id, next_incoming);
                        }
                        _ => break,
                    }
                }
            }
            
            paths.push(path);
            //paths.push(Vec::from_iter(path));
        }

        paths
    
    }

    pub fn iter_max_path_comp<F, F2>(&self, score: F, solid_path: F2) -> PathCompIter<K, D, F, F2> 
    where 
    F: Fn(&D) -> f32,
    F2: Fn(&D) -> bool
    {
        let component_iterator = self.iter_components();
        PathCompIter { graph: self, component_iterator, graph_pos: 0, score, solid_path }
    }

    /// write the paths from `iter_max_path_comp` to a fasta file
    pub fn path_to_fasta<F, F2>(&self, f: &mut dyn std::io::Write, path_iter: PathCompIter<K, D, F, F2>, return_lens: bool) -> (Vec<usize>, Vec<usize>)
    where 
    F: Fn(&D) -> f32,
    F2: Fn(&D) -> bool
    {
        // width of fasta lines
        let columns = 80;

        // sizes of components and of paths
        let mut comp_sizes = Vec::new();
        let mut path_lens = Vec::new();

        for (seq_counter, (component, path)) in path_iter.enumerate() {
            // get dna sequence from path
            let seq = self.sequence_of_path(path.iter());

            //write header with length & start node
            writeln!(f, ">path{} len={} start_node={}", seq_counter, seq.len(), path[0].0).unwrap();

            // calculate how sequence has to be split up
            let slices = (seq.len() / columns) + 1;
            let mut ranges = Vec::with_capacity(slices);

            let mut start = 0;
            while start < seq.len() {
                ranges.push(start..start + columns);
                start += columns;
            }

            let last_start = ranges.pop().expect("no kmers in parallel ranges").start;
            ranges.push(last_start..seq.len());

            // split up sequence and write to file accordingly
            for range in ranges {
                writeln!(f, "{:?}", seq.slice(range.start, range.end)).unwrap();
            }

            if return_lens {
                comp_sizes.push(component.len());
                path_lens.push(path.len());
            }
        }    

        (comp_sizes, path_lens)
        
    }


    /// Get the sequence of a path through the graph. The path is given as a sequence of node_id integers
    pub fn sequence_of_path<'a, I: 'a + Iterator<Item = &'a (usize, Dir)>>(
        &self,
        path: I,
    ) -> DnaString {
        let mut seq = DnaString::new();

        for (idx, &(node_id, dir)) in path.enumerate() {
            let start = if idx == 0 { 0 } else { K::k() - 1 };

            let node_seq = match dir {
                Dir::Left => self.get_node(node_id).sequence(),
                Dir::Right => self.get_node(node_id).sequence().rc(),
            };

            for p in start..node_seq.len() {
                seq.push(node_seq.get(p))
            }
        }

        seq
    }

    /// write a node to a dot file
    /// 
    /// ### Arguments: 
    /// 
    /// * `node`: [`Node<K, D>`] which will be written to a dot file
    /// * `node_label`: closure taking [`Node<K, D>`] and returning a string containing commands for dot nodes 
    /// * `edge_label`: closure taking [`Node<K, D>`], the base as a [`u8`], the incoming [`Dir`] of the edge 
    ///    and if the neighbor is flipped - returns a string containing commands for dot edges, 
    /// * `f`: writer
    fn node_to_dot<FN: Fn(&Node<K, D>) -> String, FE: Fn(&Node<K, D>, u8, Dir, bool) -> String>(
        &self,
        node: &Node<'_, K, D>,
        node_label: &FN,
        edge_label: &FE,
        f: &mut dyn Write,
    ) {
        writeln!(f, "n{} {}", node.node_id, node_label(node)).unwrap();
        assert_eq!(node.exts().val.count_ones() as usize, node.l_edges().len() + node.r_edges().len());

        for (base, id, incoming_dir, flipped) in node.l_edges() {
            writeln!(f, "n{} -> n{} {}", id, node.node_id, edge_label(node, base, incoming_dir, flipped)).unwrap();
        }

        for (base, id, incoming_dir, flipped) in node.r_edges() {
            writeln!(f, "n{} -> n{} {}", node.node_id, id, edge_label(node, base, incoming_dir, flipped)).unwrap();
        }
    }

    /// Write the graph to a dot file
    /// 
    /// ### Arguments: 
    /// 
    /// * `path`: path to the output file
    /// * `node_label`: closure taking [`Node<K, D>`] and returning a string containing commands for dot nodes 
    /// * `edge_label`: closure taking [`Node<K, D>`], the base as a [`u8`], the incoming [`Dir`] of the edge 
    ///    and if the neighbor is flipped - returns a string containing commands for dot edges, 
    pub fn to_dot<P, FN, FE>(&self, path: P, node_label: &FN, edge_label: &FE) 
    where 
    P: AsRef<Path>,
    FN: Fn(&Node<K, D>) -> String,
    FE: Fn(&Node<K, D>, u8, Dir, bool) -> String,
    {
        let mut f = BufWriter::with_capacity(BUF, File::create(path).expect("error creating dot file"));

        let pb = ProgressBar::new(self.len() as u64);
        pb.set_style(ProgressStyle::with_template(PROGRESS_STYLE).unwrap().progress_chars("#/-"));
        pb.set_message(format!("{:<32}", "writing graph to DOT file"));

        writeln!(&mut f, "digraph {{\nrankdir=\"LR\"\nmodel=subset\noverlap=scalexy").unwrap();
        for i in (0..self.len()).progress_with(pb) {
            self.node_to_dot(&self.get_node(i), node_label, edge_label, &mut f);
        }
        writeln!(&mut f, "}}").unwrap();
        
        f.flush().unwrap();
        debug!("large to dot loop: {}", self.len());
    }

    /// Write the graph to a dot file in parallel
    /// Will write in to n_threads files simultaniously,
    /// then go though the files and add the contents to a larger file, 
    /// and delete the small files.
    /// 
    /// ### Arguments: 
    /// 
    /// * `path`: path to the output file
    /// * `node_label`: closure taking [`Node<K, D>`] and returning a string containing commands for dot nodes 
    /// * `edge_label`: closure taking [`Node<K, D>`], the base as a [`u8`], the incoming [`Dir`] of the edge 
    ///    and if the neighbor is flipped - returns a string containing commands for dot edges, 
    pub fn to_dot_parallel<P, FN, FE>(&self, path: P, node_label: &FN, edge_label: &FE) 
    where 
        D: Sync,
        K: Sync,
        P: AsRef<Path> + Display + Sync,
        FN: Fn(&Node<K, D>) -> String + Sync,
        FE: Fn(&Node<K, D>, u8, Dir, bool) -> String + Sync,
    {        
        let slices = current_num_threads();
        let n_nodes = self.len();
        let sz = n_nodes / slices + 1;

        debug!("n_nodes: {}", n_nodes);
        debug!("sz: {}", sz);

        let mut parallel_ranges = Vec::with_capacity(slices);
        let mut start = 0;
        while start < n_nodes {
            parallel_ranges.push(start..start + sz);
            start += sz;
        }

        let last_start = parallel_ranges.pop().expect("no kmers in parallel ranges").start;
        parallel_ranges.push(last_start..n_nodes);
        debug!("parallel ranges: {:?}", parallel_ranges);

        let mut files = Vec::with_capacity(current_num_threads());

        for i in 0..parallel_ranges.len() {
            files.push(format!("{}-{}.dot", path, i));
        } 

        let pb = ProgressBar::new(self.len() as u64);
        pb.set_style(ProgressStyle::with_template(PROGRESS_STYLE).unwrap().progress_chars("#/-"));
        pb.set_message(format!("{:<32}", "writing graph to DOT files"));
    
        parallel_ranges.into_par_iter().enumerate().for_each(|(i, range)| {
            let mut f = BufWriter::with_capacity(BUF, File::create(&files[i]).expect("error creating parallel dot file"));

            for i in range {
                self.node_to_dot(&self.get_node(i), node_label, edge_label, &mut f);
                pb.inc(1);
            }

            f.flush().unwrap();
        });
        pb.finish_and_clear();

        let mut out_file = BufWriter::with_capacity(BUF, File::create(path).expect("error creating combined dot file"));

        writeln!(&mut out_file, "digraph {{\nrankdir=\"LR\"\nmodel=subset\noverlap=scalexy").unwrap();

        let pb = ProgressBar::new(files.len() as u64);
        pb.set_style(ProgressStyle::with_template(PROGRESS_STYLE).unwrap().progress_chars("#/-"));
        pb.set_message(format!("{:<32}", "combining files"));

        for file in files.iter().progress_with(pb) {
            let open_file = File::open(file).expect("error opening parallel dot file");
            let mut reader = BufReader::new(open_file);
            let mut buffer = [0; BUF];

            loop {
                let linecount = reader.read(&mut buffer).unwrap();
                if linecount == 0 { break }
                out_file.write_all(&buffer[..linecount]).unwrap();
            }

            remove_file(file).unwrap();
        }

        writeln!(&mut out_file, "}}").unwrap();

        out_file.flush().unwrap();


    }

    /// Write part of the graph to a dot file
    /// 
    /// ### Arguments: 
    /// 
    /// * `path`: path to the output file
    /// * `node_label`: closure taking [`Node<K, D>`] and returning a string containing commands for dot nodes 
    /// * `edge_label`: closure taking [`Node<K, D>`], the base as a [`u8`], the incoming [`Dir`] of the edge 
    ///    and if the neighbor is flipped - returns a string containing commands for dot edges, 
    /// * `nodes`: [`Vec<usize>`] listing all IDs of nodes which should be included
    pub fn to_dot_partial<P, FN, FE>(&self, path: P, node_label: &FN, edge_label: &FE, nodes: Vec<usize>) 
    where 
        P: AsRef<Path>,
        FN: Fn(&Node<K, D>) -> String,
        FE: Fn(&Node<K, D>, u8, Dir, bool) -> String,
    {
        let mut f = BufWriter::with_capacity(BUF, File::create(path).expect("error creating dot file"));

        let pb = ProgressBar::new(nodes.len() as u64);
        pb.set_style(ProgressStyle::with_template(PROGRESS_STYLE).unwrap().progress_chars("#/-"));
        pb.set_message(format!("{:<32}", "writing graph to DOT file"));

        writeln!(&mut f, "# {:?}", nodes).unwrap();
        writeln!(&mut f, "digraph {{\nrankdir=\"LR\"\nmodel=subset\noverlap=scalexy").unwrap();
        for i in nodes.into_iter().progress_with(pb) {
            self.node_to_dot(&self.get_node(i), node_label, edge_label, &mut f);
        }
        writeln!(&mut f, "}}").unwrap();

        f.flush().unwrap();

        debug!("large to dot loop: {}", self.len());
    }

    fn node_to_gfa<F: Fn(&Node<'_, K, D>) -> String>(
        &self,
        node: &Node<'_, K, D>,
        w: &mut dyn Write,
        tag_func: Option<&F>,
    ) -> Result<(), Error> {
        match tag_func {
            Some(f) => {
                let tags = (f)(node);
                writeln!(
                    w,
                    "S\t{}\t{}\t{}",
                    node.node_id,
                    node.sequence().to_dna_string(),
                    tags
                )?;
            }
            _ => writeln!(
                w,
                "S\t{}\t{}",
                node.node_id,
                node.sequence().to_dna_string()
            )?,
        }

        for (_, target, dir, _) in node.l_edges() {
            if target >= node.node_id {
                let to_dir = match dir {
                    Dir::Left => "+",
                    Dir::Right => "-",
                };
                writeln!(
                    w,
                    "L\t{}\t-\t{}\t{}\t{}M",
                    node.node_id,
                    target,
                    to_dir,
                    K::k() - 1
                )?;
            }
        }

        for (_, target, dir, _) in node.r_edges() {
            if target > node.node_id {
                let to_dir = match dir {
                    Dir::Left => "+",
                    Dir::Right => "-",
                };
                writeln!(
                    w,
                    "L\t{}\t+\t{}\t{}\t{}M",
                    node.node_id,
                    target,
                    to_dir,
                    K::k() - 1
                )?;
            }
        }

        Ok(())
    }

    /// Write the graph to GFA format
    pub fn to_gfa<P: AsRef<Path>>(&self, gfa_out: P) -> Result<(), Error> {
        let wtr = BufWriter::with_capacity(BUF, File::create(gfa_out).expect("error creating gfa file"));
        self.write_gfa(&mut std::io::BufWriter::new(wtr))
    }

    pub fn write_gfa(&self, wtr: &mut impl Write) -> Result<(), Error> {
        writeln!(wtr, "H\tVN:Z:debruijn-rs")?;

        type DummyFn<K, D> = fn(&Node<'_, K, D>) -> String;
        let dummy_opt: Option<&DummyFn<K, D>> = None;

        let pb = ProgressBar::new(self.len() as u64);
        pb.set_style(ProgressStyle::with_template(PROGRESS_STYLE).unwrap().progress_chars("#/-"));
        pb.set_message(format!("{:<32}", "writing graph to GFA file"));

        for i in (0..self.len()).progress_with(pb) {
            let n = self.get_node(i);
            self.node_to_gfa(&n, wtr, dummy_opt)?;
        }

        wtr.flush().unwrap();

        Ok(())
    }

    /// Write the graph to GFA format
    pub fn to_gfa_with_tags<P: AsRef<Path>, F: Fn(&Node<'_, K, D>) -> String>(
        &self,
        gfa_out: P,
        tag_func: F,
    ) -> Result<(), Error> {
        let mut wtr = BufWriter::with_capacity(BUF, File::create(gfa_out).expect("error creatinf gfa file"));
        writeln!(wtr, "H\tVN:Z:debruijn-rs")?;

        let pb = ProgressBar::new(self.len() as u64);
        pb.set_style(ProgressStyle::with_template(PROGRESS_STYLE).unwrap().progress_chars("#/-"));
        pb.set_message(format!("{:<32}", "writing graph to GFA file"));

        for i in (0..self.len()).progress_with(pb) {
            let n = self.get_node(i);
            self.node_to_gfa(&n, &mut wtr, Some(&tag_func))?;
        }

        wtr.flush().unwrap();

        Ok(())
    }

    /// Write the graph to GFA format, with multithreading, 
    /// pass `tag_func=None` to write without tags
    pub fn to_gfa_otags_parallel<P: AsRef<Path> + Display + Sync, F: Fn(&Node<'_, K, D>) -> String + Sync>(
        &self,
        gfa_out: P,
        tag_func: Option<&F>,
    ) -> Result<(), Error> 
    where 
    K: Sync,
    D: Sync,
    {
        // split into ranges according to thread count
        let slices = current_num_threads();
        let n_nodes = self.len();
        let sz = n_nodes / slices + 1;

        debug!("n_nodes: {}", n_nodes);
        debug!("sz: {}", sz);

        let mut parallel_ranges = Vec::with_capacity(slices);
        let mut start = 0;
        while start < n_nodes {
            parallel_ranges.push(start..start + sz);
            start += sz;
        }

        let last_start = parallel_ranges.pop().expect("no kmers in parallel ranges").start;
        parallel_ranges.push(last_start..n_nodes);
        debug!("parallel ranges: {:?}", parallel_ranges);

        let mut files = Vec::with_capacity(current_num_threads());

        for i in 0..parallel_ranges.len() {
            files.push(format!("{}-{}.gfa", gfa_out, i));
        } 

        let pb = ProgressBar::new(self.len() as u64);
        pb.set_style(ProgressStyle::with_template(PROGRESS_STYLE).unwrap().progress_chars("#/-"));
        pb.set_message(format!("{:<32}", "writing graph to GFA file"));
        
        
        parallel_ranges.into_par_iter().enumerate().for_each(|(i, range)| {
            let mut wtr = BufWriter::with_capacity(BUF, File::create(&files[i]).expect("error creating parallel gfa file"));

            for i in range {
                let n = self.get_node(i);
                self.node_to_gfa(&n, &mut wtr, tag_func).unwrap();
                pb.inc(1);
            }

            wtr.flush().unwrap();
        });

        pb.finish_and_clear();

        // combine files
        let mut out_file = BufWriter::with_capacity(BUF, File::create(format!("{}.gfa", gfa_out)).expect("error creating combined gfa file"));
        writeln!(out_file, "H\tVN:Z:debruijn-rs")?;

        let pb = ProgressBar::new(files.len() as u64);
        pb.set_style(ProgressStyle::with_template(PROGRESS_STYLE).unwrap().progress_chars("#/-"));
        pb.set_message(format!("{:<32}", "combining files"));

        for file in files.iter() {
            let open_file = File::open(file).expect("error opening parallel gfa file");
            let mut reader = BufReader::new(open_file);
            let mut buffer = [0; BUF];

            loop {
                let linecount = reader.read(&mut buffer).unwrap();
                if linecount == 0 { break }
                out_file.write_all(&buffer[..linecount]).unwrap();
            }

            remove_file(file).unwrap();
        }

        out_file.flush().unwrap();

        Ok(())
    }

    /// Write the graph to GFA format
    pub fn to_gfa_partial<P: AsRef<Path>, F: Fn(&Node<'_, K, D>) -> String>(&self, gfa_out: P, tag_func: Option<&F>, nodes: Vec<usize>) -> Result<(), Error> {
        let mut wtr = BufWriter::with_capacity(BUF, File::create(gfa_out).expect("error creating gfa file"));
        writeln!(wtr, "H\tVN:Z:debruijn-rs")?;

        let pb = ProgressBar::new(self.len() as u64);
        pb.set_style(ProgressStyle::with_template(PROGRESS_STYLE).unwrap().progress_chars("#/-"));
        pb.set_message(format!("{:<32}", "writing graph to GFA file"));

        for i in nodes.into_iter().progress_with(pb) {
            let n = self.get_node(i);
            self.node_to_gfa(&n, &mut wtr, tag_func)?;
        }

        wtr.flush().unwrap();

        Ok(())    
    }

    pub fn to_json_rest<W: Write, F: Fn(&D) -> Value>(
        &self,
        fmt_func: F,
        mut writer: &mut W,
        rest: Option<Value>,
    ) {
        writeln!(writer, "{{\n\"nodes\": [").unwrap();
        for i in 0..self.len() {
            let node = self.get_node(i);
            node.to_json(&fmt_func, writer);
            if i == self.len() - 1 {
                writeln!(writer).unwrap();
            } else {
                writeln!(writer, ",").unwrap();
            }
        }
        writeln!(writer, "],").unwrap();

        writeln!(writer, "\"links\": [").unwrap();
        for i in 0..self.len() {
            let node = self.get_node(i);
            match node.edges_to_json(writer) {
                true => {
                    if i == self.len() - 1 {
                        writeln!(writer).unwrap();
                    } else {
                        writeln!(writer, ",").unwrap();
                    }
                }
                _ => continue,
            }
        }
        writeln!(writer, "]").unwrap();

        match rest {
            Some(Value::Object(v)) => {
                for (k, v) in v.iter() {
                    writeln!(writer, ",").expect("io error");
                    write!(writer, "\"{}\": ", k).expect("io error");
                    serde_json::to_writer(&mut writer, v).expect("io error");
                    writeln!(writer).expect("io error");
                }
            }
            _ => {
                writeln!(writer).expect("io error");
            }
        }

        writeln!(writer, "}}").expect("io error");
    }

    /// Write the graph to JSON
    pub fn to_json<W: Write, F: Fn(&D) -> Value, RF: Fn(&mut W)>(
        &self,
        fmt_func: F,
        writer: &mut W,
    ) {
        self.to_json_rest(fmt_func, writer, None);
    }

    /// Print a text representation of the graph.
    pub fn print(&self) {
        println!("DebruijnGraph {{ len: {}, K: {} }} :", self.len(), K::k());
        for node in self.iter_nodes() {
            println!("{:?}", node);
        }
    }

    pub fn print_with_data(&self) {
        println!("DebruijnGraph {{ len: {}, K: {} }} :", self.len(), K::k());
        for node in self.iter_nodes() {
            println!("{:?} ({:?})", node, node.data());
        }
    }

    pub fn max_path_beam<F, F2>(&self, beam: usize, score: F, _solid_path: F2) -> Vec<(usize, Dir)>
    where
        F: Fn(&D) -> f32,
        F2: Fn(&D) -> bool,
    {
        if self.is_empty() {
            return Vec::default();
        }

        let mut states = Vec::new();

        for i in 0..self.len() {
            let node = self.get_node(i);

            // Initialize beam search on terminal nodes
            if node.exts().num_exts_l() == 0 || node.exts().num_exts_r() == 0 {
                let dir = if node.exts().num_exts_l() > 0 {
                    Dir::Right
                } else {
                    Dir::Left
                };

                let status = if node.exts().num_exts_l() == 0 && node.exts().num_exts_r() == 0 {
                    Status::End
                } else {
                    Status::Active
                };

                let mut path = SmallVec8::new();
                path.push((i as u32, dir));

                let s = State {
                    path,
                    status,
                    score: score(node.data()),
                };
                states.push(s);
            }
        }

        // No end nodes -- just start on the first node
        if states.is_empty() {
            // Make a start
            let node = self.get_node(0);
            let mut path = SmallVec8::new();
            path.push((0, Dir::Left));
            states.push(State {
                path,
                status: Status::Active,
                score: score(node.data()),
            });
        }

        // Beam search until we can't find any more expansions
        let mut active = true;
        while active {
            let mut new_states = Vec::with_capacity(states.len());
            active = false;

            for s in states {
                if s.status == Status::Active {
                    active = true;
                    let expanded = self.expand_state(&s, &score);
                    new_states.extend(expanded);
                } else {
                    new_states.push(s)
                }
            }

            // workaround to sort by descending score - will panic if there are NaN scores
            new_states.sort_by(|a, b| (-(a.score)).partial_cmp(&-(b.score)).unwrap());
            new_states.truncate(beam);
            states = new_states;
        }

        for (i, state) in states.iter().take(5).enumerate() {
            trace!("i:{}  -- {:?}", i, state);
        }

        // convert back to using usize for node_id
        states[0]
            .path
            .iter()
            .map(|&(node, dir)| (node as usize, dir))
            .collect()
    }

    fn expand_state<F>(&self, state: &State, score: &F) -> SmallVec4<State>
    where
        F: Fn(&D) -> f32,
    {
        if state.status != Status::Active {
            panic!("only attempt to expand active states")
        }

        let (node_id, dir) = state.path[state.path.len() - 1];
        let node = self.get_node(node_id as usize);
        let mut new_states = SmallVec4::new();

        for (_, next_node_id, incoming_dir, _) in node.edges(dir.flip()) {
            let next_node = self.get_node(next_node_id);
            let new_score = state.score + score(next_node.data());

            let cycle = state
                .path
                .iter()
                .any(|&(prev_node, _)| prev_node == (next_node_id as u32));

            let status = if cycle {
                Status::Cycle
            } else if next_node.edges(incoming_dir.flip()).is_empty() {
                Status::End
            } else {
                Status::Active
            };

            let mut new_path = state.path.clone();
            new_path.push((next_node_id as u32, incoming_dir));

            let next_state = State {
                path: new_path,
                score: new_score,
                status,
            };

            new_states.push(next_state);
        }

        new_states
    }


    pub fn iter_components(&self) -> IterComponents<K, D> {
        let mut visited: Vec<bool> = Vec::with_capacity(self.len());
        let pos = 0;

        for _i in 0..self.len() {
            visited.push(false);
        }

        IterComponents { 
            graph: self, 
            visited, 
            pos }
    }


    /// iteratively returns 2D Vec with node_ids grouped according to the connected components they form
    pub fn components_i(&self) -> Vec<Vec<usize>> {
        let mut components: Vec<Vec<usize>> = Vec::with_capacity(self.len());
        let mut visited: Vec<bool> = Vec::with_capacity(self.len());

        for _i in 0..self.len() {
            visited.push(false);
        }

        for i in 0..self.len() {
            if !visited[i] {
                let comp = self.component_i(&mut visited, i);
                components.push(comp);
            }
        }

        components
    }

    /// recursively detects which nodes form separate graph components
    /// returns 2D vector with node ids per component
    /// (may lead to stack overflow)
    pub fn components_r(&self) -> Vec<Vec<usize>> {
        let mut components: Vec<Vec<usize>> = Vec::with_capacity(self.len());
        let mut visited: Vec<bool> = Vec::with_capacity(self.len());

        for _i in 0..self.len() {
            visited.push(false);
        }

        for i in 0..self.len() {
            if !visited[i] {
                let mut comp: Vec<usize> = Vec::new();
                self.component_r(&mut visited, i, &mut comp);
                components.push(comp);
            }
        }

        components

    }

    fn component_r<'a>(&'a self, visited: &'a mut Vec<bool>, i: usize, comp: &'a mut Vec<usize>) {
        
        visited[i] = true;
        comp.push(i);
        let mut edges = self.find_edges(i, Dir::Left);
        let mut r_edges = self.find_edges(i, Dir::Right);

        edges.append(&mut r_edges);

        for (_, edge, _, _) in edges.iter() {
            if !visited[*edge] {
                self.component_r(visited, *edge, comp);
            }
        }
    }

    fn component_i<'a>(&'a self, visited: &'a mut [bool], i: usize) -> Vec<usize> {
        let mut edges: Vec<usize> = Vec::new();
        let mut comp: Vec<usize> = Vec::new();

        edges.push(i);

        while let Some(current_edge) = edges.pop() {
            if !visited[current_edge] { 
                comp.push(current_edge);
                visited[current_edge] = true;

                let mut l_edges = self.find_edges(current_edge, Dir::Left);
                let mut r_edges = self.find_edges(current_edge, Dir::Right);

                l_edges.append(&mut r_edges);

                for (_, new_edge, _, _) in l_edges.into_iter() {
                    if !visited[new_edge] {
                        edges.push(new_edge);
                    }
                }
            }
        }
        comp
    }

    pub fn find_bad_nodes<F: Fn(&Node<'_, K, D>) -> bool>(&self, valid: F) -> Vec<usize> {
        let mut bad_nodes = Vec::new();

        for (i, node) in enumerate(self.iter_nodes()) {
            if !valid(&node) { bad_nodes.push(i); }
        }
        
        bad_nodes
    }
}

impl<K: Kmer, SD: Debug> DebruijnGraph<K, SD> {
    pub fn create_colors<'a, 'b: 'a, DI>(&'a self, config: &SummaryConfig, color_mode: ColorMode<'b>) -> Colors<'b, SD, DI> 
    where 
    SD: SummaryData<DI>,
    {
        Colors::new(self, config, color_mode)
    }
    
    /// edge mults will contain hanging edges if the nodes were filtered
    pub fn fix_edge_mults<DI>(&mut self) 
    where 
        SD: SummaryData<DI>
    {
        if self.get_node(0).data().edge_mults().is_some() {
            for i in 0..self.len() {
                self.base.data[i].fix_edge_mults(self.base.exts[i]);
            }
        }
    }
}

#[derive(Debug, Eq, PartialEq)]
enum Status {
    Active,
    End,
    Cycle,
}

#[derive(Debug)]
struct State {
    path: SmallVec8<(u32, Dir)>,
    score: f32,
    status: Status,
}

impl State {}

/// Iterator over nodes in a `DeBruijnGraph`
pub struct NodeIter<'a, K: Kmer + 'a, D: Debug + 'a> {
    graph: &'a DebruijnGraph<K, D>,
    node_id: usize,
}

impl<'a, K: Kmer + 'a, D: Debug + 'a> Iterator for NodeIter<'a, K, D> {
    type Item = Node<'a, K, D>;

    fn next(&mut self) -> Option<Node<'a, K, D>> {
        if self.node_id < self.graph.len() {
            let node = self.graph.get_node(self.node_id);
            self.node_id += 1;
            Some(node)
        } else {
            None
        }
    }
}

impl<'a, K: Kmer + 'a, D: Debug + 'a> IntoIterator for &'a DebruijnGraph<K, D> {
    type Item = NodeKmer<'a, K, D>;
    type IntoIter = NodeIntoIter<'a, K, D>;

    fn into_iter(self) -> Self::IntoIter {
        NodeIntoIter {
            graph: self,
            node_id: 0,
        }
    }
}

/// Iterator over nodes in a `DeBruijnGraph`
pub struct NodeIntoIter<'a, K: Kmer + 'a, D: Debug + 'a> {
    graph: &'a DebruijnGraph<K, D>,
    node_id: usize,
}

impl<'a, K: Kmer + 'a, D: Debug + 'a> Iterator for NodeIntoIter<'a, K, D> {
    type Item = NodeKmer<'a, K, D>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.node_id < self.graph.len() {
            let node_id = self.node_id;
            let node = self.graph.get_node(node_id);
            let node_seq = node.sequence();

            self.node_id += 1;
            Some(NodeKmer {
                node_id,
                node_seq_slice: node_seq,
                phantom_d: PhantomData,
                phantom_k: PhantomData,
            })
        } else {
            None
        }
    }
}

/// A `DebruijnGraph` node with a reference to the sequence of the node.
#[derive(Clone)]
pub struct NodeKmer<'a, K: Kmer + 'a, D: Debug + 'a> {
    pub node_id: usize,
    node_seq_slice: DnaStringSlice<'a>,
    phantom_k: PhantomData<K>,
    phantom_d: PhantomData<D>,
}

/// An iterator over the kmers in a `DeBruijn graph node`
pub struct NodeKmerIter<'a, K: Kmer + 'a, D: Debug + 'a> {
    kmer_id: usize,
    kmer: K,
    num_kmers: usize,
    node_seq_slice: DnaStringSlice<'a>,
    phantom_k: PhantomData<K>,
    phantom_d: PhantomData<D>,
}

impl<'a, K: Kmer + 'a, D: Debug + 'a> IntoIterator for NodeKmer<'a, K, D> {
    type Item = K;
    type IntoIter = NodeKmerIter<'a, K, D>;

    fn into_iter(self) -> Self::IntoIter {
        let num_kmers = self.node_seq_slice.len() - K::k() + 1;

        let kmer = if num_kmers > 0 {
            self.node_seq_slice.get_kmer::<K>(0)
        } else {
            K::empty()
        };

        NodeKmerIter {
            kmer_id: 0,
            kmer,
            num_kmers,
            node_seq_slice: self.node_seq_slice,
            phantom_k: PhantomData,
            phantom_d: PhantomData,
        }
    }
}

impl<'a, K: Kmer + 'a, D: Debug + 'a> Iterator for NodeKmerIter<'a, K, D> {
    type Item = K;

    fn next(&mut self) -> Option<Self::Item> {
        if self.num_kmers == self.kmer_id {
            None
        } else {
            let current_kmer = self.kmer;

            self.kmer_id += 1;
            if self.kmer_id < self.num_kmers {
                let next_base = self.node_seq_slice.get(self.kmer_id + K::k() - 1);
                let new_kmer = self.kmer.extend_right(next_base);
                self.kmer = new_kmer;
            }

            Some(current_kmer)
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.num_kmers, Some(self.num_kmers))
    }

    /// Provide a 'fast-forward' capability for this iterator
    /// MPHF will use this to reduce the number of kmers that
    /// need to be produced.
    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        if n <= 4 {
            // for small skips forward, shift one base at a time
            for _ in 0..n {
                self.next();
            }
        } else {
            self.kmer_id += n;
            self.kmer = self.node_seq_slice.get_kmer::<K>(self.kmer_id);
        }

        self.next()
    }
}

/// Marker signifying that NodeKmerIter has a known size.
impl<'a, K: Kmer + 'a, D: Debug + 'a> ExactSizeIterator for NodeKmerIter<'a, K, D> {}

/// Unbranched sequence in the DeBruijn graph
pub struct Node<'a, K: Kmer + 'a, D: 'a> {
    pub node_id: usize,
    pub graph: &'a DebruijnGraph<K, D>,
}

impl<'a, K: Kmer, D: Debug> Node<'a, K, D> {
    /// Length of the sequence of this node
    pub fn len(&self) -> usize {
        self.graph.base.sequences.get(self.node_id).len()
    }

    pub fn is_empty(&self) -> bool {
        self.graph.base.sequences.get(self.node_id).is_empty()
    }

    /// Sequence of the node
    pub fn sequence(&self) -> DnaStringSlice<'a> {
        self.graph.base.sequences.get(self.node_id)
    }

    /// Reference to auxiliarly data associated with the node
    pub fn data(&self) -> &'a D {
        &self.graph.base.data[self.node_id]
    }

    /// Extension bases from this node
    pub fn exts(&self) -> Exts {
        self.graph.base.exts[self.node_id]
    }

    /// Edges leaving the left side of the node in the format
    /// (base, target_node id, incoming side of target node, whether target node is flipped)
    pub fn l_edges(&self) -> SmallVec4<(u8, usize, Dir, bool)> {
        self.graph.find_edges(self.node_id, Dir::Left)
    }

    /// Edges leaving the right side of the node in the format
    /// (base, target_node id, incoming side of target node, whether target node is flipped)
    pub fn r_edges(&self) -> SmallVec4<(u8, usize, Dir, bool)> {
        self.graph.find_edges(self.node_id, Dir::Right)
    }

    /// Edges leaving the 'dir' side of the node in the format
    /// (base, target_node id, incoming side of target node, whether target node is flipped)
    pub fn edges(&self, dir: Dir) -> SmallVec4<(u8, usize, Dir, bool)> {
        self.graph.find_edges(self.node_id, dir)
    }

    fn to_json<F: Fn(&D) -> Value>(&self, func: &F, f: &mut dyn Write) {
        write!(
            f,
            "{{\"id\":\"{}\",\"L\":{},\"D\":{},\"Se\":\"{:?}\"}}",
            self.node_id,
            self.sequence().len(),
            (func)(self.data()),
            self.sequence(),
        )
        .unwrap();
    }

    fn edges_to_json(&self, f: &mut dyn Write) -> bool {
        let mut wrote = false;
        let edges = self.r_edges();
        for (idx, &(_, id, incoming_dir, _)) in edges.iter().enumerate() {
            write!(
                f,
                "{{\"source\":\"{}\",\"target\":\"{}\",\"D\":\"{}\"}}",
                self.node_id,
                id,
                match incoming_dir {
                    Dir::Left => "L",
                    Dir::Right => "R",
                }
            )
            .unwrap();

            if idx < edges.len() - 1 {
                write!(f, ",").unwrap();
            }

            wrote = true;
        }
        wrote
    }
}

// TODO make generic instead of u8 (u8 is sufficient for dbg)
impl<K: Kmer, SD: Debug> Node<'_, K, SD>  {
    /// get default format for dot edges based on node data
    pub fn edge_dot_default<DI>(&self, colors: &Colors<SD, DI>, base: u8, incoming_dir: Dir, flipped: bool) -> String 
    where SD: SummaryData<DI>
    {
        // set color based on dir
        let color = match incoming_dir {
            Dir::Left => "blue",
            Dir::Right => "red"
        };
        
        if let Some(em) = self.data().edge_mults() {
            
            let dir = if flipped { 
                incoming_dir 
            } else {
                incoming_dir.flip()
            };

            // set penwidth based on count
            let count = em.edge_mult(base, dir);
            let penwidth = colors.edge_width(count);

            format!("[color={color}, penwidth={penwidth}, label=\"{}: {count}\"]", bits_to_base(base))
        } else {
            format!("[color={color}]")
        }
    }

    /// get default format for dot nodes, based on node data
    pub fn node_dot_default<DI>(&self, colors: &Colors<SD, DI>, config: &SummaryConfig, tag_translator: &bimap::BiHashMap<String, DI> , outline: bool) -> String
    where SD: SummaryData<DI>
    {
        // set color based on labels/fold change/p-value
        let color = colors.node_color(self.data(), config, outline);

        let data_info = self.data().print(tag_translator, config);
        const MIN_TEXT_WIDTH: usize = 40;
        let wrap = if self.len() > MIN_TEXT_WIDTH { self.len() } else { MIN_TEXT_WIDTH };

        let label = textwrap::fill(&format!("id: {}, len: {}, exts: {:?}, seq: {}\n{}", 
            self.node_id,
            self.len(),
            self.exts(),
            self.sequence(),
            data_info
        ), wrap);

        format!("[style=filled, {color}, label=\"{label}\"]")
    }
}

impl<K: Kmer, D> fmt::Debug for Node<'_, K, D>
where
    D: Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Node {{ id:{}, Exts: {:?}, L:{:?} R:{:?}, Seq: {:?}, Data: {:?} }}",
            self.node_id,
            self.exts(),
            self.l_edges(),
            self.r_edges(),
            self.sequence().len(),
            self.data()
        )
    }
}

pub struct IterComponents<'a, K: Kmer, D> {
    graph: &'a DebruijnGraph<K, D>,
    visited: Vec<bool>,
    pos: usize,
}

impl<K: Kmer, D: Debug> Iterator for IterComponents<'_, K, D> {
    type Item = Vec<usize>;
    fn next(&mut self) -> Option<Self::Item> {
        while self.pos < self.graph.len() {
            if !self.visited[self.pos] {
                let comp = self.graph.component_i(&mut self.visited, self.pos);
                self.pos += 1;
                return Some(comp)
            } else {
                self.pos += 1;
            }
        }
        assert!(self.visited.iter().map(|x| *x as usize).sum::<usize>() == self.graph.len());
        None
    }
    
}

pub struct PathCompIter<'a, K: Kmer, D: Debug, F, F2> 
where 
F: Fn(&D) -> f32,
F2: Fn(&D) -> bool
{
    graph: &'a DebruijnGraph<K, D>,
    component_iterator: IterComponents<'a, K, D>,
    graph_pos: usize,
    score: F,
    solid_path: F2,
}

/// returns the component and the "best" path in the component
impl<K: Kmer, D: Debug, F, F2> Iterator for PathCompIter<'_, K, D, F, F2> 
where 
F: Fn(&D) -> f32,
F2: Fn(&D) -> bool
{
    type Item = (Vec<usize>, VecDeque<(usize, Dir)>,);
    fn next(&mut self) -> Option<Self::Item> {
        match self.component_iterator.next() {
            Some(component) => {
                let current_comp = component;
                
    
                let mut best_node = current_comp[0];
                let mut best_score = f32::MIN;
                for c in current_comp.iter() {
                    let node = self.graph.get_node(*c);
                    let node_score = (self.score)(node.data());
    
                    if node_score > best_score {
                        best_node = *c;
                        best_score = node_score;
                    }
                }
    
                let oscore = |state| match state {
                    None => 0.0,
                    Some((id, _)) => (self.score)(self.graph.get_node(id).data()),
                };
    
                let osolid_path = |state| match state {
                    None => false,
                    Some((id, _)) => (self.solid_path)(self.graph.get_node(id).data()),
                };
    
                // Now expand in each direction, greedily taking the best path. Stop if we hit a node we've
                // already put into the path
                let mut used_nodes = HashSet::new();
                let mut path = VecDeque::new();
    
                // Start w/ initial state
                used_nodes.insert(best_node);
                path.push_front((best_node, Dir::Left));
    
                for init in [(best_node, Dir::Left, false), (best_node, Dir::Right, true)].iter() {
                    let &(start_node, dir, do_flip) = init;
                    let mut current = (start_node, dir);
    
                    loop {
                        let mut next = None;
                        let (cur_id, incoming_dir) = current;
                        let node = self.graph.get_node(cur_id);
                        let edges = node.edges(incoming_dir.flip());
    
                        let mut solid_paths = 0;
                        for (_, id, dir, _) in edges {
                            let cand = Some((id, dir));
                            if osolid_path(cand) {
                                solid_paths += 1;

                                // second if clause is outside of first in original code (see max_path) 
                                // but would basically ignore path validity.
                                if oscore(cand) > oscore(next) {
                                    next = cand;
                                }
                            }
    
                            if oscore(cand) > oscore(next) {
                                next = cand;
                            }
                        }
    
                        if solid_paths > 1 {
                            break;
                        }
    
                        match next {
                            Some((next_id, next_incoming)) if !used_nodes.contains(&next_id) => {
                                if do_flip {
                                    path.push_front((next_id, next_incoming.flip()));
                                } else {
                                    path.push_back((next_id, next_incoming));
                                }
    
                                used_nodes.insert(next_id);
                                current = (next_id, next_incoming);
                            }
                            _ => break,
                        }
                    }
                }
                
                
                Some((current_comp, path))
            }, 
            None => {
                // should technically not need graph_pos after this 
                self.graph_pos += 1;
                None
            }
        }
    }
}

#[cfg(test)]
mod test {
    use std::{fs::File, io::BufReader};

    use crate::{kmer::Kmer16, summarizer::TagsCountsSumData};

    use super::DebruijnGraph;

    #[test]
    #[cfg(not(feature = "sample128"))]
    fn test_components() {
        use crate::{summarizer::SummaryData, Dir, BUF};

        let path = "test_data/400.graph.dbg";
        let file = BufReader::with_capacity(BUF, File::open(path).unwrap());

        let (graph, _, _): (DebruijnGraph<Kmer16, TagsCountsSumData>, Vec<String>, crate::summarizer::SummaryConfig) = 
            bincode::deserialize_from(file).expect("error deserializing graph");

        let components = graph.iter_components();

        let check_components = [
            vec![3, 67, 130, 133, 59, 119, 97, 110, 68, 137, 29, 84, 131, 43, 30, 91, 14, 70, 79, 142, 136, 105, 103, 62, 
                141, 104, 134, 88, 38, 81, 108, 92, 135, 96, 116, 121, 63, 124, 106, 129, 132, 126, 93, 109, 83, 112, 118, 
                123, 125, 78, 122, 115, 75, 128, 140, 111, 26, 143, 113],
            vec![41, 138, 100, 139, 86],
            vec![53, 117, 127],
            vec![69, 144, 77, 120, 114, 107, 101],
        ];

        let mut counter = 0;

        for component in components {
            if component.len() > 1 { 
                println!("component: {:?}", component);
                assert_eq!(component, check_components[counter]);
                counter += 1;
            }
        }
        assert_eq!(vec![(139, Dir::Left)], graph.max_path(|data| data.score(), |_| true));
    }
}

