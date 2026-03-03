// Copyright 2017 10x Genomics

//! Containers for path-compressed De Bruijn graphs

use bimap::BiHashMap;
use bio::io::fasta;
use bit_set::BitSet;
use indicatif::ProgressBar;
use indicatif::ProgressIterator;
use indicatif::ProgressStyle;
use itertools::enumerate;
use log::warn;
use log::{debug, trace};
use rayon::prelude::*;
use rayon::current_num_threads;
use serde_derive::{Deserialize, Serialize};
use smallvec::SmallVec;
use std::borrow::Borrow;

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

use crate::BaseQuality;
use crate::bits_to_base;
use crate::colors::ColorMode;
use crate::colors::Colors;
use crate::compression::CompressionSpec;
use crate::dna_string::{DnaString, DnaStringSlice, PackedDnaStringSet};
use crate::summarizer::SummaryConfig;
use crate::summarizer::SummaryData;
use crate::summarizer::Translator;
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
    pub fn get_node(&'_ self, node_id: usize) -> Node<'_, K, D> {
        Node {
            node_id,
            graph: self,
        }
    }

    /// Get a node given it's `node_id`
    pub fn get_node_kmer(&'_ self, node_id: usize) -> NodeKmer<'_, K, D> {
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
    pub fn iter_nodes(&'_ self) -> NodeIter<'_, K, D> {
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

    /// mutable reference to the auxiliary data of the node node_id
    pub fn mut_data(&mut self, node_id: usize) -> &mut D {
        &mut self.base.data[node_id]
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

    pub fn iter_max_path_comp<F, F2>(&'_ self, score: F, solid_path: F2) -> PathCompIter<'_, K, D, F, F2> 
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

    /// map sequences from a fasta file to a **completely uncompressed** and **stranded** debruijn graph
    pub fn map_transcripts<P>(&self, path: P, translator: &mut Translator) -> Result<Vec<Box<[ID]>>, String> 
    where 
        P: AsRef<Path>
    {
        // return err if not stranded
        if !self.base.stranded { return Err("graph has to be stranded".to_string()) };

        let reader = fasta::Reader::new(BufReader::new(File::open(path).unwrap()));
        let mut node_transcript_ids = vec![Vec::new(); self.len()];

        let mut backup_id_tr = BiHashMap::new();

        let id_tr = if let Some(id_tr) = translator.mut_id_translator() {
            id_tr
        } else {
            &mut backup_id_tr
        };

        // go through each transcript and map to graph
        for result in reader.records() {
            let record = result.expect("error parsing transcripts fasta");

            // get gene id or make new id
            let gene_id = match id_tr.get_by_left(&record.id().to_string()) {
                Some(id) => *id,
                None => {
                    let new_id = id_tr.len() as ID;
                    id_tr.insert(record.id().to_string(), new_id);
                    new_id
                }
            };

            // iterate over k-mers in transcript and find each one in the graph
            let sequence = DnaString::from_acgt_bytes(record.seq());
            for kmer in sequence.iter_kmers::<K>() {
                if let Some(node) = self.search_kmer(kmer, Dir::Right) {
                    node_transcript_ids[node].push(gene_id);
                }
            }
        }

        // change to boxed slices to reduce memory
        let boxed_transcripts = node_transcript_ids.into_iter().map(|vec| vec.into()).collect();

        Ok(boxed_transcripts)
    }

    /// write a node to a dot file
    /// 
    /// ### Arguments: 
    /// 
    /// * `node`: [`Node<K, D>`] which will be written to a dot file
    /// * `node_label`: closure taking [`Node<K, D>`] and returning a string containing commands for dot nodes 
    /// * `edge_label`: closure taking [`Node<K, D>`], the base as a [`u8`], the incoming [`Dir`] of the edge 
    ///   and if the neighbor is flipped - returns a string containing commands for dot edges, 
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
    /// * `node_label`: closure taking [`Node<K, D>`] and returning a string containing commands for dot nodes, e.g. [`Node::node_dot_default`]
    /// * `edge_label`: closure taking [`Node<K, D>`], the base as a [`u8`], the incoming [`Dir`] of the edge, e.g. [`Node::edge_dot_default`]
    ///   and if the neighbor is flipped - returns a string containing commands for dot edges, 
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

        writeln!(&mut f, "digraph {{\nrankdir=\"LR\"\nmodel=subset").unwrap();
        for i in (0..self.len()).progress_with(pb) {
            self.node_to_dot(&self.get_node(i), node_label, edge_label, &mut f);
        }
        writeln!(&mut f, "}}").unwrap();
        
        f.flush().unwrap();
        debug!("large to dot loop: {}", self.len());
    }

    /// Write the graph to a dot file, highlight the nodes which form the 
    /// "best" path, according to [`PathCompIter`], with the number of occurences 
    /// as the score and `solid_path` always `true`.
    /// The nodes are formatted according to [`Node::node_dot_default`].
    /// 
    /// ### Arguments: 
    /// 
    /// * `path`: path to the output file
    /// * `edge_label`: closure taking [`Node<K, D>`], the base as a [`u8`], the incoming [`Dir`] of the edge, e.g. [`Node::edge_dot_default`]
    ///   and if the neighbor is flipped - returns a string containing commands for dot edges, 
    /// * `colors`: a [`Colors`] with the color settings for the graph
    /// * `translator`: a [`Translator`] which translates tags or IDs to strings
    /// * `config`: a [`SummaryConfig`] which contains settings for the graph
    pub fn to_dot_with_path<P, FE, DI>(&self, path: P, edge_label: &FE, colors: &Colors<'_, D, DI>, translator: &Translator, config: &SummaryConfig, translate_id_groups: bool)
    where 
    P: AsRef<Path>,
    D: SummaryData<DI>,
    FE: Fn(&Node<K, D>, u8, Dir, bool) -> String,
    {
        let mut f = BufWriter::with_capacity(BUF, File::create(path).expect("error creating dot file"));

        writeln!(&mut f, "digraph {{\nrankdir=\"LR\"\nmodel=subset").unwrap();

        // iterate over components
        for (component, path) in self.iter_max_path_comp(|d| d.sum().unwrap_or(1) as f32, |_| true) {
            let hashed_path = path.into_iter().map(|(id, _)| id).collect::<HashSet<usize>>();
            for node_id in component {
                self.node_to_dot(
                    &self.get_node(node_id),
                    &|node| node.node_dot_default(colors, config, translator, hashed_path.contains(&node_id), translate_id_groups), 
                    edge_label, 
                    &mut f
                );
            }
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
    ///   and if the neighbor is flipped - returns a string containing commands for dot edges, 
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

        writeln!(&mut out_file, "digraph {{\nrankdir=\"LR\"\nmodel=subset").unwrap();

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
    ///   and if the neighbor is flipped - returns a string containing commands for dot edges, 
    /// * `nodes`: [`Vec<usize>`] listing all IDs of nodes which should be included
    pub fn to_dot_partial<P, FN, FE>(&self, path: P, node_label: &FN, edge_label: &FE, nodes: &[usize]) 
    where 
        P: AsRef<Path>,
        FN: Fn(&Node<K, D>) -> String,
        FE: Fn(&Node<K, D>, u8, Dir, bool) -> String,
    {
        let mut f = BufWriter::with_capacity(BUF, File::create(path).expect("error creating dot file"));

        let pb = ProgressBar::new(nodes.len() as u64);
        pb.set_style(ProgressStyle::with_template(PROGRESS_STYLE).unwrap().progress_chars("#/-"));
        pb.set_message(format!("{:<32}", "writing graph to DOT file"));

        writeln!(&mut f, "digraph {{\nrankdir=\"LR\"\nmodel=subset").unwrap();
        for i in nodes.iter().progress_with(pb) {
            self.node_to_dot(&self.get_node(*i), node_label, edge_label, &mut f);
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

    fn node_to_tsv<W: Write, F>(&self, writer: &mut W, node_id: usize, data_format: F) -> Result<(), Box<dyn std::error::Error>> 
    where 
        F: Fn(&Node<'_, K, D>) -> String,
    {

        let node = self.get_node(node_id);
        let l_e = node.l_edges();
        let r_e = node.r_edges();

        
        if self.base.stranded {
            // format: node id    l nb    r nb    seq    data
            let l_nb = l_e.iter().map(|(_b, nb, _d, _f)| *nb).collect::<Vec<_>>();
            let r_nb = r_e.iter().map(|(_b, nb, _d, _f)| *nb).collect::<Vec<_>>();
            writeln!(writer, "{node_id}\t{:?}\t{:?}\t{}\t{}", l_nb, r_nb, node.sequence(), data_format(&node))?
        } else {
            // format: node id    l nb    l inc dir    r nb    r inc dir    seq    data
            let l_nb = l_e.iter().map(|(_b, nb, _d, _f)| *nb).collect::<Vec<_>>();
            let r_nb = r_e.iter().map(|(_b, nb, _d, _f)| *nb).collect::<Vec<_>>();
            let l_nb_dirs = l_e.iter().map(|(_b, _nb, dir, _f)| *dir).collect::<Vec<_>>();
            let r_nb_dirs = r_e.iter().map(|(_b, _nb, dir, _f)| *dir).collect::<Vec<_>>();
            writeln!(writer, "{node_id}\t{:?}\t{:?}\t{:?}\t{:?}\t{}\t{}", l_nb, l_nb_dirs, r_nb, r_nb_dirs, node.sequence(), data_format(&node))?
        }

        Ok(())
    }

    /// save the graph as a tsv file with custom formatting for the node data
    pub fn to_tsv<P, F>(&self, path: P, data_format: F) -> Result<(), Box<dyn std::error::Error>> 
    where 
        F: Fn(&Node<'_, K, D>) -> String,
        P: AsRef<Path> + Display,
    { 
        let mut writer = BufWriter::new(File::create(path)?);

        // different format if stranded vs unstranded
        if self.base.stranded {
            writeln!(writer, "node id\tleft neighbors\tright neighbors\tsequence\tdata")?;
        } else {
            writeln!(writer, "node id\tleft neighbors\tleft nb incoming dirs\tright neighbors\tright nb incoming dirs\tsequence\tdata")?;
        }
            
        for i in 0..self.len() {
            self.node_to_tsv(&mut writer, i, &data_format)?
        }
        
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


    pub fn iter_components(&'_ self) -> IterComponents<'_, K, D> {
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

    /// iterate over all edges of the graph, item: (node, ext base, ext dir, target node)
    pub fn iter_edges(&self) -> EdgeIter<'_, K, D> {
        EdgeIter::new(self)
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
    
    /// [`crate::EdgeMult`] will contain hanging edges if the nodes were filtered
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

    /// if there are [`crate::EdgeMult`]s in the data, prune the graph by removing edges that have a low coverage
    pub fn filter_edges<DI>(&mut self, min: u32) -> Result<(), String>
    where 
        SD: SummaryData<DI>
    {
        // return if there is no edge coverage available
        if self.get_node(0).data().edge_mults().is_none() { return Err(String::from("no edge mults available")) };

        for i in 0..self.len() {
            let em = self.get_node(i).data().edge_mults().expect("shold have em").clone();
            let edges = [(Dir::Left, 0), (Dir::Left, 1), (Dir::Left, 2), (Dir::Left, 3), 
                (Dir::Right, 0), (Dir::Right, 1), (Dir::Right, 2), (Dir::Right, 3)];
            
            for (dir, base) in edges {
                if min > em.edge_mult(base, dir) {
                    // remove invalid ext from node
                    let ext = self.base.exts[i].remove(dir, base);
                    self.base.exts[i] = ext;
                }
            }
        }

        Ok(())
    }

    /// if a node has a connection to a high quality node and low quality nodes 
    /// in the same direction, remove the connections to the low quality ndoes
    pub fn remove_lq_splits<DI>(&mut self, min_quality: BaseQuality) -> Result<(), String>
    where 
        SD: SummaryData<DI>
    {
        // check we do indeed have quality and graph is stranded
        if self.get_node(0).data().quality().is_none() { return Err(String::from("no quality scores available")); }
        if !self.base.stranded { return Err(String::from("graph must be stranded to remove ladders")) };

        // iterate over all nodes and look in both directions if there is a low quality node splitting off
        for (node_id, out_dir) in (0..self.len()).flat_map(|id| [(id, Dir::Right), (id, Dir::Left)]) {
            let out_edges = self.get_node(node_id).edges(out_dir);

            // check for split
            if out_edges.len() < 2 {
                // no split, move on
                continue;
            }

            // get qualities 
            let nb_qualities = out_edges
                .iter()
                .map(|(_, nb_id, _, _)| (*nb_id, self
                    .get_node(*nb_id)
                    .data()
                    .quality()
                    .unwrap()
                )).collect::<Vec<_>>();

            // check for good neighbor
            let has_good_nb = nb_qualities.iter()
                .any(|(_, quality)| *quality >= min_quality);

            if !has_good_nb {
                // split but no good path, move on
                continue;
            }

            // remove split with low quality
            for (nb_id, _) in nb_qualities.iter().filter(|(_, quality)| *quality < min_quality ) {
                let path = vec![node_id, *nb_id];
                if self.remove_path(path, out_dir).is_err() {
                    warn!("lq tip path could not be removed")
                }
            }
        }

        Ok(())
    }

    /// remove bubbles/ladders and tips in which one path has a quality lower than the given `min_quality`
    pub fn remove_lq_ladders_tips<DI>(&mut self, min_quality: BaseQuality) -> Result<(), String>
    where
        SD: SummaryData<DI>
    {
        // check we do indeed have quality and graph is stranded
        if self.get_node(0).data().quality().is_none() { return Err(String::from("no quality scores available")); }
        if !self.base.stranded { return Err(String::from("graph must be stranded to remove ladders")) };

        let min_path = 2 * K::k() - 1;
        //let max_path = 4 * K::k() - 1;

        // iterate over nodes
        for (node_id, out_dir) in (0..self.len()).flat_map(|id| [(id, Dir::Right), (id, Dir::Left)]) {
            let in_dir = out_dir.flip();

            // check if node has multile outs to the right, at least one with bad quality and one with good quality
            let node_out_edges = self.get_node(node_id).edges(out_dir);

            let good_neighbors = node_out_edges.iter()
                .map(|(_, target_id, _, _)| *target_id)
                .filter(|target_id| self.get_node(*target_id)
                    .data()
                    .quality()
                    .unwrap() >= min_quality
                ).collect::<Vec<_>>();

            let bad_neighbors = node_out_edges.iter()
                .map(|(_, target_id, _, _)| *target_id)
                .filter(|target_id| self.get_node(*target_id)
                    .data()
                    .quality()
                    .unwrap() < min_quality
                ).collect::<Vec<_>>();

            if good_neighbors.is_empty() | bad_neighbors.is_empty() { continue; }

            // follow the bad quality paths until we reach nodes with high quality again (or max search radius)
            let mut possible_paths = Vec::new();
            let mut tips = Vec::new();

            for bn in bad_neighbors {
                let mut current_node_id = bn;
                let mut path_groups = vec![vec![node_id]];
                let mut path_length = K::k() - 1;

                let mut state = LadderState::Singular;
                let mut q_state_high = false;

                loop {
                    // check if path has reached max length -> interrupt
                    /* if path_length > max_path {
                        break;
                    } */

                    let current_node = self.get_node(current_node_id);
                    let out_edges = current_node.edges(out_dir);

                    // add current node length to path
                    path_length += current_node.len() - K::k() + 1;
                    
                    // check state and add current node to path
                    let path_index = path_groups.len() - 1;

                    // if in singular state now or before, add node to path
                    if matches!(state, LadderState::Singular) { 
                        path_groups[path_index].push(current_node_id);
                    }

                    // check if we increase ladder state
                    let q_increase = (current_node.data().quality().unwrap() >= min_quality) & !q_state_high;
                    let mult_increase = current_node.edges(in_dir).len() > 1;
                    
                    if q_increase {
                        q_state_high = true
                    }

                    if q_increase | mult_increase {
                        match state {
                            LadderState::Singular => {
                                state = LadderState::Double;
                            }
                            LadderState::Double => () // ignore
                        };
                    }

                    // check if we decrease ladder state
                    let q_decrease = (current_node.data().quality().unwrap() < min_quality) & q_state_high;
                    let mult_decrease = out_edges.len() > 1;

                    if q_decrease {
                        q_state_high = false;
                    }

                    if q_decrease | mult_decrease {
                        match state {
                            LadderState::Singular => (), // ignore
                            LadderState::Double => {
                                state = LadderState::Singular;
                                path_groups.push(vec![current_node_id]); // start new path group
                            }
                        }
                    }

                    // check if we have met end criterium -> save path
                    let quality_req =  current_node.data().quality().unwrap() >= min_quality;
                    let len_req = path_length >= min_path;
                    // path is a simple tip -> save as tip
                    let is_tip = out_edges.is_empty() & (path_groups.len() == 1); // TODO maybe remove req 2 in future

                    if quality_req & len_req {
                        possible_paths.push(path_groups);
                        break;
                    } else if is_tip {
                        tips.push(path_groups);
                        break;
                    }

                    // we have not met the conditions and keep moving

                    // find next node
                    let next_node_id = if out_edges.len() == 1 {
                        let (_, next_node_id, _, _) = out_edges[0];
                        next_node_id
                    } else if out_edges.is_empty() {
                        // dead end which has not qualified as tip
                        break;
                    } else {
                        // choose the lowest quality path

                        let worst_neighbor = out_edges.iter()
                            .map(|(_, target_id, _, _)| (*target_id, self.get_node(*target_id)
                                .data()
                                .quality()
                                .unwrap()))
                            .min_by(|(_, q_a), (_, q_b)| q_a.cmp(q_b));

                        if let Some((worst_nb_id, _)) = worst_neighbor {
                            worst_nb_id
                        } else {
                            break;
                        }
                    };

                    current_node_id = next_node_id;
                }
            }

            let possible_targets = possible_paths.iter().map(|p| p.last().unwrap().last().unwrap()).collect::<Vec<_>>();
            let mut confirmed_targets = Vec::new();
            // follow the good quality paths until we reach a possible target node (or max search radius)
            // if we reached a target node, send bad path to be removed from graph
            for gn in good_neighbors {
                let mut path_length = K::k() - 1;
                let mut current_node_id = gn;

                loop {
                    // check if we have exceeded the search radius
                    /* if path_length > max_path {
                        break;
                    } */

                    // add current node length to path length
                    let current_node = self.get_node(current_node_id);
                    path_length += current_node.len() - K::k() + 1;

                    // check if we have found a target
                    if possible_targets.contains(&&current_node_id) {
                        confirmed_targets.push(current_node_id);
                        break;
                    }

                    // look for next node
                    let out_edges = current_node.edges(out_dir);

                    let next_node_id = if out_edges.len() == 1 {
                        let (_, next_node_id, _, _) = out_edges[0];
                        next_node_id
                    } else if out_edges.is_empty() {
                        // dead end, break
                        break;
                    } else {
                        // use node with highest quality
                        // TODO future: use read mapping?
                        let good_neighbors = out_edges.iter()
                            .map(|(_, target_id, _, _)| (*target_id, self.get_node(*target_id)
                                .data()
                                .quality()
                                .unwrap()))
                            .max_by(|(_, q_a), (_, q_b)| q_a.cmp(q_b));
                        if let Some((best_nb_id, _)) = good_neighbors {
                            best_nb_id
                        } else {
                            break;
                        }
                    };

                    current_node_id = next_node_id;
                }
            }

            // check if we have found end nodes of possible paths by following good quality paths
            // if so, remove path
            for path_group in possible_paths {
                let target = path_group.last().unwrap().last().unwrap();

                if confirmed_targets.contains(target) {
                    for path in path_group {
                        if let Err(err) = self.remove_path(path.clone(), out_dir) {
                            warn!("lq ladder partial path could not be removed, likely cause: loop, edges were already removed. parital path: {:?}\nerr: {err}", path)
                        }
                    }
                }
            }

            // remove tip paths
            for path_group in tips {
                let path = path_group.into_iter().next().expect("empty tip path found");
                if let Err(err) = self.remove_path(path.clone(), out_dir) {
                    warn!("lq tip partial path could not be removed, likely cause: loop, edges were already removed. parital path: {:?}\nerr: {err}", path)
                }
            }


        }


        Ok(())
    }

    /// remove simple ladder structures (bubbles) caused by 1-base sequencing errors from the graph
    /// 
    /// graph must contain edge mults and be stranded
    /// this function will likely leave tips on the graph, so it is recommended to run
    /// [`DebruijnGraph::remove_tips`] afterwards
    /// 
    /// ladder structure refers to bubbles where one side has been compressed into one node
    /// but the other has low compression, due to differences in coverage and thus data variance,
    /// leading to a ladder-like appearance
    pub fn remove_ladders<DI, P>(&mut self, min_diff_factor: u32, max_avg_low_cov: f32, out_path: Option<P>) -> Result<(), String> 
    where 
        SD: SummaryData<DI>,
        P: AsRef<Path>
    {
        // interrupt if we dont have edge mults
        if self.get_node(0).data().edge_mults().is_none() { return Err(String::from("no edge mults available")) };
        if !self.base.stranded { return Err(String::from("must be stranded")) };

        let mut writer = out_path.map(|path| BufWriter::new(File::create(path).expect("error creating ladder stats file")));
        if let Some(wtr) = writer.as_mut() {
            writeln!(wtr, "avg high cov,avg low cov,high truth ratio,low truth ratio").unwrap();
        }

        // iterate over nodes
        for node_id in 0..self.len() {
            // check if node has outgoing edge(s) with both high and low coverage
            let outs = self.get_node(node_id).data().edge_mults().expect("should have em").right();

            let Some((out_max_base, &out_max_cov)) = outs.iter().rev().enumerate().filter(|&(_, &c)| c  > 0).max_by(|&(_b1, &c1), &(_b2, c2)| c1.cmp(c2)) else { continue };
            let smaller_outs = outs.iter().copied().rev().enumerate().filter(|&(_b, c)| (c > 0) & (out_max_cov > c)).collect::<Vec<_>>();

            // TODO add factor back in, remove writer

            if smaller_outs.is_empty() { continue; }

            // follow path with highest coverage until target length is reached
            let Some((target_node, avg_high_cov, high_cc)) = self.follow_ladder_path_high(node_id, out_max_base as u8, out_max_cov) else { continue; };
            
            // check all small outs
            for (s_base, s_cov) in smaller_outs {
                let Some((target_paths, avg_low_cov, low_cc)) = self.follow_ladder_path_low(node_id, s_base as u8, s_cov) else { continue; };
                let other_target_node = *target_paths.last().expect("should have at least one element").last().expect("should have at least two elements");

                if target_node == other_target_node {
                    // the two paths landed on the same node -> remove all edges in the low coverage path
                    if (s_cov * min_diff_factor <= out_max_cov) & (avg_low_cov <= max_avg_low_cov) {
                        for path in target_paths.iter() {
                            if self.remove_path(path.clone(), Dir::Right).is_err() {
                                warn!("removing ladders: partial path could not be removed, likely cause: loop, edges were already removed. parital path: {:?}", path)
                            }
                        } 
                    }
                    
                    if let Some(wtr) = writer.as_mut() {
                        writeln!(wtr, "{},{},{},{}", avg_high_cov, avg_low_cov, high_cc, low_cc).unwrap();
                    }
                }
            }
        }

        Ok(())
    }

    /// follow the presumably erroneous ladder path with low coverage
    fn follow_ladder_path_low<DI>(&self, start_node_id: usize, start_ext: u8, start_cov: u32) -> Option<(Vec<Vec<usize>>, f32, f32)> 
    where SD: SummaryData<DI>
    {
        // state is switched if coverage rises by 20% + 2 (so it's at least 2 more)
        // nodes with higher state are kept connected together
        // tests have shown that increase in coverage does not necessarily indicate 
        // a multiplicity change
        // however, we are trying to avoid false positives - false negatives will 
        // most likely be cut off from the component or would be removed by remove_tips
        // as a next step
        // increasing these values will increase the number of removed edges
        const COV_STATE_FACTOR: f32 = 0.2; 
        const COV_STATE_ADD: f32 = 2.;

        // path, including start and target node
        let mut paths = Vec::new();
        paths.push(Vec::new());
        paths[0].push(start_node_id);

        let target_length = K::k();
        let mut path_length = 0;
        let mut n_correct_edges = 0;

        // get fist next node
        let sequence = self.base.sequences.get(start_node_id);
        let term_kmer: K = sequence.term_kmer(Dir::Right);
        let next_kmer = term_kmer.extend(start_ext, Dir::Right);
        let (mut current_node_id, _, _) = self.find_link(next_kmer, Dir::Right).expect("link should exist"); 
        paths[0].push(current_node_id);

        if self.check_edge_truth(start_node_id, current_node_id) {
            n_correct_edges += 1;
        }

        let mut current_cov = start_cov as f32;
        let mut sum_path_cov = start_cov as f32;
        let mut coverage_counter = 1; 

        // state to track if nodes should be removed or we're suspecting there's another path here
        // only track edges for removal while state is Singular
        // only track one additional path, else too complicated, return None
        // increase state if there is another incoming edge or if the coverage suddenly increases drastically
        // decrease state if there is another outgoing edge ot if the coverage suddenly decreases drastically
        // if state cannot be increased or decreased further, return None
        let mut state = LadderState::Singular;

        loop {
            // check if current node is "target node", i.e., the path has reached the desired length
            match path_length {
                len if len == target_length => return Some((paths, (sum_path_cov / coverage_counter as f32), (n_correct_edges as f32 / (path_length + 1) as f32))), // path has target length
                len if len > target_length => return None, // path has surpassed desired length => invalid
                _ => () // path has not yet reached desired length, continue
            }

            // since current node is not target node, add length to path
            let len = self.get_node(current_node_id).len() - (K::k() - 1);
            path_length += len;

            // get current node
            let current_node = self.get_node(current_node_id);
            
            // get edges
            let in_edges = current_node.l_edges();
            let out_edges = current_node.r_edges();

            // get coverage
            // edge cov must be available
            let edge_coverages = current_node.data().edge_mults().expect("must have edge mults");

            // choose next edge by choosing coverage closest to current coverage
            let mut out_ext = None;
            let mut out_node_id = None;
            let mut cov_diff = i32::MAX;
            for (e, id, _, _) in out_edges.iter() {
                let cov = edge_coverages.edge_mult(*e, Dir::Right) as i32;
                let new_cov_diff = (current_cov as i32 - cov).abs();
                if new_cov_diff < cov_diff {
                    (out_ext, out_node_id, cov_diff) = (Some(*e), Some(*id), new_cov_diff);
                }
            }
            let (Some(out_ext), Some(out_node_id)) = (out_ext, out_node_id) else { return None; };

            // ideally, low path nodes should have one incoming and one outgoing edge
            // if two incoming edges, increase state from singular to double if possible
            match in_edges.len() {
                1 => (),
                2 => match state {
                    LadderState::Singular => state = LadderState::Double,
                    LadderState::Double => return None // getting too complicated
                }
                _ => return None
            }

            // if two outgoing edges, decrease state from double to singular if possible
            // if already singular, 
            match out_edges.len() {
                1 => (),
                2 => match state {
                    LadderState::Singular => (), // could return None, but would kill if overlap is only one node long
                    LadderState::Double => {
                        state = LadderState::Singular;
                        paths.push(vec![current_node_id]); // start new path
                    }
                }
                _ => return None
            }

            // check coverage, in theoretical ladder, coverage should be uniform
            let coverage = edge_coverages.edge_mult(out_ext, Dir::Right) as f32;
            match state {
                LadderState::Singular => {
                    // single state: check if next coverage is similar enough to current coverage
                    // if way bigger, increase state, else adapt current coverage
                    if  coverage > current_cov + current_cov * COV_STATE_FACTOR + COV_STATE_ADD {
                        // higher by too much, change state
                        state = LadderState::Double;
                    } else {
                        // in acceptable frame
                        current_cov = coverage;
                    }
                }
                LadderState::Double => {
                    // double state: if way smaller (back in acceptable frame), decrease state, else dont treat as current coverage
                    if coverage < current_cov + current_cov * COV_STATE_FACTOR + COV_STATE_ADD {
                        // back in acceptable range
                        state = LadderState::Singular;
                        paths.push(vec![current_node_id]); // start new path
                        current_cov = coverage;
                    } // else continue on 
                    // TODO check if better to also interrupt if coverage increases further
                }
            }

            // if state singular try to check if edge is "correct", add extra length of current node to it to account for compression
            match state {
                LadderState::Singular => {
                    if self.check_edge_truth(current_node_id, out_node_id) {
                        n_correct_edges += len;
                    }
                }
                LadderState::Double => ()
            }

            // set current node id to next node to be visited
            current_node_id = out_node_id;

            // add current node to path, unless state is double
            match state {
                LadderState::Double => (),
                LadderState::Singular => {
                    sum_path_cov += current_cov;
                    coverage_counter += 1;
                    let last_path = paths.len() - 1;
                    paths[last_path].push(current_node_id); // add node to latest path
                }
            }
        }
    }

    /// follow a path of the length 2*k - 1 by choosing the edges with the hightest coverage
    /// requres the graph to have edge mults
    fn follow_ladder_path_high<DI>(&self, start_node_id: usize, start_ext: u8, start_cov: u32) -> Option<(usize, f32, f32)> 
    where SD: SummaryData<DI>
    {

        let target_length = K::k();
        let mut path_length = 0;
        let mut sum_path_cov = start_cov; 
        let mut coverage_counter = 1;
        let mut n_correct_edges = 0; 

        // get fist next node
        let sequence = self.base.sequences.get(start_node_id);
        let term_kmer: K = sequence.term_kmer(Dir::Right);
        let next_kmer = term_kmer.extend(start_ext, Dir::Right);
        let (mut current_node_id, _, _) = self.find_link(next_kmer, Dir::Right).expect("link should exist"); 

        if self.check_edge_truth(start_node_id, current_node_id) {
            n_correct_edges += 1;
        }

        loop {
            // check if current node is "target node", i.e., the path has reached the desired length
            match path_length {
                pl if pl == target_length => return Some((current_node_id, (sum_path_cov as f32 / coverage_counter as f32), (n_correct_edges as f32 / (path_length + 1) as f32))), // path has target length
                pl if pl > target_length => return None, // path has surpassed desired length => invalid
                _ => () // path has not yet reached desired length, continue
            }

            // since current node is not target node, add length to path
            let len = self.get_node(current_node_id).len() - (K::k() - 1);
            path_length += len;

            // get current node
            let current_node = self.get_node(current_node_id);
        
            // find outgoing edge with highest coverage
            let em = current_node.data().edge_mults().expect("must have edge mults").right();
            let (max_cov_base, &max_cov) = em.iter().rev().enumerate().filter(|&(_, &c)| c > 0).max_by(|&(_b1, &c1), &(_b2, c2)| c1.cmp(c2))?;

            // get next node id
            let sequence = self.base.sequences.get(current_node_id);
            let term_kmer: K = sequence.term_kmer(Dir::Right);
            let next_kmer = term_kmer.extend(max_cov_base as u8, Dir::Right);
            let (out_node_id, _, _) = self.find_link(next_kmer, Dir::Right).expect("link should exist"); 

            // try to check if edge is "correct", add extra length of current node to it to account for compression
            if self.check_edge_truth(current_node_id, out_node_id) {
                n_correct_edges += len;
            }

            // set current node id to next node to be visited
            current_node_id = out_node_id;
            sum_path_cov += max_cov;
            coverage_counter += 1;
        }
    }

    /// remove tips from the graph, reqires the graoh to have edge mults and be stranded
    /// it is recommended to use this function after [`DebruijnGraph::remove_ladders`], since 
    /// the latter will likely leave tips in the graph
    pub fn remove_tips<DI, P>(&mut self, min_diff_factor: u32, max_avg_tip_cov: f32, out_path: Option<P>) -> Result<(), String> 
    where 
        SD: SummaryData<DI>,
        P: AsRef<Path>
    {
        // interrupt if we dont have edge mults
        if self.get_node(0).data().edge_mults().is_none() { return Err(String::from("no edge mults available")) };
        if !self.base.stranded { return Err(String::from("must be stranded")) };

        let mut writer = out_path.map(|path| BufWriter::new(File::create(path).expect("error creating ladder stats file")));
        if let Some(wtr) = writer.as_mut() {
            writeln!(wtr, "high cov,avg low cov,max true,tip truth ratio,tip len,dir").unwrap();
        }

        let max_len = K::k();

        // iterate over nodes
        for node_id in 0..self.len() {
            // do for both left and right as "outgoing" direction 
            // right for regular outgoing and tips from the read-end
            // left for the reverse -> tips from the read-start
            for dir in [Dir::Left, Dir::Right] {
                let current_node = self.get_node(node_id);
                // check if node has outgoing edge(s) with both high and low coverage
                let outs = current_node.data().edge_mults().expect("should have em").single_dir(dir).edge_mults;

                let Some((out_max_base, &out_max_cov)) = outs.iter().rev().enumerate().filter(|&(_, &c)| c  > 0).max_by(|&(_b1, &c1), &(_b2, c2)| c1.cmp(c2)) else { continue };
                let smaller_outs = outs.iter().copied().rev().enumerate().filter(|&(_b, c)| (c > 0) & (out_max_cov > c)).collect::<Vec<_>>();
                
                if smaller_outs.is_empty() { continue; }

                // try and check if max cov edge is correct
                let max_connection_correct = {
                    let sequence = self.base.sequences.get(node_id);
                    let term_kmer: K = sequence.term_kmer(dir);
                    let next_kmer = term_kmer.extend(out_max_base as u8, dir);
                    let (next_node_id, _, _) = self.find_link(next_kmer, dir).expect("missing link");

                    self.check_edge_truth(node_id, next_node_id)
                }; // if current node is incorrect, connection is incorrect in any case
                
                // check all small outs
                for (s_base, s_cov) in smaller_outs {
                    // path has to be short + unambiguous + consistently low coverage
                    let Some((tip_path, avg_tip_coverage, truth_ratio, tip_len)) = self.follow_tip_path(node_id, s_base as u8, s_cov, dir) else { continue; };

                    // remove path if coverage below threshold
                    if (s_cov * min_diff_factor <= out_max_cov) & (avg_tip_coverage <= max_avg_tip_cov) & (tip_len <= max_len) {
                        self.remove_path(tip_path, dir)?;
                    }
                        
                    if let Some(wtr) = writer.as_mut() {
                        writeln!(wtr, "{},{},{},{},{},{:?}", out_max_cov, avg_tip_coverage, max_connection_correct, truth_ratio, tip_len, dir).unwrap();
                    }
                }
            }
        }
        
        Ok(())
    }

    /// follow a tip path, returns path, average coverage and ration of correct to incorrect connections
    fn follow_tip_path<DI>(&self, start_node_id: usize, start_ext: u8, start_cov: u32, dir: Dir) -> Option<(Vec<usize>, f32, f32, usize)> 
    where 
        SD: SummaryData<DI>
    {
        const COV_MARGIN: f32 = 0.3;
        const COV_ADD_MARGIN: f32 = 4.;

        // path, including start and target node
        let mut path = Vec::new();
        path.push(start_node_id);

        let mut path_length = 0;
        let mut n_correct_edges = 0;

        // get fist next node
        let sequence = self.base.sequences.get(start_node_id);
        let term_kmer: K = sequence.term_kmer(dir);
        let next_kmer = term_kmer.extend(start_ext, dir);
        let (mut current_node_id, _, _) = self.find_link(next_kmer, dir).expect("link should exist"); 
        path.push(current_node_id);

        if self.check_edge_truth(start_node_id, current_node_id) {
            n_correct_edges += 1;
        }

        let mut current_cov = start_cov as f32;
        let mut sum_path_cov = start_cov as f32;
        let mut coverage_counter = 1;

        loop {     
            // get current node
            let current_node = self.get_node(current_node_id);

            // add length to path
            let len = current_node.len() - (K::k() - 1);
            path_length += len;

            // low path nodes should have only one incoming and one outgoing edge
            let in_edges = current_node.edges(dir.flip());
            if in_edges.len() != 1 { return None; }

            let out_edges = current_node.edges(dir);
            match out_edges.len() {
                oe if oe > 1 => return None,
                0 => return Some((path, (sum_path_cov / coverage_counter as f32), (n_correct_edges as f32 / path_length as f32), path_length)),
                _ => ()
            }

            // if edge cov avalable, check for similarity
            let (out_ext, out_node_id, _, _) = out_edges[0];
            if let Some(em) = current_node.data().edge_mults() {
                let coverage = em.edge_mult(out_ext, dir) as f32;
                if coverage < current_cov - current_cov * COV_MARGIN - COV_ADD_MARGIN // extra two for small values
                    || coverage > current_cov + current_cov * COV_MARGIN + COV_ADD_MARGIN {
                    return None;
                } else {
                    current_cov = coverage;
                }
            }

            // try to check if edge is "correct", add extra length of current node to it to account for compression
            if self.check_edge_truth(current_node_id, out_node_id) {
                n_correct_edges += len;
            }

            // set current node id to next node to be visited
            current_node_id = out_node_id;
            sum_path_cov += current_cov;
            coverage_counter += 1;
            // add current node to path
            path.push(current_node_id);

        }
    }

    /// remove the edges (not the nodes) of a path in the graph
    fn remove_path<DI>(&mut self, path: Vec<usize>, base_dir: Dir) -> Result<(), String> 
    where SD: SummaryData<DI>
    {
        let mut path_iter = path.clone().into_iter();
        let Some(mut current_node_id) = path_iter.next() else { return Ok(()); };

        loop {
            // get next node id
            let Some(next_node_id) = path_iter.next() else { return Ok(()); }; // whole path has been covered
            // remove ext to the right of current node
            let Some((out_base, _, _, _)) = self.get_node(current_node_id).edges(base_dir).iter().find(|&(_, id, _, _)| *id == next_node_id).copied()
                else { return Err(format!(
"no edge to remove (base dir)
path: {:?}
dir: {:?}
node1:
    id: {current_node_id}
    left e: {:?}
    right e: {:?}
    exts: {:?}
    data: {:?}
node2:
    id: {next_node_id}
    left e: {:?}
    right e: {:?}
    exts: {:?}
    data: {:?}",
                path,
                base_dir,
                self.get_node(current_node_id).l_edges(),
                self.get_node(current_node_id).r_edges(),
                self.get_node(current_node_id).exts(),
                self.get_node(current_node_id).data(),
                self.get_node(next_node_id).l_edges(),
                self.get_node(next_node_id).r_edges(),
                self.get_node(next_node_id).exts(),
                self.get_node(next_node_id).data(),
                ))
            };

            // remove ext 
            self.base.exts[current_node_id] = self.base.exts[current_node_id].remove(base_dir, out_base);
        

            // remove ext to the left of the next node
            let Some((in_base, _, _, _)) = self.get_node(next_node_id).edges(base_dir.flip()).iter().find(|&(_, id, _, _)| *id == current_node_id).copied() 
                else { return Err( format!(
"no edge to remove (other dir)
path: {:?}
dir: {:?}
node1:
    id: {current_node_id}
    left e: {:?}
    right e: {:?}
    exts: {:?}
    data: {:?}
node2:
    id: {next_node_id}
    left e: {:?}
    right e: {:?}
    exts: {:?}
    data: {:?}",
                path,
                base_dir,
                self.get_node(current_node_id).l_edges(),
                self.get_node(current_node_id).r_edges(),
                self.get_node(current_node_id).exts(),
                self.get_node(current_node_id).data(),
                self.get_node(next_node_id).l_edges(),
                self.get_node(next_node_id).r_edges(),
                self.get_node(next_node_id).exts(),
                self.get_node(next_node_id).data(),
                ))
            };

            // remove ext
            self.base.exts[next_node_id] = self.base.exts[next_node_id].remove(base_dir.flip(), in_base); 

            // use new exts to fix edge mults
            self.base.data[current_node_id].fix_edge_mults(self.base.exts[current_node_id]);
            self.base.data[next_node_id].fix_edge_mults(self.base.exts[next_node_id]);


            current_node_id = next_node_id;
        }
    }

    /// use mapped ids to check if the nodes of an edge were mapped to the same id
    /// returns false if mapped ids are not available
    pub fn check_edge_truth<DI>(&self, node_id_1: usize, node_id_2: usize) -> bool
    where 
        SD: SummaryData<DI>
    {
        let mapped_ids_1 = self.get_node(node_id_1).data().mapped_ids();
        let mapped_ids_2 = self.get_node(node_id_2).data().mapped_ids();

        if let (Some(mids1), Some(mids2)) = (mapped_ids_1, mapped_ids_2) {
            let hashed_mi2 = mids2.iter().copied().collect::<HashSet<_>>();
            for mid in mids1 {
                if hashed_mi2.contains(mid) {
                    return true;
                }
            }
        }

        false
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

#[derive(Debug, PartialEq, Eq)]
pub enum LadderState {
    Singular,
    Double
}

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
            format!("[color={color}, penwidth={}]", colors.edge_width(1)) // since there should be no edge mults, this will return default value
        }
    }

    /// get default format for dot nodes, based on node data
    pub fn node_dot_default<DI>(&self, colors: &Colors<SD, DI>, config: &SummaryConfig, translator: &Translator, outline: bool, translate_id_groups: bool) -> String
    where SD: SummaryData<DI>
    {
        // set color based on labels/fold change/p-value
        let color = colors.node_color(self.data(), config, outline);
        let translate_id_groups = if translate_id_groups { colors.id_group_ids() } else { None };

        let data_info = self.data().print(translator, config, translate_id_groups);
        const MIN_TEXT_WIDTH: usize = 40;
        let wrap = if self.len() > MIN_TEXT_WIDTH { self.len() } else { MIN_TEXT_WIDTH };

        let label = textwrap::fill(&format!("id: {}, len: {}, exts: {:?}, seq: {}\n{}", 
            self.node_id,
            self.len(),
            self.exts(),
            self.sequence(),
            data_info
        ), wrap);

        format!("[{color}, label=\"{label}\"]")
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
                        }
                        
                        // break if multiple solid paths are available
                        /* if solid_paths > 1 {
                            break;
                        } */
    
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


/// iterator over the edges of the de bruijn graph
pub struct EdgeIter<'a, K: Kmer, D: Debug> {
    graph: &'a DebruijnGraph<K, D>,
    visited_edges: HashSet<(usize, usize)>,
    current_node: usize,
    current_dir: Dir,
    node_edge_iter: smallvec::IntoIter<[(u8, usize, Dir, bool); 4]>
}

impl<K: Kmer, D: Debug> EdgeIter<'_, K, D> {
    pub fn new(graph: &DebruijnGraph<K, D>) -> EdgeIter<'_, K, D>{
        let node_edge_iter = graph.get_node(0).l_edges().into_iter();

        EdgeIter { 
            graph, 
            visited_edges: HashSet::new(), 
            current_node: 0, 
            current_dir: Dir::Left, 
            node_edge_iter
        }
    }
}

impl<K: Kmer, D: Debug> Iterator for EdgeIter<'_, K, D> {
    type Item = (usize, Dir, u8, usize); // node, direction leaving node, base, target node

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some((base, nb_node_id, _, _)) = self.node_edge_iter.next() {
                let edge = if self.current_node > nb_node_id { (nb_node_id, self.current_node) } else { (self.current_node, nb_node_id) };

                if self.visited_edges.insert(edge) { return Some((self.current_node, self.current_dir, base, nb_node_id)); } // else simply skip and move on

            } else {
                match self.current_dir {
                    Dir::Left => {
                        // no left edges, switch to right edges
                        self.current_dir = Dir::Right;
                        self.node_edge_iter = self.graph.get_node(self.current_node).r_edges().into_iter();
                        
                    }
                    Dir::Right => {
                        // no right edges, switch to next node left edges
                        self.current_node += 1;

                        // quit if end of graph is reached
                        if self.current_node == self.graph.len() { return None }

                        self.current_dir = Dir::Left;
                        self.node_edge_iter = self.graph.get_node(self.current_node).l_edges().into_iter();
                    }
                }
            }
            
        }
    }
}

#[cfg(test)]
mod test {
    use std::fs::remove_file;

    use crate::{BaseQuality, Exts, colors::Colors, compression::{CheckCompress, ScmapCompress, compress_kmers_with_hash, uncompressed_graph}, dna_string::DnaString, filter::filter_kmers, kmer::{Kmer6, Kmer16, Kmer22}, reads::{Reads, ReadsPaired}, serde::SerKmers, summarizer::{IDMapEMData, IDMapEMQualityData, IDTag, SampleInfo, SummaryConfig, TagsCountsData, TagsCountsSumData, Translator}, test::random_dna};

    use crate::{summarizer::SummaryData, Dir};


    #[test]
    #[cfg(not(feature = "sample128"))]
    fn test_components() {

        // cargo run -- -i data/test_400.fastq.gz -o ../rust-debruijn/test_data/400 -s tags-counts-sum -k 18 -t t

        use crate::{kmer::Kmer16, serde::SerGraph};
        let path = "test_data/400.graph.dbg";

        let ser_graph: SerGraph<Kmer16, TagsCountsSumData> = SerGraph::deserialize_from(path);
        let graph = ser_graph.graph();

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
        assert_eq!(vec![(139, Dir::Left)], graph.max_path(|data| data.sum().unwrap_or(1) as f32, |_| true));
    }

    #[test]
    fn test_iter_edges() {
        use crate::{compression::uncompressed_graph, filter::filter_kmers, reads::{Reads, ReadsPaired}, summarizer::{SampleInfo, SummaryConfig, TagsData}};

        let read1 = "CAGCATCGATGCGACGAGCGCTCGCATCGA".as_bytes();
        let read2 = "ACGATCGTACGTAGCTAGCTGACTGAGC".as_bytes();

        let mut reads = Reads::new(crate::reads::Strandedness::Forward);
        reads.add_from_bytes(read1, None, 0u8);
        reads.add_from_bytes(read2, None, 1);

        let reads_paired = ReadsPaired::Unpaired { reads };

        let sample_info = SampleInfo::new(0b1, 0b10, vec![12, 12]);
        let summary_config = SummaryConfig::new(sample_info);
        let (kmers, _) = filter_kmers::<TagsData, Kmer16, _>(&reads_paired, &summary_config, false, 1., false);

        let graph = uncompressed_graph(&kmers, true).finish();

        let check_edges: Vec<(usize, Dir, u8, usize)> = vec![(0, Dir::Left, 2, 16), (0, Dir::Right, 1, 2), (1, Dir::Left, 0, 11), 
        (1, Dir::Right, 0, 16), (3, Dir::Left, 0, 13), (3, Dir::Right, 3, 10), (4, Dir::Left, 2, 21), (4, Dir::Right, 1, 17), (5, Dir::Left, 1, 27), 
        (5, Dir::Right, 2, 19), (6, Dir::Left, 1, 24), (6, Dir::Right, 2, 12), (7, Dir::Left, 2, 23), (7, Dir::Right, 3, 11), (8, Dir::Left, 0, 12), 
        (9, Dir::Left, 3, 17), (9, Dir::Right, 3, 24), (10, Dir::Right, 1, 21), (13, Dir::Left, 1, 18), (14, Dir::Right, 0, 26), 
        (15, Dir::Left, 2, 20), (15, Dir::Right, 3, 22), (18, Dir::Left, 2, 19), (20, Dir::Left, 1, 26), (22, Dir::Right, 2, 25), (23, Dir::Left, 1, 25)];

        let edges = graph.iter_edges().collect::<Vec<_>>();

        assert_eq!(check_edges, edges);
    }

    // dbg -c ../marbel_datasets/sim_reads_100.csv -s sum --stranded -o ../rust-debruijn/test_data/marbel_100_sum --checkpoint -k 22
    #[cfg(not(feature = "sample128"))]
    const TEST_GRAPH: &str = "test_data/marbel_100_sum.kmers.dbg";

    // cargo run --features sample128 -- -c ../marbel_datasets/sim_reads_100.csv -s sum --stranded -o ../rust-debruijn/test_data/marbel_100_sum_128 --checkpoint -k 22
    #[cfg(feature = "sample128")]
    const TEST_GRAPH: &str = "test_data/marbel_100_sum_128.kmers.dbg";

    #[test]
    fn test_map_transcripts() {
        // dbg -c ../marbel_datasets/sim_reads_100.csv -s sum --stranded -o ../rust-debruijn/test_data/marbel_100_sum --checkpoint -k 22
        let graph_path = TEST_GRAPH;
        let t_ref_path = "test_data/marbel_100_tr_ref.fasta";
        let (kmers, mut translator, _) = SerKmers::<Kmer22, u32>::deserialize_from(graph_path).dissolve();

        let unc_graph = uncompressed_graph(&kmers, true).finish();

        let t_map = unc_graph.map_transcripts(t_ref_path, &mut translator).unwrap();
        assert_eq!(t_map.len(), unc_graph.len());
        assert_eq!(t_map.iter().filter(|&ids| !ids.is_empty()).collect::<Vec<_>>().len(), 10720);
    }

    fn build_reads_quality_test() -> ReadsPaired<IDTag> {
        let correct1 = "CGATGCTGCTGATGCTGAGTCTGACGTATGCGATCGATCGACGATCGTACTAGCTGACTGTGCAGCTAGCTGACTGATCGTAGCTAGCTACGTGCTAGCTACTAGCACTGATGC";
        let qu_corr1 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
        let incorrect1 =                                     "CGACGATCGTACTAGCTGACTGTGCAGCTAGCTGACTGATCGTGGCTAGCTACGTGCTAGCTA";
        let qu_incorr1 =                                     "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCC";
        let correct2 =         "GCATCGATCGACGACGTTACGTACGATCTACGTAGCTAGCTAGCTGACATGCTAGCTAGCTCTGACTGATCGTGGCTAGCTGACTGACTGTAGCT";
        let qu_corr2 =         "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
        let incorrect2 = "ATGCTGCTGATGCTGAGTCTGACGTAGGCGATCGATCGACGATCGTACTAGCTGACT";
        let qu_incorr2 = "CCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
        let incorrect3 = "ATGCTGCTGATGCTGACTCTGACGTAGGCGATCGATCGATGATCGTACTAGCTGACT";
        let qu_incorr3 = "CCCCCCCCCCCCCCCC-CCCCCCCCC-CCCCCCCCCCCC-CCCCCCCCCCCCCCCCC";

        let incorrect4_ins =                 "GACGTATGCGATCGATCGACGATCGTACTAGCTGACTTGTGCAGCTAGCTGACTGAT";
        let qu_incorr4_ins =                 "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCC";

        let incorrect5_tip =                                                             "AGCTGACTGATCGTAGCTAGCTACGTGCTAGCTACTATCACTGATGC";
        let qu_incorr5_tip =                                                             "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCC";


        let mut reads = Reads::new_with_quality(crate::reads::Strandedness::Forward);

        for _i in 0..20 {
            reads.add_read(DnaString::from_acgt_bytes(correct1.as_bytes()), None, IDTag::new(1, 0), Some(qu_corr1.as_bytes()));
        }

        for _i in 0..40 {
            reads.add_read(DnaString::from_acgt_bytes(correct2.as_bytes()), None, IDTag::new(2, 0), Some(qu_corr2.as_bytes()));
        }

        reads.add_read(DnaString::from_acgt_bytes(incorrect1.as_bytes()), None, IDTag::new(3, 1), Some(qu_incorr1.as_bytes()));
        reads.add_read(DnaString::from_acgt_bytes(incorrect2.as_bytes()), None, IDTag::new(4, 1), Some(qu_incorr2.as_bytes()));
        reads.add_read(DnaString::from_acgt_bytes(incorrect3.as_bytes()), None, IDTag::new(5, 1), Some(qu_incorr3.as_bytes()));
        reads.add_read(DnaString::from_acgt_bytes(incorrect4_ins.as_bytes()), None, IDTag::new(6, 1), Some(qu_incorr4_ins.as_bytes()));
        reads.add_read(DnaString::from_acgt_bytes(incorrect5_tip.as_bytes()), None, IDTag::new(7, 1), Some(qu_incorr5_tip.as_bytes()));


        ReadsPaired::Unpaired { reads }
    }

    #[test]
    fn test_remove_lq_splits() {
        let print = false;

        type K = Kmer16;

        let seqs = build_reads_quality_test();
        let sample_info = SampleInfo::new(1, 0b111110, vec![1000, 10, 20, 20, 20, 20]);
        let summary_config = SummaryConfig::new(sample_info);
        let (kmers, _) = filter_kmers::<IDMapEMQualityData, K, IDTag>(&seqs, &summary_config, false, 1., false);


        // make uncompressed graph
        let mut unc_graph = uncompressed_graph(&kmers, true).finish();

        // add "mapped" ids to graph
        for i in 0..unc_graph.len() {
            let data = unc_graph.mut_data(i);
            if let Some(ids) = data.ids() {
                if ids.contains(&1) {
                    data.set_mapped_ids(vec![1].into());
                }
                else if ids.contains(&2) {
                    data.set_mapped_ids(vec![2].into());
                }
            }
        }

        let colors = Colors::new(&unc_graph, &summary_config, crate::colors::ColorMode::IDS { n_ids: 7 });

        if print { unc_graph.to_dot("uncompressed_bf.dot", &|node| node.node_dot_default(&colors, &summary_config, &Translator::empty(), false, false), &|node, base, dir, flip| node.edge_dot_default(&colors, base, dir, flip)); }
        let n_edges = unc_graph.iter_edges().count();
        unc_graph.remove_lq_splits(BaseQuality::Medium).unwrap();
        if print { unc_graph.to_dot("uncompressed_af.dot", &|node| node.node_dot_default(&colors, &summary_config, &Translator::empty(), false, false), &|node, base, dir, flip| node.edge_dot_default(&colors, base, dir, flip)); }
    
        let n_edges_af = unc_graph.iter_edges().count();
        assert_eq!(n_edges, n_edges_af + 11);

        // make compressed graph
        let spec = CheckCompress::new(|d: IDMapEMQualityData, _| d, |d, d1| d.join_test(d1));
        let mut c_graph = compress_kmers_with_hash(true, &spec, &kmers, false, false).finish();
    
         // add "mapped" ids to graph
        for i in 0..c_graph.len() {
            let data = c_graph.mut_data(i);
            if let Some(ids) = data.ids() {
                if ids.contains(&1) {
                    data.set_mapped_ids(vec![1].into());
                }
                else if ids.contains(&2) {
                    data.set_mapped_ids(vec![2].into());
                }
            }
        }

        let colors = Colors::new(&c_graph, &summary_config, crate::colors::ColorMode::IDS { n_ids: 7 });

        if print { c_graph.to_dot("compressed_bf.dot", &|node| node.node_dot_default(&colors, &summary_config, &Translator::empty(), false, false), &|node, base, dir, flip| node.edge_dot_default(&colors, base, dir, flip)); }
        let n_edges = c_graph.iter_edges().count();
        c_graph.remove_lq_splits(BaseQuality::Medium).unwrap();
        if print { c_graph.to_dot("compressed_af.dot", &|node| node.node_dot_default(&colors, &summary_config, &Translator::empty(), false, false), &|node, base, dir, flip| node.edge_dot_default(&colors, base, dir, flip)); }
    
        let n_edges_af = c_graph.iter_edges().count();
        assert_eq!(n_edges, n_edges_af + 11);
    }

    #[test]
    fn test_remove_lq_ladders_tips() {
        let print = false;

        type K = Kmer16;

        let seqs = build_reads_quality_test();
        let sample_info = SampleInfo::new(1, 0b111110, vec![1000, 10, 20, 20, 20, 20]);
        let summary_config = SummaryConfig::new(sample_info);
        let (kmers, _) = filter_kmers::<IDMapEMQualityData, K, IDTag>(&seqs, &summary_config, false, 1., false);


        // make uncompressed graph
        let mut unc_graph = uncompressed_graph(&kmers, true).finish();

        // add "mapped" ids to graph
        for i in 0..unc_graph.len() {
            let data = unc_graph.mut_data(i);
            if let Some(ids) = data.ids() {
                if ids.contains(&1) {
                    data.set_mapped_ids(vec![1].into());
                }
                else if ids.contains(&2) {
                    data.set_mapped_ids(vec![2].into());
                }
            }
        }

        let colors = Colors::new(&unc_graph, &summary_config, crate::colors::ColorMode::IDS { n_ids: 7 });

        if print { unc_graph.to_dot("uncompressed_bf.dot", &|node| node.node_dot_default(&colors, &summary_config, &Translator::empty(), false, false), &|node, base, dir, flip| node.edge_dot_default(&colors, base, dir, flip)); }
        let n_edges = unc_graph.iter_edges().count();
        unc_graph.remove_lq_ladders_tips(BaseQuality::Medium).unwrap();
        if print { unc_graph.to_dot("uncompressed_af.dot", &|node| node.node_dot_default(&colors, &summary_config, &Translator::empty(), false, false), &|node, base, dir, flip| node.edge_dot_default(&colors, base, dir, flip)); }
    
        let n_edges_af = unc_graph.iter_edges().count();
        assert_eq!(n_edges, n_edges_af + 90);

        // make compressed graph
        let spec = CheckCompress::new(|d: IDMapEMQualityData, _| d, |d, d1| d.join_test(d1));
        let mut c_graph = compress_kmers_with_hash(true, &spec, &kmers, false, false).finish();
    
         // add "mapped" ids to graph
        for i in 0..c_graph.len() {
            let data = c_graph.mut_data(i);
            if let Some(ids) = data.ids() {
                if ids.contains(&1) {
                    data.set_mapped_ids(vec![1].into());
                }
                else if ids.contains(&2) {
                    data.set_mapped_ids(vec![2].into());
                }
            }
        }

        let colors = Colors::new(&c_graph, &summary_config, crate::colors::ColorMode::IDS { n_ids: 7 });

        if print { c_graph.to_dot("compressed_bf.dot", &|node| node.node_dot_default(&colors, &summary_config, &Translator::empty(), false, false), &|node, base, dir, flip| node.edge_dot_default(&colors, base, dir, flip)); }
        let n_edges = c_graph.iter_edges().count();
        c_graph.remove_lq_ladders_tips(BaseQuality::Medium).unwrap();
        if print { c_graph.to_dot("compressed_af.dot", &|node| node.node_dot_default(&colors, &summary_config, &Translator::empty(), false, false), &|node, base, dir, flip| node.edge_dot_default(&colors, base, dir, flip)); }
    
        let n_edges_af = c_graph.iter_edges().count();
        assert_eq!(n_edges, n_edges_af + 15);
    }

    #[test]
    fn test_remove_ladders() {
        let print = false; 
        let c_csv = if print { Some("c_ladders.csv") } else { None };
        let uc_csv = if print { Some("uc_ladders.csv") } else { None };

        let   correct = "ACGATCGATCGCGATCGTAGCTGACTGCTGACGTCTGACTACTGACTGATGCTAGCTATCGTGAC".as_bytes();
        let incorrect = "ACGATCGATCGCGATCGTAGCTGACTGCTGACGGCTGACTACTGACTGATGCTAGCTATCGTGAC".as_bytes();
        let incorrec2 = "TGACAGCTGACGGCTGACTACTACGTCACTGACGATGCTGACAC".as_bytes();
        let incorrec3 = "AAAAAAAAAGCTGACTGCTGACGGCTG".as_bytes();
        let incorrec4 = "ACGGCTGACTACTGACTGAAAAAAAAAAA".as_bytes();

        let insertion = "ACGATCGATCGCGATCGATAGCTGACTGCTGACGTCTGACTACTGACTGATGCTAGCTATCGTGAC".as_bytes();

        let mut reads = Reads::new(crate::reads::Strandedness::Forward);
        for _i in 0..1000 {
            reads.add_from_bytes(correct, None, IDTag::new(0, 0));
        }

        for _i in 0..2 {
            reads.add_from_bytes(incorrect, None, IDTag::new(1, 1)); // should be removed
        }

        for _i in 0..15 {
            reads.add_from_bytes(incorrec2, None, IDTag::new(2, 2)); // should be removed
        }

        for _i in 0..25 {
            reads.add_from_bytes(incorrec3, None, IDTag::new(3, 3)); // should be removed
        }

        for _i in 0..30 {
            reads.add_from_bytes(incorrec4, None, IDTag::new(4, 4)); // should be removed
        }


        for _i in 0..1 {
            reads.add_from_bytes(insertion, None, IDTag::new(1, 3)); // should not be removed
        }

        let seqs = ReadsPaired::Unpaired { reads };
        let sample_info = SampleInfo::new(1, 0b111110, vec![1000, 10, 20, 20, 20, 20]);
        let summary_config = SummaryConfig::new(sample_info);
        let (kmers, _) = filter_kmers::<IDMapEMData, Kmer16, IDTag>(&seqs, &summary_config, false, 1., false);


        // test with uncompressed graph
        let mut unc_graph = uncompressed_graph(&kmers, true).finish();
        // add ids to graph
        for i in 0..unc_graph.len() {
            let data = unc_graph.mut_data(i);
            if let Some(ids) = data.ids() {
                if ids.contains(&0) {
                    data.set_mapped_ids(vec![0].into());
                }
            }
        }
        let colors = Colors::new(&unc_graph, &summary_config, crate::colors::ColorMode::IDS { n_ids: 5 });
        if print { unc_graph.to_dot("uncompressed_bf.dot", &|node| node.node_dot_default(&colors, &summary_config, &Translator::empty(), false, false), &|node, base, dir, flip| node.edge_dot_default(&colors, base, dir, flip)); }
        let n_edges = unc_graph.iter_edges().count();
        unc_graph.remove_ladders(10, 10., uc_csv).unwrap();
        if print { unc_graph.to_dot("uncompressed_af.dot", &|node| node.node_dot_default(&colors, &summary_config, &Translator::empty(), false, false), &|node, base, dir, flip| node.edge_dot_default(&colors, base, dir, flip)); }
        assert_eq!(n_edges - 10, unc_graph.iter_edges().count());

        // test with compressed graph
        let spec = CheckCompress::new(|d: IDMapEMData, _| d, |d, d1| d.join_test(d1));
        let mut c_graph = compress_kmers_with_hash(true, &spec, &kmers, false, false).finish();
        // add ids to graph
        for i in 0..c_graph.len() {
            let data = c_graph.mut_data(i);
            if let Some(ids) = data.ids() {
                if ids.contains(&0) {
                    data.set_mapped_ids(vec![0].into());
                }
            }
        }
        let colors = Colors::new(&c_graph, &summary_config, crate::colors::ColorMode::IDS { n_ids: 5 });
        if print { c_graph.to_dot("compressed_bf.dot", &|node| node.node_dot_default(&colors, &summary_config, &Translator::empty(), false, false), &|node, base, dir, flip| node.edge_dot_default(&colors, base, dir, flip)); }
        let n_edges = c_graph.iter_edges().count();
        c_graph.remove_ladders(10, 10., c_csv).unwrap();
        if print { c_graph.to_dot("compressed_af.dot", &|node| node.node_dot_default(&colors, &summary_config, &Translator::empty(), false, false), &|node, base, dir, flip| node.edge_dot_default(&colors, base, dir, flip)); }
        assert_eq!(n_edges - 6, c_graph.iter_edges().count());
    }

    #[test]
    fn test_remove_tips() {

        let print = false;
        let c_csv = if print { Some("c_tips.csv") } else { None };
        let uc_csv = if print { Some("uc_tips.csv") } else { None };

        let     correct = "ACGATCGATCGCGATCGTAGCTGACTGCTGACGTCTGACTACTGACTGATGCTAGCTATCGTGAC".as_bytes();
        let incorrect_r = "ACGATCGATCGCGATCGTAGCTGACTGCTGACGTCTGACTACTGACTGATGCTAGCTAACGTGAC".as_bytes();
        let incorrect_l = "ACGTTCGATCGCGATCGTAGCTGACTGCTGACGTCTGACTACTGACTGATGCTAGCTATCGTGAC".as_bytes();



        let mut reads = Reads::new(crate::reads::Strandedness::Forward);
        for _i in 0..1000 {
            reads.add_from_bytes(correct, None, IDTag::new(0, 0));
        }

        for _i in 0..10 {
            reads.add_from_bytes(incorrect_r, None, IDTag::new(1, 1)); // should be removed
        }

        for _i in 0..10 {
            reads.add_from_bytes(incorrect_l, None, IDTag::new(2, 2)); // should be removed
        }

        let seqs = ReadsPaired::Unpaired { reads };
        let sample_info = SampleInfo::new(1, 6, vec![1000, 10, 20]);
        let summary_config = SummaryConfig::new(sample_info);
        let (kmers, _) = filter_kmers::<IDMapEMData, Kmer16, IDTag>(&seqs, &summary_config, false, 1., false);


        // test with uncompressed graph
        let mut unc_graph = uncompressed_graph(&kmers, true).finish();
        // add ids to graph
        for i in 0..unc_graph.len() {
            let data = unc_graph.mut_data(i);
            if let Some(ids) = data.ids() {
                if ids.contains(&0) {
                    data.set_mapped_ids(vec![0].into());
                }
            }
        }

        let colors = Colors::new(&unc_graph, &summary_config, crate::colors::ColorMode::IDS { n_ids: 3 });
        if print { unc_graph.to_dot("uncompressed_bf.dot", &|node| node.node_dot_default(&colors, &summary_config, &Translator::empty(), false, false), &|node, base, dir, flip| node.edge_dot_default(&colors, base, dir, flip)); }
        let n_edges = unc_graph.iter_edges().count();
        
        unc_graph.remove_tips(10, 10., uc_csv).unwrap();
        if print { unc_graph.to_dot("uncompressed_af.dot", &|node| node.node_dot_default(&colors, &summary_config, &Translator::empty(), false, false), &|node, base, dir, flip| node.edge_dot_default(&colors, base, dir, flip)); }
        assert_eq!(n_edges - 11, unc_graph.iter_edges().count());

        // test with compressed graph
        let spec = CheckCompress::new(|d: IDMapEMData, _| d, |d, d1| d.join_test(d1));
        let mut c_graph = compress_kmers_with_hash(true, &spec, &kmers, false, false).finish();
        // add ids to graph
        for i in 0..c_graph.len() {
            let data = c_graph.mut_data(i);
            if let Some(ids) = data.ids() {
                if ids.contains(&0) {
                    data.set_mapped_ids(vec![0].into());
                }
            }
        }

        let colors = Colors::new(&c_graph, &summary_config, crate::colors::ColorMode::SampleGroups);
        if print { c_graph.to_dot("compressed_bf.dot", &|node| node.node_dot_default(&colors, &summary_config, &Translator::empty(), false, false), &|node, base, dir, flip| node.edge_dot_default(&colors, base, dir, flip)); }
        let n_edges = c_graph.iter_edges().count();
        
        c_graph.remove_tips(10, 10., c_csv).unwrap();
        if print { c_graph.to_dot("compressed_af.dot", &|node| node.node_dot_default(&colors, &summary_config, &Translator::empty(), false, false), &|node, base, dir, flip| node.edge_dot_default(&colors, base, dir, flip)); }
        assert_eq!(n_edges - 2, c_graph.iter_edges().count());

    }

    #[test]
    fn test_to_tsv() {
        let reads_us = Reads::from_vmer_vec(
            (0..10).map(|i| (DnaString::from_bytes(&random_dna(100)), Exts::empty(), i as u8)).collect::<Vec<_>>(), 
            crate::reads::Strandedness::Unstranded
        );

        let reads_paired = ReadsPaired::Unpaired { reads: reads_us };

        let sample_info = SampleInfo::new(0b1111100000, 0b0000011111, vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
        let summary_config = SummaryConfig::new(sample_info);
        let (kmers, _) = filter_kmers::<TagsCountsData, Kmer6, _>(&reads_paired, &summary_config, false, 1., false);

        let graph = compress_kmers_with_hash(false, &ScmapCompress::new(), &kmers, false, false).finish();

        graph.to_tsv("test_graph_unstranded.tsv", |node| node.data().print_ol(&Translator::empty(), &summary_config, None)).unwrap();


        let reads_us = Reads::from_vmer_vec(
            (0..10).map(|i| (DnaString::from_bytes(&random_dna(100)), Exts::empty(), i as u8)).collect::<Vec<_>>(), 
            crate::reads::Strandedness::Forward
        );

        let reads_paired = ReadsPaired::Unpaired { reads: reads_us };

        let sample_info = SampleInfo::new(0b1111100000, 0b0000011111, vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
        let summary_config = SummaryConfig::new(sample_info);
        let (kmers, _) = filter_kmers::<TagsCountsData, Kmer6, _>(&reads_paired, &summary_config, false, 1., false);

        let graph = compress_kmers_with_hash(true, &ScmapCompress::new(), &kmers, false, false).finish();

        graph.to_tsv("test_graph_stranded.tsv", |node| node.data().print_ol(&Translator::empty(), &summary_config, None)).unwrap();
    
        remove_file("test_graph_unstranded.tsv").unwrap();
        remove_file("test_graph_stranded.tsv").unwrap();
    }
}

