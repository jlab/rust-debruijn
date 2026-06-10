// Copyright 2017 10x Genomics

//! DeBruijn graph simplification routines. Currently tip-removal is implemented.
use crate::graph::{DebruijnGraph, Node};
use crate::Kmer;
use std::fmt::Debug;
use std::marker::PhantomData;

pub struct CleanGraph<K: Kmer, D, T1>
where
    T1: Fn(&Node<'_, K, D>) -> bool,
{
    tip_predicate: T1,
    _k: PhantomData<K>,
    _d: PhantomData<D>,
}

impl<K: Kmer, D: Debug, T1> CleanGraph<K, D, T1>
where
    T1: Fn(&Node<'_, K, D>) -> bool,
{
    pub fn new(tip_predicate: T1) -> CleanGraph<K, D, T1> {
        CleanGraph {
            tip_predicate,
            _k: PhantomData,
            _d: PhantomData,
        }
    }

    fn test_tip(&self, graph: &DebruijnGraph<K, D>, id: usize) -> Option<usize> {
        let node = graph.get_node(id);
        let exts = node.exts();
        if exts.num_exts_r() > 0 && exts.num_exts_l() > 0 {
            return None;
        }

        if ((exts.num_exts_l() == 0 && exts.num_exts_r() <= 1)
            || (exts.num_exts_r() == 0 && exts.num_exts_l() <= 1))
            && (self.tip_predicate)(&node)
        {
            return Some(id);
        }

        None
    }

    pub fn find_bad_nodes(&self, graph: &DebruijnGraph<K, D>) -> Vec<usize> {
        (0..graph.len())
            .filter_map(|i| self.test_tip(graph, i))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::{build_test_graph, clean_graph::CleanGraph, graph::Node, kmer::Kmer16, summarizer::TagsCountsData};

    #[test]
    fn test_clean_graph() {
        let (_ser_reads, _ser_kmers, ser_graph) = build_test_graph::<Kmer16, TagsCountsData, _>();
        let cleaner = CleanGraph::<Kmer16, TagsCountsData, _>::new(|node: &Node<'_, _, TagsCountsData>| (node.l_edges().is_empty() | node.r_edges().is_empty()) & (node.data().sum() == 1));
        let bad_nodes = cleaner.find_bad_nodes(ser_graph.graph());
        assert_eq!(vec![4, 18, 20, 23, 34, 44, 66, 73, 76, 88], bad_nodes);
    }
}
