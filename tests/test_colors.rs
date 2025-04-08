use std::{fs::File, io::BufReader};

use bimap::BiMap;
use debruijn::{colors::Colors, graph::DebruijnGraph, kmer::Kmer16, summarizer::TagsCountsPEMData};

#[cfg(not(feature = "sample128"))]
const TEST_FILE: &str = "test_data/sided.graph.dbg";

// use different graph file for sample128 tests
#[cfg(feature = "sample128")]
const TEST_FILE: &str = "test_data/sided-128.graph.dbg";

#[test]
#[cfg(not(feature = "sample128"))]
fn test_colors() {
    let path = TEST_FILE;
    let file = BufReader::new(File::open(path).unwrap());

    let (graph, vec_labels, config): (DebruijnGraph<Kmer16, TagsCountsPEMData>, Vec<String>, debruijn::summarizer::SummaryConfig) = 
        bincode::deserialize_from(file).expect("error deserializing graph");
    let mut hashed_labels: BiMap<String, u8> = BiMap::new();
    for (i, label) in vec_labels.iter().enumerate() {
        hashed_labels.insert(label.clone(), i as u8);
    }

    let colors = Colors::new(&graph, &config);
    println!("colors:{:?}", colors);
    
    let node_id = 0;
    let node = graph.get_node(node_id);
    assert_eq!("[style=filled, color=\"0.33333334 0.5 1\", fontcolor=\"black, label=\"id: 0, len: 16, seq: CCAGCTGTCCCAGATA,\nsamples              - counts\nC_H600               - 1\nsum: 1, p-value: 0.37390098, log2(fold\nchange): inf, edge multiplicities:\nA: 1 | 0\nC: 0 | 1\nG: 0 | 0\nT: 0 | 0\n\n\"]", node.node_dot_default(&colors, &config, &hashed_labels, false));
    assert_eq!("[color=blue, penwidth=1, label=\"1\"]", node.edge_dot_default(&colors, 0, debruijn::Dir::Left, true));

    let node_id = 3;
    let node = graph.get_node(node_id);
    assert_eq!("[style=filled, color=black, penwidth=7, fillcolor=\"0.33333334 0.5 1\", fontcolor=\"black, label=\"id: 3, len: 16, seq: CCGATACTTCTTCACG,\nsamples              - counts\nC_H200               - 1\nC_H300               - 1\nC_H500               - 1\nC_H600               - 1\nsum: 4, p-value: 0.05366346, log2(fold\nchange): inf, edge multiplicities:\nA: 0 | 0\nC: 0 | 4\nG: 0 | 0\nT: 4 | 0\n\n\"]", node.node_dot_default(&colors, &config, &hashed_labels, true));
    assert_eq!("[color=red, penwidth=8.469343, label=\"4\"]", node.edge_dot_default(&colors, 1, debruijn::Dir::Right, true));

    let node_id = 27;
    let node = graph.get_node(node_id);
    assert_eq!("[style=filled, color=black, penwidth=7, fillcolor=\"0 0.5 1\", fontcolor=\"black, label=\"id: 27, len: 16, seq: CAGGTCCGTGGCCAGG,\nsamples              - counts\nP_T400               - 1\nsum: 1, p-value: 0.3910022, log2(fold\nchange): -inf, edge multiplicities:\nA: 0 | 0\nC: 1 | 0\nG: 0 | 1\nT: 0 | 0\n\n\"]", node.node_dot_default(&colors, &config, &hashed_labels, true));
    assert_eq!("[color=red, penwidth=-inf, label=\"0\"]", node.edge_dot_default(&colors, 1, debruijn::Dir::Right, true));

    let node_id = 40;
    let node = graph.get_node(node_id);
    assert_eq!("[style=filled, color=black, penwidth=7, fillcolor=\"0.15937972 0.5 1\", fontcolor=\"black, label=\"id: 40, len: 16, seq: GCTTTGCATCAGAAGA,\nsamples              - counts\nC_H100               - 1\nC_H400               - 2\nP_T200               - 1\nP_T300               - 1\nP_T400               - 1\nsum: 6, p-value: 0.86424315,\nlog2(fold change): -0.19431013, edge\nmultiplicities:\nA: 6 | 3\nC: 0 | 0\nG: 0 | 3\nT: 0 | 0\n\n\"]", node.node_dot_default(&colors, &config, &hashed_labels, true));
    assert_eq!("[color=red, penwidth=-inf, label=\"0\"]", node.edge_dot_default(&colors, 1, debruijn::Dir::Right, true));
    
    let node_id = 2458;
    let node = graph.get_node(node_id);
    assert_eq!("[style=filled, color=black, penwidth=7, fillcolor=\"0 1 1\", fontcolor=\"black, label=\"id: 2458, len: 16, seq:\nTCTCGTACTAAGTTCA, samples              -\ncounts\nC_H600               - 1\nP_T100               - 1\nP_T200               - 1\nP_T300               - 1\nP_T400               - 2\nsum: 6, p-value: 0.030437965,\nlog2(fold change): -4.4442472, edge\nmultiplicities:\nA: 0 | 6\nC: 6 | 0\nG: 0 | 0\nT: 0 | 0\n\n\"]", node.node_dot_default(&colors, &config, &hashed_labels, true));
    assert_eq!("[color=red, penwidth=-inf, label=\"0\"]", node.edge_dot_default(&colors, 1, debruijn::Dir::Right, true));

    // write node to dot

    graph.to_dot(
        "test.dot", 
        &|node| node.node_dot_default(&colors, &config, &hashed_labels, false), 
        &|node, base, dir, flipped| node.edge_dot_default(&colors, base, dir, flipped)
    );

    graph.to_dot_parallel(
        "test.dot", 
        &|node| node.node_dot_default(&colors, &config, &hashed_labels, false), 
        &|node, base, dir, flipped| node.edge_dot_default(&colors, base, dir, flipped)
    );
}