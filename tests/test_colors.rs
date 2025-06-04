use std::{collections::HashMap, fs::{remove_file, File}, io::BufReader};

use bimap::BiMap;
use debruijn::{colors::{ColorMode, Colors}, graph::{self, DebruijnGraph}, kmer::Kmer16, serde::SerGraph, summarizer::{IDSumData, SummaryConfig, SummaryData, TagsCountsPEMData, Translator, ID}};

#[cfg(not(feature = "sample128"))]
const TEST_FILE_T: &str = "test_data/sided.graph.dbg";

// use different graph file for sample128 tests
#[cfg(feature = "sample128")]
const TEST_FILE_T: &str = "test_data/sided-128.graph.dbg";

const TEST_FILE_IDS: &str = "test_data/test_graph_ids.graph.dbg";

#[test]
fn test_colors() {
    let file_t = BufReader::new(File::open(TEST_FILE_T).unwrap());

    let (graph_tcpem, vec_labels, config_tcpem): (DebruijnGraph<Kmer16, TagsCountsPEMData>, Vec<String>, SummaryConfig) = 
        bincode::deserialize_from(file_t).expect("error deserializing graph");
    let mut hashed_labels_tcpem: BiMap<String, u8> = BiMap::new();
    for (i, label) in vec_labels.iter().enumerate() {
        hashed_labels_tcpem.insert(label.clone(), i as u8);
    }

    let (graph_ids, _, hashed_genes, config_ids) = SerGraph::<Kmer16, IDSumData>::deserialize_from(TEST_FILE_IDS).dissolve();

    let translator = Translator::new(hashed_genes, hashed_labels_tcpem);

    // test with color mode FoldChange

    let colors = Colors::new(&graph_tcpem, &config_tcpem, ColorMode::FoldChange);
    println!("colors:{:?}", colors);
    
    let node_id = 0;
    let node = graph_tcpem.get_node(node_id);
    assert_eq!("[style=filled, color=\"0.33333334 0.5 1\", fontcolor=black, label=\"id: 0, len: 16, exts: A|C, seq:\nCCAGCTGTCCCAGATA\nsamples              - counts\nC_H600               - 1\nsum: 1, p-value: 0.37390098, log2(fold\nchange): inf, edge multiplicities:\nA: 1 | 0\nC: 0 | 1\nG: 0 | 0\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, false));
    assert_eq!("[color=blue, penwidth=1, label=\"A: 1\"]", node.edge_dot_default(&colors, 0, debruijn::Dir::Left, true));

    let node_id = 3;
    let node = graph_tcpem.get_node(node_id);
    assert_eq!("[style=filled, color=black, penwidth=7, fillcolor=\"0.33333334 0.5 1\", fontcolor=black, label=\"id: 3, len: 16, exts: T|C, seq:\nCCGATACTTCTTCACG\nsamples              - counts\nC_H200               - 1\nC_H300               - 1\nC_H500               - 1\nC_H600               - 1\nsum: 4, p-value: 0.05366346, log2(fold\nchange): inf, edge multiplicities:\nA: 0 | 0\nC: 0 | 4\nG: 0 | 0\nT: 4 | 0\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, true));
    assert_eq!("[color=red, penwidth=8.469343, label=\"C: 4\"]", node.edge_dot_default(&colors, 1, debruijn::Dir::Right, true));

    let node_id = 27;
    let node = graph_tcpem.get_node(node_id);
    assert_eq!("[style=filled, color=black, penwidth=7, fillcolor=\"0 0.5 1\", fontcolor=black, label=\"id: 27, len: 16, exts: C|G, seq:\nCAGGTCCGTGGCCAGG\nsamples              - counts\nP_T400               - 1\nsum: 1, p-value: 0.3910022, log2(fold\nchange): -inf, edge multiplicities:\nA: 0 | 0\nC: 1 | 0\nG: 0 | 1\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, true));
    assert_eq!("[color=red, penwidth=-inf, label=\"C: 0\"]", node.edge_dot_default(&colors, 1, debruijn::Dir::Right, true));

    let node_id = 40;
    let node = graph_tcpem.get_node(node_id);
    assert_eq!("[style=filled, color=black, penwidth=7, fillcolor=\"0.15937972 0.5 1\", fontcolor=black, label=\"id: 40, len: 16, exts: A|AG, seq:\nGCTTTGCATCAGAAGA\nsamples              - counts\nC_H100               - 1\nC_H400               - 2\nP_T200               - 1\nP_T300               - 1\nP_T400               - 1\nsum: 6, p-value: 0.86424315,\nlog2(fold change): -0.19431013, edge\nmultiplicities:\nA: 6 | 3\nC: 0 | 0\nG: 0 | 3\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, true));
    assert_eq!("[color=red, penwidth=-inf, label=\"C: 0\"]", node.edge_dot_default(&colors, 1, debruijn::Dir::Right, true));
    
    let node_id = 2458;
    let node = graph_tcpem.get_node(node_id);
    assert_eq!("[style=filled, color=black, penwidth=7, fillcolor=\"0 1 1\", fontcolor=black, label=\"id: 2458, len: 16, exts: C|A, seq:\nTCTCGTACTAAGTTCA\nsamples              - counts\nC_H600               - 1\nP_T100               - 1\nP_T200               - 1\nP_T300               - 1\nP_T400               - 2\nsum: 6, p-value: 0.030437965,\nlog2(fold change): -4.4442472, edge\nmultiplicities:\nA: 0 | 6\nC: 6 | 0\nG: 0 | 0\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, true));
    assert_eq!("[color=red, penwidth=-inf, label=\"C: 0\"]", node.edge_dot_default(&colors, 1, debruijn::Dir::Right, true));

    // test with color mode SampleGroups

    let colors = Colors::new(&graph_tcpem, &config_tcpem, ColorMode::SampleGroups);

    let node_id = 0;
    let node = graph_tcpem.get_node(node_id);
    println!("{}", node.node_dot_default(&colors, &config_tcpem, &translator, false));

    let node_id = 40;
    let node = graph_tcpem.get_node(node_id);
    println!("{}", node.node_dot_default(&colors, &config_tcpem, &translator, false));

    let node_id = 128;
    let node = graph_tcpem.get_node(node_id);
    println!("{}", node.node_dot_default(&colors, &config_tcpem, &translator, false));

    // test with color mode IDs

    let colors_ids = Colors::new(&graph_ids, &config_ids, ColorMode::IDS { n_ids: translator.id_translator().as_ref().unwrap().len() });

    let node_id = 0;
    let node = graph_ids.get_node(node_id);
    println!("{}", node.node_dot_default(&colors_ids, &config_ids, &translator, false));

    let node_id = 40;
    let node = graph_ids.get_node(node_id);
    println!("{}", node.node_dot_default(&colors_ids, &config_ids, &translator, false));

    let node_id = 128;
    let node = graph_ids.get_node(node_id);
    println!("{}", node.node_dot_default(&colors_ids, &config_ids, &translator, false));

    // test with color mode IDGroups

    let id_ids = (0..14).zip(vec![1, 3, 0, 2, 4, 4, 4, 2, 1, 1, 2, 3, 0]).collect::<HashMap<ID, ID>>(); // 14 genes in test graph
    let colors_ids = Colors::new(&graph_ids, &config_ids, ColorMode::IDGroups { id_group_ids: &id_ids, n_id_groups: 5 });
    println!("{:?}", translator.id_translator().as_ref().unwrap());

    let node_id = 0;
    let node = graph_ids.get_node(node_id);
    println!("{}", node.node_dot_default(&colors_ids, &config_ids, &translator, false));

    let node_id = 40;
    let node = graph_ids.get_node(node_id);
    println!("{}", node.node_dot_default(&colors_ids, &config_ids, &translator, false));

    let node_id = 128;
    let node = graph_ids.get_node(node_id);
    println!("{}", node.node_dot_default(&colors_ids, &config_ids, &translator, false));

    // write node to dot

    graph_tcpem.to_dot(
        "test_dot.dot", 
        &|node| node.node_dot_default(&colors, &config_tcpem, &translator, false), 
        &|node, base, dir, flipped| node.edge_dot_default(&colors, base, dir, flipped)
    );

    graph_tcpem.to_dot_parallel(
        "test_dot_parallel.dot", 
        &|node| node.node_dot_default(&colors, &config_tcpem, &translator, false), 
        &|node, base, dir, flipped| node.edge_dot_default(&colors, base, dir, flipped)
    );

    graph_tcpem.to_dot_partial(
        "test_dot_partial.dot", 
        &|node| node.node_dot_default(&colors, &config_tcpem, &translator, false), 
        &|node, base, dir, flipped| node.edge_dot_default(&colors, base, dir, flipped),
        vec![0, 1, 2, 3]
    );

    remove_file("test_dot.dot").unwrap();
    remove_file("test_dot_parallel.dot").unwrap();
    remove_file("test_dot_partial.dot").unwrap();

    // write node to gfa

    graph_tcpem.to_gfa("test_gfa.gfa").unwrap();
    graph_tcpem.to_gfa_with_tags("test_gfa_tags.gfa", |node| node.data().print_ol(&translator, &config_tcpem)).unwrap();
    graph_tcpem.to_gfa_otags_parallel("test_gfa_parallel", Some(&|node: &graph::Node<_, TagsCountsPEMData>| node.data().print_ol(&translator, &config_tcpem))).unwrap();
    graph_tcpem.to_gfa_partial("test_gfa_partial.gfa", Some(&|node: &graph::Node<'_, debruijn::kmer::IntKmer<u32>, TagsCountsPEMData>| node.data().print_ol(&translator, &config_tcpem)), vec![0, 1, 2, 3]).unwrap();

    remove_file("test_gfa.gfa").unwrap();
    remove_file("test_gfa_tags.gfa").unwrap();
    remove_file("test_gfa_parallel.gfa").unwrap();
    remove_file("test_gfa_partial.gfa").unwrap();



}
