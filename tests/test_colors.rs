use std::{collections::HashMap, fs::remove_file};

use debruijn::{colors::{ColorMode, Colors}, compression::uncompressed_graph, filter::filter_kmers, graph::{self}, kmer::{Kmer8, Kmer16}, reads::{Reads, ReadsPaired}, serde::SerGraph, summarizer::{ID, IDMapEMData, IDSumData, IDTag, SampleInfo, SummaryConfig, SummaryData, TagsCountsPEMData, Translator}};

// cargo run --features low-ks -- -c data/test2.csv -o ../rust-debruijn/test_data/sided -s tags-counts-p-em -k 16
#[cfg(not(feature = "sample128"))]
const TEST_FILE_T: &str = "test_data/sided.graph.dbg";

// use different graph file for sample128 tests
// cargo run --features low-ks --features sample128 ---c data/test2.csv -o ../rust-debruijn/test_data/sided -s tags-counts-p-em -k 16
#[cfg(feature = "sample128")]
const TEST_FILE_T: &str = "test_data/sided-128.graph.dbg";

// marbel dataset: marbel --n-orthogroups 5 --n-species 3 --n-samples 2 2 --library-size 200 --library-size-distribution negative_binomial 
#[cfg(not(feature = "id4b"))]
const TEST_FILE_IDS: &str = "test_data/test_graph_ids.graph.dbg";

#[cfg(feature = "id4b")]
const TEST_FILE_IDS: &str = "test_data/test_graph_ids-4b.graph.dbg";

#[test]
fn test_colors() {

    let (graph_tcpem, translator_tcpem, config_tcpem) = SerGraph::<Kmer16, TagsCountsPEMData>::deserialize_from(TEST_FILE_T).dissolve();

    let hashed_labels_tcpem = translator_tcpem.tag_translator().as_ref().unwrap();

    let (graph_ids, translator, config_ids) = SerGraph::<Kmer16, IDSumData>::deserialize_from(TEST_FILE_IDS).dissolve();
    let (hashed_ids, _) = translator.dissolve();
    let translator  = Translator::new(hashed_ids.unwrap(), hashed_labels_tcpem.clone());

    // test with color mode FoldChange

    let colors = Colors::new(&graph_tcpem, &config_tcpem, ColorMode::FoldChange);
    println!("colors:{:?}", colors);
    
    let node_id = 0;
    let node = graph_tcpem.get_node(node_id);
    assert_eq!("[style=filled, color=\"0.33333334 0.5 1\", fontcolor=black, label=\"id: 0, len: 151, exts: |, seq:\nCCCCGTTTAATATCTCCCGCTCCATCTTAATGGTCGGATCCGGGATGGGAATGGCTGACAGAAGCGCCTGGGTATAGGGATGTAGCGGATGTTCAAACAGTTCATCCGGCTTGGCCCTCTCCACCAGCTGTCCCAGATACATAACCACAAT\nsamples              - counts\nC_H600               - 1\nsum: 1, p-value: 0.37390098, log2(fold change): inf, edge coverage:\nA: 0 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, false, false));

    let node_id = 3;
    let node = graph_tcpem.get_node(node_id);
    assert_eq!("[style=filled, color=black, penwidth=10, fillcolor=\"0.33333334 0.5 1\", fontcolor=black, label=\"id: 3, len: 150, exts: |, seq:\nCCCTGATATAAATCTGATATATTAGTCCGATACTTCTTCACGCATCCTCTAATATCTTTGATCTACGTGGTATTTATTCAGTATATCTGCAAATTTTAATCTCTTTACACCGGTGATTCTACTCTCCGATTTTTCCTTTATAAACTTTTA\nsamples              - counts\nC_H200               - 1\nC_H300               - 1\nC_H500               - 1\nC_H600               - 1\nsum: 4, p-value: 0.05366346, log2(fold change): inf, edge coverage:\nA: 0 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, true, false));
    
    let node_id = 456;
    let node = graph_tcpem.get_node(node_id);
    println!("{:?}", node.data());
    assert_eq!("[style=filled, color=black, penwidth=10, fillcolor=\"0.1132268 0.5 1\", fontcolor=black, label=\"id: 456, len: 17, exts: T|T, seq:\nCGCTGTTGTTCGATGAT\nsamples              - counts\nC_H200               - 1\nC_H300               - 1\nC_H500               - 1\nC_H600               - 1\nP_T100               - 1\nP_T200               - 1\nP_T300               - 1\nP_T400               - 1\nsum: 8, p-value: 0.14767084, log2(fold\nchange): -1.4249998, edge coverage:\nA: 0 | 0\nC: 0 | 0\nG: 0 | 0\nT: 8 | 4\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, true, false));
    assert_eq!("[color=red, penwidth=9.740156, label=\"T: 4\"]", node.edge_dot_default(&colors, 3, debruijn::Dir::Right, true));
    assert_eq!("[color=blue, penwidth=13.110233, label=\"T: 8\"]", node.edge_dot_default(&colors, 3, debruijn::Dir::Left, true));


    // test with color mode SampleGroups

    let colors = Colors::new(&graph_tcpem, &config_tcpem, ColorMode::SampleGroups);

    let node_id = 0;
    let node = graph_tcpem.get_node(node_id);
    assert_eq!("[style=filled, color=\"0.33333334 1 1\", fontcolor=black, label=\"id: 0, len: 151, exts: |, seq:\nCCCCGTTTAATATCTCCCGCTCCATCTTAATGGTCGGATCCGGGATGGGAATGGCTGACAGAAGCGCCTGGGTATAGGGATGTAGCGGATGTTCAAACAGTTCATCCGGCTTGGCCCTCTCCACCAGCTGTCCCAGATACATAACCACAAT\nsamples              - counts\nC_H600               - 1\nsum: 1, p-value: 0.37390098, log2(fold change): inf, edge coverage:\nA: 0 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, false, false));

    let node_id = 456;
    let node = graph_tcpem.get_node(node_id);
    assert_eq!("[style=filled, color=\"0.16666667 1 1\", fontcolor=black, label=\"id: 456, len: 17, exts: T|T, seq:\nCGCTGTTGTTCGATGAT\nsamples              - counts\nC_H200               - 1\nC_H300               - 1\nC_H500               - 1\nC_H600               - 1\nP_T100               - 1\nP_T200               - 1\nP_T300               - 1\nP_T400               - 1\nsum: 8, p-value: 0.14767084, log2(fold\nchange): -1.4249998, edge coverage:\nA: 0 | 0\nC: 0 | 0\nG: 0 | 0\nT: 8 | 4\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, false, false));

    let node_id = 128;
    let node = graph_tcpem.get_node(node_id);
    assert_eq!("[style=filled, color=\"0 1 1\", fontcolor=black, label=\"id: 128, len: 151, exts: |, seq:\nTTTACTTTTCAAGGAGTATTTCCTATGAACGAGTTAGACGGCATCAAACAGTTCACCACTGTCGTGGCAGACAGCGGCGATATTGAGTCCATTCGCCATTATCATCCCCAGGATGCCACCACCAATCCTTCGCTGTTACTCAAGGCTGCCG\nsamples              - counts\nP_T100               - 1\nP_T200               - 1\nP_T300               - 1\nP_T400               - 1\nsum: 4, p-value: 0.05123313, log2(fold change): -inf, edge coverage:\nA: 0 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, false, false));

    // test with color mode IDs

    let colors_ids = Colors::new(&graph_ids, &config_ids, ColorMode::IDS { n_ids: translator.id_translator().as_ref().unwrap().len() });

    let node_id = 0;
    let node = graph_ids.get_node(node_id);
    assert_eq!("[shape=rectangle, style=striped, color=\"0 1 1\", fontcolor=black, label=\"id: 0, len: 66, exts: C|A, seq:\nGCCGCCGCGACCCGCCGCGCGTGCCGCGCCTCCTCCAGCGCGCCGCGCAGCCCCTCCGCCGAGTGC\nIDs: ['SAM40697_RS27425'], sum: 5\"]", node.node_dot_default(&colors_ids, &config_ids, &translator, false, false));

    let node_id = 40;
    let node = graph_ids.get_node(node_id);
    assert_eq!("[shape=rectangle, style=striped, color=\"0 1 1\", fontcolor=black, label=\"id: 40, len: 68, exts: C|G, seq:\nTGGGTGCTCTCGCCGACCGGGCGGCCGGTCGCGGGCCCCAAGGACGCGGGTCCCGTGCTGCCGTCCGA\nIDs: ['SAM40697_RS27425'], sum: 5\"]", node.node_dot_default(&colors_ids, &config_ids, &translator, false, false));

    let node_id = 128;
    let node = graph_ids.get_node(node_id);
    assert_eq!("[shape=rectangle, style=striped, color=\"0.13333334 1 1\", fontcolor=black, label=\"id: 128, len: 31, exts: A|T, seq:\nCCGGTCCGCGGCCACAGCGTGCAGGTCGCGC\nIDs: ['SAM40697_RS13660'], sum: 1\"]", node.node_dot_default(&colors_ids, &config_ids, &translator, false, false));

    let node_id = 36;
    let node = graph_ids.get_node(node_id);
    assert_eq!("[shape=rectangle, style=striped, color=\"0 1 1:0.33333334 1 1:0.6666667 1 1\", fontcolor=black, label=\"id: 36, len: 22, exts: AC|CG, seq:\nTGTCGCGCCTGGAGGACAAGCT\nIDs: ['SAM40697_RS27425',\n'IE258_RS26570', 'CP976_RS34340'], sum:\n4\"]", node.node_dot_default(&colors_ids, &config_ids, &translator, false, false));

    // test with color mode IDGroups

    let id_ids = (0..15).zip(vec![0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4]).collect::<HashMap<ID, ID>>(); // 15 genes in test graph
    let colors_ids = Colors::new(&graph_ids, &config_ids, ColorMode::IDGroups { id_group_ids: &id_ids, n_id_groups: 5 });
    //println!("{:?}", translator.id_translator().as_ref().unwrap());

    let node_id = 0;
    let node = graph_ids.get_node(node_id);
    assert_eq!("[shape=rectangle, style=striped, color=\"0 1 1\", fontcolor=black, label=\"id: 0, len: 66, exts: C|A, seq:\nGCCGCCGCGACCCGCCGCGCGTGCCGCGCCTCCTCCAGCGCGCCGCGCAGCCCCTCCGCCGAGTGC\nIDs: ['SAM40697_RS27425'], sum: 5\"]", node.node_dot_default(&colors_ids, &config_ids, &translator, false, false));

    let node_id = 40;
    let node = graph_ids.get_node(node_id);
    assert_eq!("[shape=rectangle, style=striped, color=\"0 1 1\", fontcolor=black, label=\"id: 40, len: 68, exts: C|G, seq:\nTGGGTGCTCTCGCCGACCGGGCGGCCGGTCGCGGGCCCCAAGGACGCGGGTCCCGTGCTGCCGTCCGA\nIDs: ['SAM40697_RS27425'], sum: 5\"]", node.node_dot_default(&colors_ids, &config_ids, &translator, false, false));

    let node_id = 128;
    let node = graph_ids.get_node(node_id);
    assert_eq!("[shape=rectangle, style=striped, color=\"0.4 1 1\", fontcolor=black, label=\"id: 128, len: 31, exts: A|T, seq:\nCCGGTCCGCGGCCACAGCGTGCAGGTCGCGC\nIDs: ['SAM40697_RS13660'], sum: 1\"]", node.node_dot_default(&colors_ids, &config_ids, &translator, false, false));

    let node_id = 36;
    let node = graph_ids.get_node(node_id);
    assert_eq!("[shape=rectangle, style=striped, color=\"0 1 1:0 1 1:0 1 1\", fontcolor=black, label=\"id: 36, len: 22, exts: AC|CG, seq:\nTGTCGCGCCTGGAGGACAAGCT\nIDs: ['SAM40697_RS27425',\n'IE258_RS26570', 'CP976_RS34340'], sum:\n4\"]", node.node_dot_default(&colors_ids, &config_ids, &translator, false, false));

    // write node to dot

    graph_tcpem.to_dot(
        "test_dot.dot", 
        &|node| node.node_dot_default(&colors, &config_tcpem, &translator, false, false), 
        &|node, base, dir, flipped| node.edge_dot_default(&colors, base, dir, flipped)
    );

    graph_tcpem.to_dot_parallel(
        "test_dot_parallel.dot", 
        &|node| node.node_dot_default(&colors, &config_tcpem, &translator, false, false), 
        &|node, base, dir, flipped| node.edge_dot_default(&colors, base, dir, flipped)
    );

    graph_tcpem.to_dot_partial(
        "test_dot_partial.dot", 
        &|node| node.node_dot_default(&colors, &config_tcpem, &translator, false, false), 
        &|node, base, dir, flipped| node.edge_dot_default(&colors, base, dir, flipped),
        &[0, 1, 2, 3]
    );

    graph_tcpem.to_dot_with_path(
        "test_dot_with_paths.dot", 
        &|node, base, dir, flipped| node.edge_dot_default(&colors, base, dir, flipped),
        &colors,
        &translator,
        &config_tcpem,
        false
    );

    remove_file("test_dot.dot").unwrap();
    remove_file("test_dot_parallel.dot").unwrap();
    remove_file("test_dot_partial.dot").unwrap();
    remove_file("test_dot_with_paths.dot").unwrap();

    // write node to gfa

    graph_tcpem.to_gfa("test_gfa.gfa").unwrap();
    graph_tcpem.to_gfa_with_tags("test_gfa_tags.gfa", |node| node.data().print_ol(&translator, &config_tcpem, None)).unwrap();
    graph_tcpem.to_gfa_otags_parallel("test_gfa_parallel", Some(&|node: &graph::Node<_, TagsCountsPEMData>| node.data().print_ol(&translator, &config_tcpem, None))).unwrap();
    graph_tcpem.to_gfa_partial("test_gfa_partial.gfa", Some(&|node: &graph::Node<_, TagsCountsPEMData>| node.data().print_ol(&translator, &config_tcpem, None)), vec![0, 1, 2, 3]).unwrap();

    remove_file("test_gfa.gfa").unwrap();
    remove_file("test_gfa_tags.gfa").unwrap();
    remove_file("test_gfa_parallel.gfa").unwrap();
    remove_file("test_gfa_partial.gfa").unwrap();
}


#[test]
fn test_colors_mapped_ids() {
    let mut reads = Reads::new(debruijn::reads::Strandedness::Forward);
    reads.add_from_bytes("AAAAAAAAC".as_bytes(), None, IDTag::new(0, 0));
    let reads = ReadsPaired::Unpaired { reads };

    let sample_info = SampleInfo::new(0b1, 0b0, vec![2]);
    let summary_config = SummaryConfig::new(sample_info);
    let (kmers, _) = filter_kmers::<IDMapEMData, Kmer8, _>(&reads, &summary_config, false, 1., false);
    let mut graph = uncompressed_graph(&kmers, true).finish();

    graph.print();

    graph.mut_data(1).set_mapped_ids(vec![0].into());

    let id_group_ids = [(0, 0)].into_iter().collect();
    let id_translator = [("ID A".to_string(), 0)].into_iter().collect();
    let tag_translator = [("TAG B".to_string(), 0)].into_iter().collect();
    let translator = Translator::new(id_translator, tag_translator);

    // color mode id groups
    let colors = Colors::new(&graph, &summary_config, ColorMode::IDGroups { id_group_ids: &id_group_ids, n_id_groups: 1 });

    let node = graph.get_node(0);
    assert_eq!("[shape=rectangle, style=striped, color=\"0 1 1\", fontcolor=black, label=\"id: 0, len: 8, exts: |C, seq: AAAAAAAA\nIDs: ['ID A'], mapped IDs: [], edge\ncoverage: A: 0 | 0\nC: 0 | 1\nG: 0 | 0\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &summary_config, &translator, false, false));

    let node = graph.get_node(1);
    assert_eq!("[shape=rectangle, style=striped, color=\"0 1 0.6\", penwidth=30, fillcolor=\"0 1 1\", fontcolor=black, label=\"id: 1, len: 8, exts: A|, seq: AAAAAAAC\nIDs: ['ID A'], mapped IDs: ['ID A'],\nedge coverage: A: 1 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &summary_config, &translator, false, false));

    // id_group_ids
    assert_eq!(&id_group_ids, colors.id_group_ids().unwrap());

    // color mode ids
    let colors = Colors::new(&graph, &summary_config, ColorMode::IDS { n_ids: 1 });

    let node = graph.get_node(0);
    assert_eq!("[shape=rectangle, style=striped, color=\"0 1 1\", fontcolor=black, label=\"id: 0, len: 8, exts: |C, seq: AAAAAAAA\nIDs: ['ID A'], mapped IDs: [], edge\ncoverage: A: 0 | 0\nC: 0 | 1\nG: 0 | 0\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &summary_config, &translator, false, false));

    let node = graph.get_node(1);
    assert_eq!("[shape=rectangle, style=striped, color=\"0 1 0.6\", penwidth=30, fillcolor=\"0 1 1\", fontcolor=black, label=\"id: 1, len: 8, exts: A|, seq: AAAAAAAC\nIDs: ['ID A'], mapped IDs: ['ID A'],\nedge coverage: A: 1 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &summary_config, &translator, false, false));
    

}
