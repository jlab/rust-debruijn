use std::{collections::HashMap, fs::{remove_file, File}, io::BufReader};

use bimap::BiMap;
use debruijn::{colors::{ColorMode, Colors}, compression::{uncompressed_graph}, filter::filter_kmers, graph::{self, DebruijnGraph}, kmer::{Kmer16, Kmer8}, reads::{Reads, ReadsPaired}, serde::SerGraph, summarizer::{IDMapEMData, IDSumData, IDTag, SampleInfo, SummaryConfig, SummaryData, TagsCountsPEMData, Translator, ID}};

#[cfg(not(feature = "sample128"))]
const TEST_FILE_T: &str = "test_data/sided.graph.dbg";

// use different graph file for sample128 tests
#[cfg(feature = "sample128")]
const TEST_FILE_T: &str = "test_data/sided-128.graph.dbg";

// marbel dataset: marbel --n-orthogroups 5 --n-species 3 --n-samples 2 2 --library-size 200 --library-size-distribution negative_binomial 
#[cfg(not(feature = "id4b"))]
const TEST_FILE_IDS: &str = "test_data/test_graph_ids.graph.dbg";

#[cfg(feature = "id4b")]
const TEST_FILE_IDS: &str = "test_data/test_graph_ids-4b.graph.dbg";

#[test]
fn test_colors() {
    let file_t = BufReader::new(File::open(TEST_FILE_T).unwrap());

    let (graph_tcpem, vec_labels, config_tcpem): (DebruijnGraph<Kmer16, TagsCountsPEMData>, Vec<String>, SummaryConfig) = 
        bincode::deserialize_from(file_t).expect("error deserializing graph");
    let mut hashed_labels_tcpem: BiMap<String, u8> = BiMap::new();
    for (i, label) in vec_labels.iter().enumerate() {
        hashed_labels_tcpem.insert(label.clone(), i as u8);
    }

    let (graph_ids, translator, config_ids) = SerGraph::<Kmer16, IDSumData>::deserialize_from(TEST_FILE_IDS).dissolve();
    let (hashed_ids, _) = translator.dissolve();
    let translator  = Translator::new(hashed_ids.unwrap(), hashed_labels_tcpem);

    // test with color mode FoldChange

    let colors = Colors::new(&graph_tcpem, &config_tcpem, ColorMode::FoldChange);
    println!("colors:{:?}", colors);
    
    let node_id = 0;
    let node = graph_tcpem.get_node(node_id);
    assert_eq!("[style=filled, color=\"0.33333334 0.5 1\", fontcolor=black, label=\"id: 0, len: 16, exts: A|C, seq:\nCCAGCTGTCCCAGATA\nsamples              - counts\nC_H600               - 1\nsum: 1, p-value: 0.37390098, log2(fold\nchange): inf, edge coverage:\nA: 1 | 0\nC: 0 | 1\nG: 0 | 0\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, false, false));
    assert_eq!("[color=blue, penwidth=3, label=\"A: 1\"]", node.edge_dot_default(&colors, 0, debruijn::Dir::Left, true));

    let node_id = 3;
    let node = graph_tcpem.get_node(node_id);
    assert_eq!("[style=filled, color=black, penwidth=10, fillcolor=\"0.33333334 0.5 1\", fontcolor=black, label=\"id: 3, len: 16, exts: T|C, seq:\nCCGATACTTCTTCACG\nsamples              - counts\nC_H200               - 1\nC_H300               - 1\nC_H500               - 1\nC_H600               - 1\nsum: 4, p-value: 0.05366346, log2(fold\nchange): inf, edge coverage:\nA: 0 | 0\nC: 0 | 4\nG: 0 | 0\nT: 4 | 0\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, true, false));
    assert_eq!("[color=red, penwidth=9.683096, label=\"C: 4\"]", node.edge_dot_default(&colors, 1, debruijn::Dir::Right, true));

    let node_id = 27;
    let node = graph_tcpem.get_node(node_id);
    println!("{:?}", node.data());
    assert_eq!("[style=filled, color=black, penwidth=10, fillcolor=\"0 0.5 1\", fontcolor=black, label=\"id: 27, len: 16, exts: C|G, seq:\nCAGGTCCGTGGCCAGG\nsamples              - counts\nP_T400               - 1\nsum: 1, p-value: 0.3910022, log2(fold\nchange): -inf, edge coverage:\nA: 0 | 0\nC: 1 | 0\nG: 0 | 1\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, true, false));
    assert_eq!("[color=red, penwidth=3, label=\"G: 1\"]", node.edge_dot_default(&colors, 2, debruijn::Dir::Right, true));

    let node_id = 40;
    let node = graph_tcpem.get_node(node_id);
    println!("{:?}", node.data());
    assert_eq!("[style=filled, color=black, penwidth=10, fillcolor=\"0.15937972 0.5 1\", fontcolor=black, label=\"id: 40, len: 16, exts: A|AG, seq:\nGCTTTGCATCAGAAGA\nsamples              - counts\nC_H100               - 1\nC_H400               - 2\nP_T200               - 1\nP_T300               - 1\nP_T400               - 1\nsum: 6, p-value: 0.86424315, log2(fold\nchange): -0.19431013, edge coverage:\nA: 6 | 3\nC: 0 | 0\nG: 0 | 3\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, true, false));
    assert_eq!("[color=red, penwidth=8.296228, label=\"A: 3\"]", node.edge_dot_default(&colors, 0, debruijn::Dir::Right, true));
    
    let node_id = 2458;
    let node = graph_tcpem.get_node(node_id);
    println!("{:?}", node.data());
    assert_eq!("[style=filled, color=black, penwidth=10, fillcolor=\"0 1 1\", fontcolor=black, label=\"id: 2458, len: 16, exts: C|A, seq:\nTCTCGTACTAAGTTCA\nsamples              - counts\nC_H600               - 1\nP_T100               - 1\nP_T200               - 1\nP_T300               - 1\nP_T400               - 2\nsum: 6, p-value: 0.030437965, log2(fold\nchange): -4.4442472, edge coverage:\nA: 0 | 6\nC: 6 | 0\nG: 0 | 0\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, true, false));
    assert_eq!("[color=red, penwidth=11.637776, label=\"A: 6\"]", node.edge_dot_default(&colors, 0, debruijn::Dir::Right, true));

    // test with color mode SampleGroups

    let colors = Colors::new(&graph_tcpem, &config_tcpem, ColorMode::SampleGroups);

    let node_id = 0;
    let node = graph_tcpem.get_node(node_id);
    assert_eq!("[style=filled, color=\"0.33333334 1 1\", fontcolor=black, label=\"id: 0, len: 16, exts: A|C, seq:\nCCAGCTGTCCCAGATA\nsamples              - counts\nC_H600               - 1\nsum: 1, p-value: 0.37390098, log2(fold\nchange): inf, edge coverage:\nA: 1 | 0\nC: 0 | 1\nG: 0 | 0\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, false, false));

    let node_id = 40;
    let node = graph_tcpem.get_node(node_id);
    assert_eq!("[style=filled, color=\"0.16666667 1 1\", fontcolor=black, label=\"id: 40, len: 16, exts: A|AG, seq:\nGCTTTGCATCAGAAGA\nsamples              - counts\nC_H100               - 1\nC_H400               - 2\nP_T200               - 1\nP_T300               - 1\nP_T400               - 1\nsum: 6, p-value: 0.86424315, log2(fold\nchange): -0.19431013, edge coverage:\nA: 6 | 3\nC: 0 | 0\nG: 0 | 3\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, false, false));

    let node_id = 128;
    let node = graph_tcpem.get_node(node_id);
    assert_eq!("[style=filled, color=\"0 1 1\", fontcolor=black, label=\"id: 128, len: 16, exts: A|C, seq:\nAGAGATGGGGAAGCTC\nsamples              - counts\nP_T200               - 1\nP_T300               - 1\nP_T400               - 1\nsum: 3, p-value: 0.0825925, log2(fold\nchange): -inf, edge coverage:\nA: 3 | 0\nC: 0 | 3\nG: 0 | 0\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &config_tcpem, &translator, false, false));

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
        vec![0, 1, 2, 3]
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
    graph_tcpem.to_gfa_partial("test_gfa_partial.gfa", Some(&|node: &graph::Node<'_, debruijn::kmer::IntKmer<u32>, TagsCountsPEMData>| node.data().print_ol(&translator, &config_tcpem, None)), vec![0, 1, 2, 3]).unwrap();

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
    let summary_config = SummaryConfig::new(1, None, debruijn::summarizer::GroupFrac::None, 0.3, sample_info, None, debruijn::summarizer::StatTest::WelchsTTest);
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
    assert_eq!("[shape=rectangle, style=striped, color=\"0 1 1\", fontcolor=black, label=\"id: 0, len: 8, exts: |C, seq: AAAAAAAA\nIDs: ['ID A'], mapped IDs: []\"]", node.node_dot_default(&colors, &summary_config, &translator, false, false));

    let node = graph.get_node(1);
    assert_eq!("[shape=rectangle, style=striped, color=\"0 1 0.6\", penwidth=30, fillcolor=\"0 1 1\", fontcolor=black, label=\"id: 1, len: 8, exts: A|, seq: AAAAAAAC\nIDs: ['ID A'], mapped IDs: ['ID A']\"]", node.node_dot_default(&colors, &summary_config, &translator, false, false));

    // id_group_ids
    assert_eq!(&id_group_ids, colors.id_group_ids().unwrap());

    // color mode ids
    let colors = Colors::new(&graph, &summary_config, ColorMode::IDS { n_ids: 1 });

    let node = graph.get_node(0);
    assert_eq!("[shape=rectangle, style=striped, color=\"0 1 1\", fontcolor=black, label=\"id: 0, len: 8, exts: |C, seq: AAAAAAAA\nIDs: ['ID A'], mapped IDs: []\"]", node.node_dot_default(&colors, &summary_config, &translator, false, false));

    let node = graph.get_node(1);
    assert_eq!("[shape=rectangle, style=striped, color=\"0 1 0.6\", penwidth=30, fillcolor=\"0 1 1\", fontcolor=black, label=\"id: 1, len: 8, exts: A|, seq: AAAAAAAC\nIDs: ['ID A'], mapped IDs: ['ID A']\"]", node.node_dot_default(&colors, &summary_config, &translator, false, false));
    

}
