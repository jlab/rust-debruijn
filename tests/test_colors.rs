use std::{collections::HashMap, fs::remove_file};

use debruijn::{build_test_graph, colors::{ColorMode, Colors}, compression::uncompressed_graph, filter::filter_kmers, graph::{self}, kmer::{Kmer8, Kmer16}, reads::{Reads, ReadsPaired}, summarizer::{ID, IDMapEMData, IDTag, IDTagsCountsPEMData, SampleInfo, SummaryConfig, SummaryData, Translator}};

#[test]
fn test_colors() {
    let (_, _, ser_graph) = build_test_graph::<Kmer16, IDTagsCountsPEMData, _>(); 
    let (graph, translator, config) = ser_graph.dissolve();

    // with color mode FoldChange

    let colors= Colors::new(&graph, &config, ColorMode::FoldChange);
    println!("colors:{:?}", colors);
    
    let node_id = 0;
    let node = graph.get_node(node_id);
    assert_eq!("[style=filled, color=\"0.33333334 0.5 1\", fontcolor=black, label=\"id: 0, len: 16, exts: A|T, seq:\nCGTAGCGCAGGCACCA\nIDs: ['gene3'], samples              -\ncounts\nsample3              - 8\nsum: 8, p-value: 0.5, log2(fold change):\ninf, edge coverage:\nA: 8 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 7\n\"]", node.node_dot_default(&colors, &config, &translator, false, false));

    let node_id = 3;
    let node = graph.get_node(node_id);
    assert_eq!("[style=filled, color=black, penwidth=10, fillcolor=\"0 0.5 1\", fontcolor=black, label=\"id: 3, len: 21, exts: G|A, seq:\nCGGGGACGTATTATTATTAAA\nIDs: ['gene1'], samples              -\ncounts\nsample1              - 2\nsum: 2, p-value: 0.5, log2(fold change):\n-inf, edge coverage:\nA: 0 | 1\nC: 0 | 0\nG: 1 | 0\nT: 0 | 0\n\"]", node.node_dot_default(&colors, &config, &translator, true, false));
    
    let node_id = 57;
    let node = graph.get_node(node_id);
    println!("{:?}", node);
    assert_eq!("[style=filled, color=black, penwidth=10, fillcolor=\"0.21956564 0.5 1\", fontcolor=black, label=\"id: 57, len: 18, exts: T|T, seq:\nAGGCAGGACGCATTACTA\nIDs: ['gene2', 'gene4'], samples\n- counts\nsample2              - 3\nsample4              - 5\nsum: 8, p-value: 0.7896464, log2(fold\nchange): 0.7369656, edge coverage:\nA: 0 | 0\nC: 0 | 0\nG: 0 | 0\nT: 7 | 8\n\"]", node.node_dot_default(&colors, &config, &translator, true, false));
    assert_eq!("[color=red, penwidth=18.35253, label=\"T: 8\", weight=8]", node.edge_dot_default(&colors, 3, debruijn::Dir::Right, true));
    assert_eq!("[color=blue, penwidth=17.366667, label=\"T: 7\", weight=7]", node.edge_dot_default(&colors, 3, debruijn::Dir::Left, true));

    // with json
    assert_eq!("\"id\": 57, \"len\": 18, \"seq\": \"AGGCAGGACGCATTACTA\", \"hue\": 80, \"ids\": [\"gene2\", \"gene4\"], \"sum\": 8, \"samples\": [\"sample2\", \"sample4\"], \"counts\": [3, 5], \"p_value\": 0.7896464, \"fold_change\": 0.7369656", node.node_json_default(&colors, &config, &translator, true));
    assert_eq!("\"source\": 57, \"target\": 83, \"source_b\": \"T\", \"target_b\": \"G\", \"dir\": 0, \"strength\": 8", node.edge_json_default(83, 3, debruijn::Dir::Right, true));

    // test with color mode SampleGroups

    let colors = Colors::new(&graph, &config, ColorMode::SampleGroups);

    let node_id = 0;
    let node = graph.get_node(node_id);
    assert_eq!("[style=filled, color=\"0.33333334 1 1\", fontcolor=black, label=\"id: 0, len: 16, exts: A|T, seq:\nCGTAGCGCAGGCACCA\nIDs: ['gene3'], samples              -\ncounts\nsample3              - 8\nsum: 8, p-value: 0.5, log2(fold change):\ninf, edge coverage:\nA: 8 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 7\n\"]", node.node_dot_default(&colors, &config, &translator, false, false));

    let node_id = 57;
    let node = graph.get_node(node_id);
    assert_eq!("[style=filled, color=\"0.16666667 1 1\", fontcolor=black, label=\"id: 57, len: 18, exts: T|T, seq:\nAGGCAGGACGCATTACTA\nIDs: ['gene2', 'gene4'], samples\n- counts\nsample2              - 3\nsample4              - 5\nsum: 8, p-value: 0.7896464, log2(fold\nchange): 0.7369656, edge coverage:\nA: 0 | 0\nC: 0 | 0\nG: 0 | 0\nT: 7 | 8\n\"]", node.node_dot_default(&colors, &config, &translator, false, false));

    let node_id = 28;
    let node = graph.get_node(node_id);
    assert_eq!("[style=filled, color=\"0 1 1\", fontcolor=black, label=\"id: 28, len: 27, exts: T|A, seq:\nTATATATCATTATTTTTTTCTATAAAA\nIDs: ['gene2'], samples              -\ncounts\nsample2              - 5\nsum: 5, p-value: 0.5, log2(fold change):\n-inf, edge coverage:\nA: 0 | 4\nC: 0 | 0\nG: 0 | 0\nT: 5 | 0\n\"]", node.node_dot_default(&colors, &config, &translator, false, false));

    // with json
    assert_eq!("\"id\": 28, \"len\": 27, \"seq\": \"TATATATCATTATTTTTTTCTATAAAA\", \"hue\": 0, \"ids\": [\"gene2\"], \"sum\": 5, \"samples\": [\"sample2\"], \"counts\": [5], \"p_value\": 0.5, \"fold_change\": -inf", node.node_json_default(&colors, &config, &translator, false));

    // test with color mode IDs

    let colors_ids = Colors::new(&graph, &config, ColorMode::IDS { n_ids: translator.id_translator().as_ref().unwrap().len() });

    let node_id = 0;
    let node = graph.get_node(node_id);
    assert_eq!("[shape=rectangle, style=striped, color=\"0.5 1 1\", fontcolor=black, label=\"id: 0, len: 16, exts: A|T, seq:\nCGTAGCGCAGGCACCA\nIDs: ['gene3'], samples              -\ncounts\nsample3              - 8\nsum: 8, p-value: 0.5, log2(fold change):\ninf, edge coverage:\nA: 8 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 7\n\"]", node.node_dot_default(&colors_ids, &config, &translator, false, false));

    let node_id = 40;
    let node = graph.get_node(node_id);
    assert_eq!("[shape=rectangle, style=striped, color=\"0 1 1\", fontcolor=black, label=\"id: 40, len: 20, exts: A|G, seq:\nCGTATTATTATTAAAATTGC\nIDs: ['gene1'], samples              -\ncounts\nsample1              - 1\nsum: 1, p-value: 0.5, log2(fold change):\n-inf, edge coverage:\nA: 1 | 0\nC: 0 | 0\nG: 0 | 1\nT: 0 | 0\n\"]", node.node_dot_default(&colors_ids, &config, &translator, false, false));

    let node_id = 28;
    let node = graph.get_node(node_id);
    assert_eq!("[shape=rectangle, style=striped, color=\"0.25 1 1\", fontcolor=black, label=\"id: 28, len: 27, exts: T|A, seq:\nTATATATCATTATTTTTTTCTATAAAA\nIDs: ['gene2'], samples              -\ncounts\nsample2              - 5\nsum: 5, p-value: 0.5, log2(fold change):\n-inf, edge coverage:\nA: 0 | 4\nC: 0 | 0\nG: 0 | 0\nT: 5 | 0\n\"]", node.node_dot_default(&colors_ids, &config, &translator, false, false));

    let node_id = 36;
    let node = graph.get_node(node_id);
    assert_eq!("[shape=rectangle, style=striped, color=\"0.5 1 1\", fontcolor=black, label=\"id: 36, len: 20, exts: C|C, seq:\nTATATTACGCGATAAAGAGC\nIDs: ['gene3'], samples              -\ncounts\nsample3              - 6\nsum: 6, p-value: 0.5, log2(fold change):\ninf, edge coverage:\nA: 0 | 0\nC: 5 | 5\nG: 0 | 0\nT: 0 | 0\n\"]", node.node_dot_default(&colors_ids, &config, &translator, false, false));

    // with json
    assert_eq!("\"id\": 36, \"len\": 20, \"seq\": \"TATATTACGCGATAAAGAGC\", \"hue\": 121, \"ids\": [\"gene3\"], \"sum\": 6, \"samples\": [\"sample3\"], \"counts\": [6], \"p_value\": 0.5, \"fold_change\": inf", node.node_json_default(&colors, &config, &translator, false));

    // test with color mode IDGroups

    let id_ids = (0..4).zip(vec![0, 0, 1, 1]).collect::<HashMap<ID, ID>>(); // 4 genes in test graph
    let colors_ids = Colors::new(&graph, &config, ColorMode::IDGroups { id_group_ids: &id_ids, n_id_groups: 5 });
    //println!("{:?}", translator.id_translator().as_ref().unwrap());

    let node_id = 0;
    let node = graph.get_node(node_id);
    assert_eq!("[shape=rectangle, style=striped, color=\"0.2 1 1\", fontcolor=black, label=\"id: 0, len: 16, exts: A|T, seq:\nCGTAGCGCAGGCACCA\nIDs: [1], samples              - counts\nsample3              - 8\nsum: 8, p-value: 0.5, log2(fold change):\ninf, edge coverage:\nA: 8 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 7\n\"]", node.node_dot_default(&colors_ids, &config, &translator, false, true));

    let node_id = 40;
    let node = graph.get_node(node_id);
    assert_eq!("[shape=rectangle, style=striped, color=\"0 1 1\", fontcolor=black, label=\"id: 40, len: 20, exts: A|G, seq:\nCGTATTATTATTAAAATTGC\nIDs: [0], samples              - counts\nsample1              - 1\nsum: 1, p-value: 0.5, log2(fold change):\n-inf, edge coverage:\nA: 1 | 0\nC: 0 | 0\nG: 0 | 1\nT: 0 | 0\n\"]", node.node_dot_default(&colors_ids, &config, &translator, false, true));

    let node_id = 28;
    let node = graph.get_node(node_id);
    assert_eq!("[shape=rectangle, style=striped, color=\"0 1 1\", fontcolor=black, label=\"id: 28, len: 27, exts: T|A, seq:\nTATATATCATTATTTTTTTCTATAAAA\nIDs: [0], samples              - counts\nsample2              - 5\nsum: 5, p-value: 0.5, log2(fold change):\n-inf, edge coverage:\nA: 0 | 4\nC: 0 | 0\nG: 0 | 0\nT: 5 | 0\n\"]", node.node_dot_default(&colors_ids, &config, &translator, false, true));

    let node_id = 36;
    let node = graph.get_node(node_id);
    assert_eq!("[shape=rectangle, style=striped, color=\"0.2 1 1\", fontcolor=black, label=\"id: 36, len: 20, exts: C|C, seq:\nTATATTACGCGATAAAGAGC\nIDs: [1], samples              - counts\nsample3              - 6\nsum: 6, p-value: 0.5, log2(fold change):\ninf, edge coverage:\nA: 0 | 0\nC: 5 | 5\nG: 0 | 0\nT: 0 | 0\n\"]", node.node_dot_default(&colors_ids, &config, &translator, false, true));

    // with json
    assert_eq!("\"id\": 36, \"len\": 20, \"seq\": \"TATATTACGCGATAAAGAGC\", \"hue\": 121, \"ids\": [\"gene3\"], \"sum\": 6, \"samples\": [\"sample3\"], \"counts\": [6], \"p_value\": 0.5, \"fold_change\": inf", node.node_json_default(&colors, &config, &translator, false));
    
    // write node to dot

    graph.to_dot(
        "test_dot.dot", 
        &|node| node.node_dot_default(&colors, &config, &translator, false, false), 
        &|node, base, dir, flipped| node.edge_dot_default(&colors, base, dir, flipped)
    );

    graph.to_dot_parallel(
        "test_dot_parallel.dot", 
        &|node| node.node_dot_default(&colors, &config, &translator, false, false), 
        &|node, base, dir, flipped| node.edge_dot_default(&colors, base, dir, flipped)
    );

    graph.to_dot_partial(
        "test_dot_partial.dot", 
        &|node| node.node_dot_default(&colors, &config, &translator, false, false), 
        &|node, base, dir, flipped| node.edge_dot_default(&colors, base, dir, flipped),
        &[0, 1, 2, 3]
    );

    graph.to_dot_with_path(
        "test_dot_with_paths.dot", 
        &|node, base, dir, flipped| node.edge_dot_default(&colors, base, dir, flipped),
        &colors,
        &translator,
        &config,
        false
    );

    remove_file("test_dot.dot").unwrap();
    remove_file("test_dot_parallel.dot").unwrap();
    remove_file("test_dot_partial.dot").unwrap();
    remove_file("test_dot_with_paths.dot").unwrap();

    // write node to gfa

    graph.to_gfa("test_gfa.gfa").unwrap();
    graph.to_gfa_with_tags("test_gfa_tags.gfa", |node| node.data().print_ol(&translator, &config, None)).unwrap();
    graph.to_gfa_otags_parallel("test_gfa_parallel", Some(&|node: &graph::Node<_, IDTagsCountsPEMData>| node.data().print_ol(&translator, &config, None))).unwrap();
    graph.to_gfa_partial("test_gfa_partial.gfa", Some(&|node: &graph::Node<_, IDTagsCountsPEMData>| node.data().print_ol(&translator, &config, None)), vec![0, 1, 2, 3]).unwrap();

    remove_file("test_gfa.gfa").unwrap();
    remove_file("test_gfa_tags.gfa").unwrap();
    remove_file("test_gfa_parallel.gfa").unwrap();
    remove_file("test_gfa_partial.gfa").unwrap();

    // write to json (3d)
    graph.to_json_3d(
        "test_json_partial.json", 
        &|node| node.node_json_default(&colors, &config, &translator, false), 
        &|node, target_id, base, dir, flipped| node.edge_json_default(target_id, base, dir, flipped),
        Some(&vec![0, 1, 2, 3])
    ).unwrap();

    graph.to_json_3d(
        "test_json.json", 
        &|node| node.node_json_default(&colors, &config, &translator, false), 
        &|node, target_id, base, dir, flipped| node.edge_json_default(target_id, base, dir, flipped),
        None
    ).unwrap();

    remove_file("test_json.json").unwrap();
    remove_file("test_json_partial.json").unwrap();
}


#[test]
fn test_colors_mapped_ids() {
    let mut reads = Reads::new(debruijn::reads::Strandedness::Forward);
    reads.add_from_bytes("AAAAAAAAC".as_bytes(), None, IDTag::new(0, 0));
    let reads = ReadsPaired::Unpaired { reads };

    let sample_info = SampleInfo::new(0b1, 0b0, vec![2]);
    let summary_config = SummaryConfig::new(sample_info);
    let (kmers, _) = filter_kmers::<IDMapEMData, Kmer8, _>(&reads, &summary_config, false, 1., false);
    let mut graph = uncompressed_graph(kmers, true).finish();

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
