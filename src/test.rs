// Copyright 2017 10x Genomics

//! Generate random genomes (with lots of re-used sustrings), reassemble them, and check sanity

use crate::complement;
use crate::Kmer;
use crate::Vmer;

use rand::distributions::{Distribution, Gamma, Range};
use rand::{self, Rng, RngCore};
use std::cmp::{max, min};

/// Generate a uniformly random base
pub fn random_base(r: &mut impl Rng) -> u8 {
    (r.next_u64() % 4) as u8
}

/// Generate uniformly random DNA sequences
pub fn random_dna(len: usize) -> Vec<u8> {
    let mut r = rand::thread_rng();
    (0..len)
        .map(|_| (r.next_u64() % 4) as u8)
        .collect()
}

/// Randomly mutate each base with probability `p`
pub fn edit_dna<R: Rng>(seq: &mut [u8], p: f64, r: &mut R) {
    for b in seq.iter_mut() {
        if r.gen_range(0.0, 1.0) < p {
            *b = random_base(r);
        }
    }
}

pub fn random_kmer<K: Kmer>() -> K {
    let mut r = rand::thread_rng();
    let mut kmer = K::empty();
    for pos in 0..K::k() {
        let b = (r.next_u64() % 4) as u8;
        kmer.set_mut(pos, b);
    }
    kmer
}

pub fn random_vmer<K: Kmer, V: Vmer>() -> V {
    let mut r = rand::thread_rng();
    let len = r.gen_range(K::k(), min(200, V::max_len()));
    let mut lmer = V::new(len);

    for pos in 0..len {
        let b = (r.next_u64() % 4) as u8;
        lmer.set_mut(pos, b);
    }
    lmer
}

pub fn simple_random_contigs() -> Vec<Vec<u8>> {
    let p1 = random_dna(40);
    let p2 = random_dna(30);

    let pc = random_dna(100);

    let p3 = random_dna(30);
    let p4 = random_dna(40);

    // Simulate contigs
    let mut c1 = Vec::new();
    c1.extend(p1);
    c1.extend(pc.clone());
    c1.extend(p3);

    let mut c2 = Vec::new();
    c2.extend(p2);
    c2.extend(pc);
    c2.extend(p4);

    let mut c3 = Vec::new();
    c3.extend(random_dna(30));

    // Stick a palindrome in the middle
    let palindrome1 = random_dna(33);
    let mut palindrome2 = palindrome1.clone();
    palindrome2.reverse();
    for v in &mut palindrome2 {
        *v = complement(*v);
    }

    c3.extend(palindrome1);
    c3.extend(palindrome2);
    c3.extend(random_dna(50));

    let contigs = vec![c1, c2, c3];
    contigs
}

// Generate random contigs with complicated repeats
pub fn random_contigs() -> Vec<Vec<u8>> {
    // Generate a bunch of sequence chunks
    let mut rng = rand::thread_rng();

    let gamma_dist = Gamma::new(0.6, 25.0);

    let nchunks = max(5, gamma_dist.sample(&mut rng) as u32);
    let chunk_sample = Range::new(0, nchunks);

    let length_dist = Gamma::new(1.5, 200.0);

    let mut chunks: Vec<Vec<u8>> = Vec::with_capacity(nchunks as usize);
    for _ in 0..nchunks {
        let len = max(10, length_dist.sample(&mut rng) as usize);
        let seq = random_dna(len);
        chunks.push(seq);
    }

    // Now make a bunch of chromosomes by pasting together chunks
    let nchrom = max(4, gamma_dist.sample(&mut rng) as u32);

    let mut chroms = Vec::with_capacity(nchrom as usize);
    for _ in 0..nchrom {
        let chrom_chunks = max(4, gamma_dist.sample(&mut rng) as u32);

        let mut chrom_seq = Vec::new();
        for _ in 0..chrom_chunks {
            let chunk_idx = chunk_sample.sample(&mut rng) as usize;
            chrom_seq.extend(chunks[chunk_idx].clone());
        }
        chroms.push(chrom_seq);
    }

    chroms
}

#[cfg(test)]
mod tests {

    use crate::compression::{compress_graph, compress_kmers_with_hash, ScmapCompress, SimpleCompress};
    use crate::graph::{self, BaseGraph};
    use crate::reads::{Reads, Strandedness};
    use crate::{DnaBytes, Tags};
    use crate::{Dir, Exts, Kmer};
    use bimap::BiMap;
    use boomphf::hashmap::BoomHashMap2;
    use boomphf::Mphf;
    use std::collections::{HashMap, HashSet};
    use std::fs::File;
    use std::iter::FromIterator;
    use std::marker::PhantomData;
    use std::time;

    use crate::dna_string::DnaString;
    use crate::filter::{self, filter_kmers};
    use crate::kmer::Kmer6;
    use crate::kmer::{IntKmer, VarIntKmer, K31};
    use crate::msp;
    use std::ops::Sub;
    use crate::summarizer::{GroupFrac, SampleInfo, SummaryConfig, TagsCountsEMData, TagsCountsSumData, TagsSumData};

    use super::*;
    use pretty_assertions::assert_eq;

    #[test]
    fn simple_kmer_compress() {
        let contigs = simple_random_contigs();
        reassemble_contigs::<IntKmer<u64>, DnaString>(contigs, false);
    }

    #[test]
    fn simple_sharded() {
        let contigs = simple_random_contigs();
        reassemble_sharded::<IntKmer<u64>, DnaString>(contigs, false);
    }

    #[test]
    fn degen_seq_asm() {
        let ctg = "AAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAA";
        let seq: Vec<u8> = ctg
            .as_bytes()
            .iter()
            .cloned()
            .map(crate::base_to_bits)
            .collect();

        reassemble_contigs::<VarIntKmer<u64, K31>, DnaString>(vec![seq.clone(), seq], false);
    }

    #[test]
    fn degen_seq_asm_sharded() {
        let ctg = "AAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAA";
        let seq: Vec<u8> = ctg
            .as_bytes()
            .iter()
            .cloned()
            .map(crate::base_to_bits)
            .collect();

        reassemble_sharded::<VarIntKmer<u64, K31>, DnaString>(vec![seq.clone(), seq], false);
    }

    #[test]
    fn complex_kmer_compress() {
        for _ in 0..10 {
            let contigs = random_contigs();
            reassemble_contigs::<IntKmer<u64>, DnaString>(contigs, false);
        }
    }

    #[test]
    fn complex_sharded() {
        for _ in 0..10 {
            let contigs = random_contigs();
            reassemble_sharded::<IntKmer<u64>, DnaString>(contigs, false);
        }
    }

    #[test]
    fn simple_path_compress() {
        let contigs = simple_random_contigs();
        simplify_from_kmers::<IntKmer<u64>>(contigs, false);
    }

    #[test]
    fn complex_path_compress_k31() {
        for _ in 0..100 {
            let contigs = random_contigs();
            simplify_from_kmers::<VarIntKmer<u64, K31>>(contigs, false);
        }
    }

    #[test]
    fn complex_path_compress() {
        for _ in 0..10 {
            let contigs = random_contigs();
            simplify_from_kmers::<IntKmer<u64>>(contigs, false);
        }
    }

    fn simplify_from_kmers<K: Kmer + Send + Sync>(mut contigs: Vec<Vec<u8>>, stranded: bool) {
        let seqs: Vec<(DnaBytes, Exts, ())> = contigs
            .drain(..)
            .map(|x| (DnaBytes(x), Exts::empty(), ()))
            .collect();


        let sample_info = SampleInfo::new(0, 0, 0, 0, Vec::new());
        let config = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info, None, crate::summarizer::StatTest::StudentsTTest);


        let (valid_kmers, _): (BoomHashMap2<K, Exts, u32>, _) = filter::filter_kmers(
            &crate::reads::ReadsPaired::Unpaired { reads: Reads::from_vmer_vec(seqs, Strandedness::Unstranded) },
            &config,
            stranded,
            4,
            true,
        );

        let spec =
            SimpleCompress::new(|d1: u32, d2: &u32| ((d1 + *d2) % 65535));
        let from_kmers = compress_kmers_with_hash(stranded, &spec, &valid_kmers, true, false, true).finish();
        let is_cmp = from_kmers.is_compressed(&spec);
        if is_cmp.is_some() {
            println!("not compressed: nodes: {:?}", is_cmp);
            from_kmers.print();
        }
        assert!(from_kmers.is_compressed(&spec).is_none());

        // Create a DBG with one node per input kmer
        let mut base_graph: BaseGraph<K, u16> = BaseGraph::new(stranded);

        for (kmer, exts, _) in valid_kmers.iter() {
            base_graph.add(kmer.iter(), *exts, 1);
        }
        let uncompressed_dbg = base_graph.finish();

        // Canonicalize the graph with
        let spec = SimpleCompress::new(|d1: u16, d2: &u16| d1 + d2);
        let simp_dbg = compress_graph(stranded, &spec, uncompressed_dbg, None);

        let is_cmp = simp_dbg.is_compressed(&spec);
        if is_cmp.is_some() {
            println!("not compressed: nodes: {:?}", is_cmp);
            simp_dbg.print();
        }

        assert!(simp_dbg.is_compressed(&spec).is_none());

        let total_kmers = valid_kmers.len();

        // Test the Boomphf DBG indexing machinery
        // Make an MPHF of the kmers in the DBG.
        // Each kmer should hash to a unique slot.
        let mphf =
            Mphf::from_chunked_iterator_parallel(1.7, &simp_dbg, None, valid_kmers.len() as u64, 2);

        let mut got_slot = vec![false; total_kmers];

        for n in simp_dbg.iter_nodes() {
            for kmer in n.sequence().iter_kmers() {
                let r = mphf.try_hash(&kmer).unwrap() as usize;
                assert!(!got_slot[r]);
                got_slot[r] = true;
            }
        }

        assert!(got_slot.iter().all(|x| *x));
    }


    // Take some input contig, which likely form a complicated graph,
    // and test the kmer, bsp, sedge and edge construction machinery
    fn reassemble_contigs<K: Kmer + Copy + Send + Sync, V: Vmer + Clone>(contigs: Vec<Vec<u8>>, stranded: bool) {
        let ctg_lens: Vec<_> = contigs.iter().map(std::vec::Vec::len).collect();
        println!("Reassembling contig sizes: {:?}", ctg_lens);

        // kmer vector
        let mut kmers = Vec::new();
        for c in contigs.iter() {
            let mut _kmers = K::kmers_from_bytes(c);
            kmers.extend(_kmers.iter().map(Kmer::min_rc));
        }

        // True kmer set
        let mut kmer_set = HashSet::new();
        kmer_set.extend(kmers.iter());

        // Bsps of kmers
        let p = 6;

        let mut seqs: Vec<(V, Exts, u8)> = Vec::new();
        let permutation: Vec<usize> = (0..1 << (2 * p)).collect();
        for c in contigs.iter() {
            let msps =
                msp::msp_sequence::<Kmer6, V>(K::k(), c.as_slice(), Some(&permutation), true);
            seqs.extend(msps.clone().into_iter().map(|(_, e, v)| (v, e, 0u8)));
            seqs.extend(msps.into_iter().map(|(_, e, v)| (v, e, 1u8)));
        }

        // kmer set from bsps
        let mut msp_kmers = HashSet::new();
        for (v, _, _) in seqs.iter() {
            for k in v.iter_kmers::<K>() {
                msp_kmers.insert(k.min_rc());
            }
        }

        // Kmer extensions from BSP match raw kmers
        if kmer_set != msp_kmers {
            println!("{:?}", kmer_set);
            println!("{:?}", msp_kmers);
        }

        // Raw kmers and BSP kmers match
        assert!(kmer_set == msp_kmers);

        let sample_info = SampleInfo::new(0, 0, 0, 0,Vec::new());
        let config = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info, None, crate::summarizer::StatTest::StudentsTTest);



        // Check the correctness of the process_kmer_shard kmer filtering function
        let (valid_kmers, _): (BoomHashMap2<K, Exts, u32>, _) = filter::filter_kmers(
            &crate::reads::ReadsPaired::Unpaired { reads: Reads::from_vmer_vec(seqs, Strandedness::Unstranded) },
            &config,
            stranded,
            4,
            true,
        );
        let mut process_kmer_set: HashSet<K> = HashSet::new();
        for k in valid_kmers.iter().map(|x| x.0) {
            process_kmer_set.insert(*k);
        }
        assert_eq!(process_kmer_set, kmer_set);

        // Every kmer should be reachable as the extension of a kmer.
        // No new kmers should be reachable
        let mut extension_kmer_set: HashSet<K> = HashSet::new();
        for (kmer, exts, _) in &valid_kmers {
            for e in kmer.get_extensions(*exts, Dir::Left) {
                extension_kmer_set.insert(e.min_rc());
            }

            for e in kmer.get_extensions(*exts, Dir::Right) {
                extension_kmer_set.insert(e.min_rc());
            }
        }

        // Kmer extensions from BSP match raw kmers
        if kmer_set != extension_kmer_set {
            println!("n:{}, {:?}", kmer_set.len(), kmer_set);
            println!("n:{}, {:?}", extension_kmer_set.len(), extension_kmer_set);

            if extension_kmer_set.is_superset(&kmer_set) {
                let invented = extension_kmer_set.sub(&kmer_set);
                println!("Invented kmers: {:?}", invented);
            }
        }

        assert!(kmer_set == extension_kmer_set);

        let spec = SimpleCompress::new(|d1: u32, d2: &u32| d1.saturating_add(*d2));

        // Generate compress DBG for these kmers
        let graph = compress_kmers_with_hash(stranded, &spec, &valid_kmers, true, false, true);

        // Check that all the lines have valid kmers,
        // and have extensions into other valid kmers
        let mut all_contig_set = HashSet::new();

        for seq_id in 0..graph.len() {
            let seq = graph.sequences.get(seq_id);
            let exts = graph.exts[seq_id];

            let contig_set = HashSet::from_iter(seq.iter_kmers().map(|km: K| km.min_rc()));
            all_contig_set.extend(contig_set.iter());
            assert!(kmer_set.is_superset(&contig_set));

            for l_ext in exts.get(Dir::Left) {
                let f: K = seq.first_kmer();
                let ext_kmer: K = f.extend_left(l_ext);
                assert!(kmer_set.contains(&(ext_kmer.min_rc())));
            }

            for r_ext in exts.get(Dir::Right) {
                let f: K = seq.last_kmer();
                let ext_kmer: K = f.extend_right(r_ext);
                assert!(kmer_set.contains(&(ext_kmer.min_rc())));
            }
        }

        assert_eq!(kmer_set, all_contig_set);
    }

    // Take some input contig, which likely form a complicated graph,
    // and the msp / shard_asm / main_asm loop
    fn reassemble_sharded<K: Kmer + Copy + Sync + Send, V: Vmer + Clone>(
        contigs: Vec<Vec<u8>>,
        stranded: bool,
    ) {
        let ctg_lens: Vec<_> = contigs.iter().map(std::vec::Vec::len).collect();
        println!("Reassembling contig sizes: {:?}", ctg_lens);

        // kmer vector
        let mut kmer_set = HashSet::new();
        for c in contigs.iter() {
            let mut _kmers = K::kmers_from_bytes(c);
            kmer_set.extend(_kmers.iter().map(Kmer::min_rc));
        }

        // Bsps of kmers
        let mut shards = HashMap::new();

        for ctg in contigs.iter() {
            let msps = msp::msp_sequence::<Kmer6, V>(K::k(), ctg.as_slice(), None, true);

            for (shard, exts, seq) in msps {
                let shard_vec = shards.entry(shard).or_insert_with(Vec::new);
                shard_vec.push((seq.clone(), exts, 0u8));
                shard_vec.push((seq, exts, 1u8));
            }
        }

        let mut shard_asms = Vec::new();

        let sample_info = SampleInfo::new(0, 0, 0, 0,Vec::new());
        let config = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), None, crate::summarizer::StatTest::StudentsTTest);



        // Do a subassembly in each shard
        for seqs in shards.into_values() {
            // Check the correctness of the process_kmer_shard kmer filtering function
            let (valid_kmers, _): (BoomHashMap2<K, Exts, u32>, _) = filter::filter_kmers(
                &crate::reads::ReadsPaired::Unpaired { reads: Reads::from_vmer_vec(seqs, Strandedness::Unstranded) },
                &config,
                stranded,
                4,
                true,
            );

            // Generate compress DBG for this shard
            let spec = SimpleCompress::new(|d1: u32, d2: &u32| d1.saturating_add(*d2));

            //print!("{:?}", valid_kmers);
            let graph = compress_kmers_with_hash(stranded, &spec, &valid_kmers, true, false, true);
            shard_asms.push(graph.clone());
            //graph.finish().print();
        }

        // Shove the subassemblies into a partially compress base graph
        let combined_graph = BaseGraph::combine(shard_asms.into_iter()).finish();
        let cmp = SimpleCompress::new(|a: u32, b: &u32| max(a, *b));
        let dbg_graph = compress_graph(false, &cmp, combined_graph, None);

        // Switch on for debugging
        //dbg_graph.print();
        //dbg_graph.write_gfa(&mut std::io::stdout().lock());

        let graph = dbg_graph.base;

        // Check that all the lines have valid kmers,
        // and have extensions into other valid kmers
        let mut all_contig_set = HashSet::new();

        for seq_id in 0..graph.len() {
            let seq = graph.sequences.get(seq_id);
            let exts = graph.exts[seq_id];

            let contig_set = HashSet::from_iter(seq.iter_kmers().map(|km: K| km.min_rc()));
            all_contig_set.extend(contig_set.iter());
            assert!(kmer_set.is_superset(&contig_set));

            for l_ext in exts.get(Dir::Left) {
                let f: K = seq.first_kmer();
                let ext_kmer: K = f.extend_left(l_ext);
                assert!(kmer_set.contains(&(ext_kmer.min_rc())));
            }

            for r_ext in exts.get(Dir::Right) {
                let f: K = seq.last_kmer();
                let ext_kmer: K = f.extend_right(r_ext);
                assert!(kmer_set.contains(&(ext_kmer.min_rc())));
            }
        }

        assert_eq!(kmer_set, all_contig_set);
    }

    #[test]
    fn simple_tip_clean() {
        let contigs = vec![random_dna(100), random_dna(100)];
        test_tip_cleaning::<IntKmer<u64>>(contigs, false);
    }

    // Take some input contig, which likely form a complicated graph,
    // and test the kmer, bsp, sedge and edge construction machinery
    fn test_tip_cleaning<K: Kmer + Sync + Send>(contigs: Vec<Vec<u8>>, stranded: bool) {

        let mut rng = rand::thread_rng();
        
        let mut clean_seqs = Reads::new(Strandedness::Unstranded);
        let mut all_seqs = Reads::new(Strandedness::Unstranded);

        // Generate 5x coverage of the main sequences & add some tips
        for c in contigs {
            if c.len() < K::k() * 3 {
                continue;
            }

            for _i in 0..5 {
                let read = DnaString::from_bytes(&c);
                clean_seqs.add_read(read.clone(), Exts::empty(), rng.gen_range(0, 10) as u8);
                all_seqs.add_read(read, Exts::empty(), rng.gen_range(0, 10) as u8);
            }

            let junk = random_dna(5);
            let mut err_ctg = c.clone();
            let l = err_ctg.len();
            err_ctg.truncate(l / 2);
            err_ctg.extend(junk);
            all_seqs.add_read(DnaString::from_bytes(&err_ctg), Exts::empty(), 3u8);
            all_seqs.add_read(DnaString::from_bytes(&err_ctg), Exts::empty(), 3u8);
        }

        let sample_info = SampleInfo::new(0, 0, 0, 0,Vec::new());
        let config = SummaryConfig::new(2, None, GroupFrac::None, 0.33, sample_info.clone(), None, crate::summarizer::StatTest::StudentsTTest);

        // Assemble w/o tips
        let (valid_kmers_clean, _): (BoomHashMap2<K, Exts, u32>, _) = filter::filter_kmers(
            &crate::reads::ReadsPaired::Unpaired { reads: clean_seqs },
            &config,
            stranded,
            4,
            true,
        );
        let spec = SimpleCompress::new(|d1: u32, d2: &u32| d1 + d2);
        let graph = compress_kmers_with_hash(stranded, &spec, &valid_kmers_clean, true, false, true);
        let graph1 = graph.finish();
        graph1.print();
        println!("components: {:?}", graph1.components_r());

        let all_seqs_p = crate::reads::ReadsPaired::Unpaired { reads: all_seqs };
        // Assemble w/ tips
        let (valid_kmers_errs, _): (BoomHashMap2<K, Exts, TagsCountsSumData>, _) = filter::filter_kmers(
            &all_seqs_p,
            &config,
            stranded,
            4,
            true,
        );
        let (valid_kmers_errs2, _): (BoomHashMap2<K, Exts, TagsCountsSumData>, _) = filter::filter_kmers_parallel(
            &all_seqs_p,
            &config,
            stranded,
            4,
            true,
        );

        println!("1: {:?}", valid_kmers_errs);
        println!("2: {:?}", valid_kmers_errs2);

        let (_valid_kmers_errs3, _): (BoomHashMap2<K, Exts, u32>, _) = filter::filter_kmers_parallel(
            &all_seqs_p,
            &config,
            stranded,
            4,
            true,
        );
        let (_valid_kmers_errs4, _): (BoomHashMap2<K, Exts, TagsCountsEMData>, _) = filter::filter_kmers_parallel(
            &all_seqs_p,
            &config,
            stranded,
            4,
            true,
        );

        println!("3: {:?}", _valid_kmers_errs3);
        println!("4: {:?}", _valid_kmers_errs4);

        //let spec = SimpleCompress::new(|d1: u16, d2: &u16| d1 + d2);
        let spec = ScmapCompress::new();
        let graph = compress_kmers_with_hash(stranded, &spec, &valid_kmers_errs, true, true, true);
        println!("graph: {:?}", graph);
        let graph = compress_kmers_with_hash(stranded, &spec, &valid_kmers_errs, true, false, true);
        println!("graph: {:?}", graph);

        let mut graph = graph.finish();
        graph.fix_exts(None);
        graph.fix_edge_mults();
        graph.print();
        /* graph.to_dot("test_out", &|d| format!("{:?}", d));
        graph.to_dot_parallel("test_out_par", &|d  format!("{:?}", d)); */
        
        println!("components i: {:?}", graph.components_i());
        println!("components r: {:?}", graph.components_r());

        for path in graph.max_path_comp(|_| 1., |_| true) {
            println!("path m1: {:?}", graph.sequence_of_path(path.iter()))
        }
        
        let path_iter = graph.iter_max_path_comp(|_| 1., |_| true);

        let mut file = File::create("test_fasta_out.fasta").unwrap();
        graph.path_to_fasta(&mut file, path_iter, false);

        for x in graph.iter_components() {
            println!("component: {:?}", x);
        }

        for path in graph.iter_max_path_comp(|_| 1., |_| true) {
            println!("path seq: {:?}", path);
            println!("path seq: {:?}", graph.sequence_of_path(path.1.iter()));
        }

        graph.to_gfa_with_tags("gfa_out_seq.gfa", |node| format!("{:?}", node.data())).unwrap();
        graph.to_gfa_otags_parallel("gfa_out_par.gfa", Some(&|node: &graph::Node<K, TagsCountsSumData>| format!("{:?}", node.data()))).unwrap();

        //let graph2 = graph.finish();
        //graph2.print();
        //graph.print();
        //println!("components: {:?}", graph.components_r());
        //let max_path = graph2.max_path(|d| *d as f32, |_| true);
        //println!("one graph: {:?}", max_path); 
        //let max_path_c = graph2.max_path_comp(|d| *d as f32, |_| true);
        //println!("all graphs: {:?}", max_path_c); 
        /* for i in 0..graph2.len() {
            println!("node {}: {}", i,  graph2.get_node(i).sequence());
            println!("rc   {}: {}", i,  graph2.get_node(i).sequence().rc());
        } */
        //let u_graph: BaseGraph<K, (Vec<()>, i32)> = uncompressed_graph(&valid_kmers_errs);
        //let u_graph2 = u_graph.finish();
        //println!("uncompressed graph: {:?}", u_graph2);
        //u_graph2.print();

        // Now try to clean the tips.
        //let cleaner = CleanGraph::new(|node| node.len() < K::k() * 2);
        //let nodes_to_censor = cleaner.find_bad_nodes(&graph2);

        //println!("censor: {:?}", nodes_to_censor);
        //let spec = SimpleCompress::new(|d1: u16, d2: &u16| d1 + d2);
        //let fixed = compress_graph(stranded, &spec, graph2, Some(nodes_to_censor));
        //fixed.print();
    }

    #[test]
    fn test_tags() {
        // translation hash map
        let str_vec = vec!["tag1", "tag2", "tag3", "tag4", "tag5", "tag6", "tag7"];
        let mut str_map = BiMap::new();
        for (i, str) in str_vec.into_iter().enumerate() {
            str_map.insert(str.to_string(), i as u8);
        }

        // build Tags from val
        let tag = Tags::new(83);
        assert_eq!(tag.val, 83);
        assert_eq!(tag.to_u8_vec(), vec![0, 1, 4, 6]);
        assert!(format!("{:064b}", tag.val).ends_with("1010011"));
        assert_eq!(tag.to_string_vec(&str_map), vec!["tag1", "tag2", "tag5", "tag7"]);
        
        // build tags from u8 vec
        let vec = Tags::from_u8_vec(vec![0, 1, 4, 6]);
        let vec2 = Tags::from_u8_vec(vec![1, 2, 3, 4, 6]);

        assert_eq!(vec.val, 83);
        assert_eq!(vec2.val, 94);
        assert_eq!(vec.to_u8_vec(), vec![0, 1, 4, 6]);
        assert_eq!(vec2.to_u8_vec(), vec![1, 2, 3, 4, 6]);
        assert_eq!(vec.to_string_vec(&str_map), vec!["tag1", "tag2", "tag5", "tag7"]);
        assert_eq!(vec2.to_string_vec(&str_map), vec!["tag2", "tag3", "tag4", "tag5", "tag7"]);
    }

    #[test]
    #[should_panic]
    #[cfg(not(feature = "sample128"))]
    fn test_tags_overflow() {
        let _tags = Tags::from_u8_vec(vec![5, 64]);
    }

    #[test]
    #[cfg(not(feature = "sample128"))]
    fn test_tags_no_overflow() {
        let _tags = Tags::from_u8_vec(vec![5, 63]);
    }

    #[test]
    #[should_panic]
    #[cfg(feature = "sample128")]
    fn test_tags_overflow() {
        let _tags = Tags::from_u8_vec(vec![5, 128]);
    }

    #[test]
    #[cfg(feature = "sample128")]
    fn test_tags_no_overflow() {
        let _tags = Tags::from_u8_vec(vec![5, 127]);
    }

    #[test]
    fn kmer_experiment() {
        let mut counts = [0; 256];
        let max = (4usize.pow(6) - 1) as u16;
        for i in 0..=max {
            let kmer = crate::kmer::Kmer6::from(VarIntKmer { storage: i, phantom: PhantomData });
            let min_rc = kmer.min_rc();
            counts[filter::bucket(min_rc)] += 1;
        }
        println!("{:?}", counts);
    }

    #[test]
    fn test_speed_calculation() {
        let max_len = 150;
        const WIDTH: i32 = 2;
        let start = time::Instant::now();
        for _i in 0..100000000 {
            let _ = (max_len as f64 / 32.).ceil() as usize;
        }
        let after_simple = start.elapsed();
        let start = time::Instant::now();
        for _i in 0..100000000 {
            let _ = ((max_len * WIDTH) >> 6) + (if (max_len * WIDTH) & 0x3F > 0 { 1 } else { 0 });
        }
        let after_complicated = start.elapsed();
        println!("simple method: {} s\ncomplicated method: {} s", after_simple.as_secs_f32(), after_complicated.as_nanos());

    }


    #[test]
    fn test_filter_kmers() {
        let fastq = [
            (DnaString::from_dna_string("ACGATCGATCGTAGCTACG"), Exts::empty(), 6u8),
            (DnaString::from_dna_string("ATGCAGCTATTGCGAGCTACTTCAG"), Exts::empty(), 6u8),
            (DnaString::from_dna_string("GTCGCAGATCGACTGAGCTAGC"), Exts::empty(), 7u8),
            (DnaString::from_dna_string("GCGATCTACTGACGGATCTATCGGAGCTA"), Exts::empty(), 3u8),
            (DnaString::from_dna_string("GCGATCTAGCGGATCTGCGAGCTATGC"), Exts::empty(), 6u8),
        ];

        let mut reads = Reads::new(Strandedness::Unstranded);

        for (read, exts, data) in fastq {
            reads.add_read(read, exts, data);
        }

        let sample_info = SampleInfo::new(0, 0, 0, 0, Vec::new());
        let config = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info, None, crate::summarizer::StatTest::StudentsTTest);

        let hm: (BoomHashMap2<Kmer6, Exts, TagsSumData>, Vec<_>) = filter_kmers(
            &crate::reads::ReadsPaired::Unpaired { reads }, 
            &config,
            false, 
            1,
            false,
         );

         println!("{:?}", hm);

    }
}

