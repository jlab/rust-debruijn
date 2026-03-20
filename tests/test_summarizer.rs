
use std::mem;

use bimap::BiHashMap;
use debruijn::{BaseQuality, EdgeMult, Exts, Kmer, KmerDataItem, Tags, kmer::Kmer8, reads::ReadData, summarizer::{GroupCountData, GroupFrac, ID, IDData, IDEMData, IDMapEMData, IDMapEMQualityData, IDSumData, IDTagsCountsData, IDTagsCountsPEMData, RelCountData, SampleInfo, Summarizers, SummaryConfig, SummaryData, Tag, TagsCountsData, TagsCountsEMData, TagsCountsPData, TagsCountsPEMData, TagsCountsPEMQualityData, TagsCountsSumData, TagsData, TagsSumData, Translator}};

#[derive(Debug, PartialEq)]
struct SummaryTest {
    print: String,
    print_ol: String,
    print_json: String, 
    tags: Option<Tags>,
    mem: usize,
    sum: Option<usize>,
    ids: Option<Vec<ID>>,
    p_value: Option<f32>,
    fold_change: Option<f32>,
    sample_count: Option<usize>,
    edge_mults: Option<EdgeMult>,
    quality: Option<BaseQuality>,
    mapped_ids: Option<Vec<ID>>, // also include set_ma
    valid: bool,
    summarizer: Summarizers
}

fn test_summarize<'a, SD: SummaryData<DI>, K: Kmer, DI, F>(items: F, config: &'a SummaryConfig, translator: &'a Translator ) 
    -> SummaryTest
where 
    F: Iterator<Item = KmerDataItem<K, DI>>,
{
    let (valid, _, mut data) = SD::summarize(items, config);

    let sum = data.sum();
    let tags = data.tags();
    let sample_count = data.sample_count();

    if let Some(t) = tags {
        assert_eq!(t.to_tag_vec().len(), sample_count.unwrap());
    }

    assert_eq!(data.valid(config), valid);

    let edge_mults = data.edge_mults().cloned();

    data.fix_edge_mults(Exts::new(0));
    if let Some(e) = data.edge_mults() {
        assert_eq!(e.edge_mults(), [0; 8]);
    }

    data.set_edge_mults(Some(EdgeMult::new_from([0, 0, 0, 0, 0, 0, 0, 1])));
    if let Some(e) = data.edge_mults() {
        assert_eq!(e.edge_mults(), [0, 0, 0, 0, 0, 0, 0, 1]);
    }

    assert!(data.join_test(&data));

    data.set_mapped_ids([1, 2, 3].into());
    if let Some(m) = data.mapped_ids() {
        assert_eq!(m, vec![1, 2, 3])
    }

    let mapped_ids = data.mapped_ids().map(|m| m.into());
    let ids = data.ids().map(|ids| ids.to_vec());

    SummaryTest { 
        print: data.print(translator, config, None),
        print_ol: data.print_ol(translator, config, None), 
        print_json: data.print_json(translator, config, None), 
        tags, 
        mem: data.mem(), 
        sum, 
        ids, 
        p_value: data.p_value(config), 
        fold_change: data.fold_change(config), 
        sample_count, 
        edge_mults:  edge_mults.clone(), 
        quality: data.quality(), 
        mapped_ids, 
        valid, 
        summarizer: SD::summarizer() 
    }
}

/// generate summary data to test
fn get_summary_input<DI: ReadData>() -> ([KmerDataItem<Kmer8, DI>; 6], SummaryConfig, Translator) {
    let marker0 = 0b0000000001111;
    let marker1 = 0b1111111110000;
    let sample_kmers = vec![2834, 2343, 12, 1234, 345345, 122, 234, 23455, 231, 2, 3564, 12344, 34555];
    let sample_info = SampleInfo::new(marker0, marker1, sample_kmers);
    let summary_config = SummaryConfig::new(sample_info.clone());
    let mut tag_translator = BiHashMap::new();
    (0..15).for_each(|i| { tag_translator.insert(format!("{i}"), i as u8); } );
    let mut id_translator = BiHashMap::new();
    (0..15).for_each(|i| { id_translator.insert(format!("{i}"), i as ID); } );
    let translator = Translator::new(id_translator, tag_translator);

    let input = [
        KmerDataItem::new(Kmer8::from_u64(12), Exts::new(1), DI::new(0, 0), Some(debruijn::BaseQuality::Medium)),
        KmerDataItem::new(Kmer8::from_u64(12), Exts::new(1), DI::new(1, 1), Some(debruijn::BaseQuality::Medium)),
        KmerDataItem::new(Kmer8::from_u64(12), Exts::new(1), DI::new(2, 2), Some(debruijn::BaseQuality::Medium)),
        KmerDataItem::new(Kmer8::from_u64(12), Exts::new(1), DI::new(3, 3), Some(debruijn::BaseQuality::Medium)),
        KmerDataItem::new(Kmer8::from_u64(12), Exts::new(1), DI::new(7, 7), Some(debruijn::BaseQuality::Medium)),
        KmerDataItem::new(Kmer8::from_u64(12), Exts::new(1), DI::new(8, 8), Some(debruijn::BaseQuality::Medium)),           
    ];

    (input, summary_config, translator)
}

#[test]
fn test_summary_data() {
    // input

    let (input_tags, mut config, translator) = get_summary_input();
    let (input_id_tags, _, _) = get_summary_input();

    let sum = Some(input_tags.len());
    let s_m = 8; // usize
    
    let tags = Some(Tags::from_tag_vec(vec![0, 1, 2, 3, 7, 8]));
    let t_m = tags.as_ref().unwrap().mem(); // Marker

    let c_m_s = 16; // Boxed slice
    let c_m_h = 6 * 4; 
    
    let sample_count = Some(6);
    
    let p_value = Some(0.39023498);
    let p_m = 4; // f32
    
    let fold_change = Some(5.4498405);

    let edge_mults = Some(EdgeMult::new_from([0, 0, 0, 0, 0, 0, 0, 6]));
    let em_m = 8 * 4; // EdgeMult


    let ids = Some(vec![0, 1, 2, 3, 7, 8]);
    let i_m_s = 16; // Boxed slice
    let i_m_h = 6 * mem::size_of::<ID>(); // ids, depend on feature
    let mapped_ids = Some(vec![1, 2, 3]);
    let mi_m_s = 16; // Boxed slice
    let mi_m_h = 3 * mem::size_of::<ID>(); // ids, depend on feature

    let quality = Some(BaseQuality::Medium);
    let q_m = 1; // BaseQuality (enum -> u8)

    let valid = true;

    // u32

    let test_data = test_summarize::<u32, Kmer8, _, _>(input_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "sum: 6".to_string(),
        print_ol: "sum: 6".to_string(),
        print_json: "\"sum\": 6".to_string(),
        tags: None,
        mem: 4,
        sum,
        ids: None,
        p_value: None,
        fold_change: None,
        sample_count: None,
        edge_mults: None,
        quality: None,
        mapped_ids: None,
        valid,
        summarizer: Summarizers::Sum,
    };

    assert_eq!(test_data, compare_data);

    // Vec<Tag>

    let test_data = test_summarize::<Vec<Tag>, Kmer8, _, _>(input_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "samples: ['0', '1', '2', '3', '7', '8']".to_string(),
        print_ol: "samples: ['0', '1', '2', '3', '7', '8']".to_string(),
        print_json: "\"samples\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"]".to_string(),
        tags,
        mem: size_aligned(3*8, 6, 8),
        sum: None,
        ids: None,
        p_value: None,
        fold_change: None,
        sample_count,
        edge_mults: None,
        quality: None,
        mapped_ids: None,
        valid,
        summarizer: Summarizers::VecTags,
    };

    assert_eq!(test_data, compare_data);

    // IDData

    let test_data = test_summarize::<IDData, Kmer8, _, _>(input_id_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "IDs: ['0', '1', '2', '3', '7', '8']".to_string(),
        print_ol: "IDs: ['0', '1', '2', '3', '7', '8']".to_string(),
        print_json: "\"ids\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"]".to_string(),
        tags: None,
        mem: size_aligned(i_m_s, i_m_h, 8),
        sum: None,
        ids:  ids.clone(),
        p_value: None,
        fold_change: None,
        sample_count: None,
        edge_mults: None,
        quality: None,
        mapped_ids: None,
        valid,
        summarizer: Summarizers::ID,
    };

    assert_eq!(test_data, compare_data);

    // IDSumData

    let test_data = test_summarize::<IDSumData, Kmer8, _, _>(input_id_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "IDs: ['0', '1', '2', '3', '7', '8'], sum: 6".to_string(),
        print_ol: "IDs: ['0', '1', '2', '3', '7', '8'], sum: 6".to_string(),
        print_json: "\"ids\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"], \"sum\": 6".to_string(),
        tags: None,
        mem: size_aligned(i_m_s + s_m, i_m_h, 8),
        sum,
        ids:  ids.clone(),
        p_value: None,
        fold_change: None,
        sample_count: None,
        edge_mults: None,
        quality: None,
        mapped_ids: None,
        valid,
        summarizer: Summarizers::IDSum,
    };

    assert_eq!(test_data, compare_data);
    
    // TagsData

    let test_data = test_summarize::<TagsData, Kmer8, _, _>(input_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "samples:\n0\n1\n2\n3\n7\n8\n".to_string(),
        print_ol: "samples: ['0', '1', '2', '3', '7', '8']".to_string(),
        print_json: "\"samples\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"]".to_string(),
        tags,
        mem: t_m,
        sum: None,
        ids: None,
        p_value: None,
        fold_change: None,
        sample_count,
        edge_mults: None,
        quality: None,
        mapped_ids: None,
        valid,
        summarizer: Summarizers::Tags,
    };

    assert_eq!(test_data, compare_data);

    // TagsSumData

    let test_data = test_summarize::<TagsSumData, Kmer8, _, _>(input_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "samples:\n0\n1\n2\n3\n7\n8\nsum: 6".to_string(),
        print_ol: "samples: ['0', '1', '2', '3', '7', '8'], sum: 6".to_string(),
        print_json: "\"samples\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"], \"sum\": 6".to_string(),
        tags,
        mem: size_aligned(t_m + s_m, 0, t_m),
        sum,
        ids: None,
        p_value: None,
        fold_change: None,
        sample_count,
        edge_mults: None,
        quality: None,
        mapped_ids: None,
        valid,
        summarizer: Summarizers::TagsSum,
    };

    assert_eq!(test_data, compare_data);

    // TagsCountsSumData

    let test_data = test_summarize::<TagsCountsSumData, Kmer8, _, _>(input_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(),
        print_ol: "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(),
        print_json: "\"sum\": 6, \"samples\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"], \"counts\": [1, 1, 1, 1, 1, 1], \"p_value\": 0.39023498, \"fold_change\": 5.4498405".to_string(),
        tags,
        mem: size_aligned(t_m + c_m_s + s_m, c_m_h, t_m),
        sum,
        ids: None,
        p_value,
        fold_change,
        sample_count,
        edge_mults: None,
        quality: None,
        mapped_ids: None,
        valid,
        summarizer: Summarizers::TagsCountsSum,
    };

    assert_eq!(test_data, compare_data);

    // TagsCountsData

    let test_data = test_summarize::<TagsCountsData, Kmer8, _, _>(input_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(),
        print_ol: "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(),
        print_json: "\"sum\": 6, \"samples\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"], \"counts\": [1, 1, 1, 1, 1, 1], \"p_value\": 0.39023498, \"fold_change\": 5.4498405".to_string(),
        tags,
        mem: size_aligned(t_m + c_m_s, c_m_h, t_m),
        sum,
        ids: None,
        p_value,
        fold_change,
        sample_count,
        edge_mults: None,
        quality: None,
        mapped_ids: None,
        valid,
        summarizer: Summarizers::TagsCounts,
    };

    assert_eq!(test_data, compare_data);

    // TagsCountsPData

    let test_data = test_summarize::<TagsCountsPData, Kmer8, _, _>(input_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(),
        print_ol: "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(),
        print_json: "\"sum\": 6, \"samples\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"], \"counts\": [1, 1, 1, 1, 1, 1], \"p_value\": 0.39023498, \"fold_change\": 5.4498405".to_string(),
        tags,
        mem: size_aligned(t_m + c_m_s + p_m, c_m_h, t_m),
        sum,
        ids: None,
        p_value,
        fold_change,
        sample_count,
        edge_mults: None,
        quality: None,
        mapped_ids: None,
        valid,
        summarizer: Summarizers::TagsCountsP,
    };

    assert_eq!(test_data, compare_data);

    // TagsCountsEMData

    let test_data = test_summarize::<TagsCountsEMData, Kmer8, _, _>(input_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge coverage: \nA: 1 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        print_ol: "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge coverage: A: 1, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string(),
        print_json: "\"sum\": 6, \"samples\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"], \"counts\": [1, 1, 1, 1, 1, 1], \"p_value\": 0.39023498, \"fold_change\": 5.4498405".to_string(),
        tags,
        mem: size_aligned(t_m + c_m_s + em_m, c_m_h, t_m),
        sum,
        ids: None,
        p_value,
        fold_change,
        sample_count,
        edge_mults:  edge_mults.clone(),
        quality: None,
        mapped_ids: None,
        valid,
        summarizer: Summarizers::TagsCountsEM,
    };

    assert_eq!(test_data, compare_data);

    // TagsCountsPEMData

    let test_data = test_summarize::<TagsCountsPEMData, Kmer8, _, _>(input_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge coverage: \nA: 1 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        print_ol: "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge coverage: A: 1, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string(),
        print_json: "\"sum\": 6, \"samples\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"], \"counts\": [1, 1, 1, 1, 1, 1], \"p_value\": 0.39023498, \"fold_change\": 5.4498405".to_string(),
        tags,
        mem: size_aligned(t_m + c_m_s + p_m + em_m, c_m_h, t_m),
        sum,
        ids: None,
        p_value,
        fold_change,
        sample_count,
        edge_mults:  edge_mults.clone(),
        quality: None,
        mapped_ids: None,
        valid,
        summarizer: Summarizers::TagsCountsPEM,
    };

    assert_eq!(test_data, compare_data);

    // TagsCountsPEMQualityData

    let test_data = test_summarize::<TagsCountsPEMQualityData, Kmer8, _, _>(input_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, quality: medium, edge coverage: \nA: 1 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        print_ol: "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, quality: medium, edge coverage: A: 1, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string(),
        print_json: "\"sum\": 6, \"samples\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"], \"counts\": [1, 1, 1, 1, 1, 1], \"p_value\": 0.39023498, \"fold_change\": 5.4498405, \"quality\": 2".to_string(),
        tags,
        mem: size_aligned(t_m + c_m_s + p_m + em_m + q_m, c_m_h, t_m),
        sum,
        ids: None,
        p_value,
        fold_change,
        sample_count,
        edge_mults:  edge_mults.clone(),
        quality,
        mapped_ids: None,
        valid,
        summarizer: Summarizers::TagsCountsPEMQuality,
    };

    assert_eq!(test_data, compare_data);

    // IDTagsCountsData

    let test_data = test_summarize::<IDTagsCountsData, Kmer8, _, _>(input_id_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "IDs: ['0', '1', '2', '3', '7', '8'], samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(), 
        print_ol: "IDs: ['0', '1', '2', '3', '7', '8'], samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(),
        print_json: "\"ids\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"], \"sum\": 6, \"samples\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"], \"counts\": [1, 1, 1, 1, 1, 1], \"p_value\": 0.39023498, \"fold_change\": 5.4498405".to_string(),
        tags,
        mem: size_aligned(t_m + c_m_s + i_m_s, c_m_h + i_m_h, t_m),
        sum,
        ids:  ids.clone(),
        p_value,
        fold_change,
        sample_count,
        edge_mults: None,
        quality: None,
        mapped_ids: None,
        valid,
        summarizer: Summarizers::IDTagsCounts,
    };

    assert_eq!(test_data, compare_data);

    // IDTagsCountsPEMData

    let test_data = test_summarize::<IDTagsCountsPEMData, Kmer8, _, _>(input_id_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "IDs: ['0', '1', '2', '3', '7', '8'], samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge coverage: \nA: 1 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        print_ol: "IDs: ['0', '1', '2', '3', '7', '8'], samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge coverage: A: 1, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string(),
        print_json: "\"ids\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"], \"sum\": 6, \"samples\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"], \"counts\": [1, 1, 1, 1, 1, 1], \"p_value\": 0.39023498, \"fold_change\": 5.4498405".to_string(),
        tags,
        mem: size_aligned(t_m + c_m_s + i_m_s + p_m + em_m, c_m_h + i_m_h, t_m),
        sum,
        ids:  ids.clone(),
        p_value,
        fold_change,
        sample_count,
        edge_mults:  edge_mults.clone(),
        quality: None,
        mapped_ids: None,
        valid,
        summarizer: Summarizers::IDTagsCountsPEM,
    };

    assert_eq!(test_data, compare_data);

    // IDEMData

    let test_data = test_summarize::<IDEMData, Kmer8, _, _>(input_id_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "IDs: ['0', '1', '2', '3', '7', '8'], edge coverage: \nA: 1 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        print_ol: "IDs: ['0', '1', '2', '3', '7', '8'], edge coverage: A: 1, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string(),
        print_json: "\"ids\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"]".to_string(),
        tags: None,
        mem: size_aligned(i_m_s + em_m, i_m_h, 8),
        sum: None,
        ids:  ids.clone(),
        p_value: None,
        fold_change: None,
        sample_count: None,
        edge_mults:  edge_mults.clone(),
        quality: None,
        mapped_ids: None,
        valid,
        summarizer: Summarizers::IDEM,
    };

    assert_eq!(test_data, compare_data);

    // IDMapEMData

    let test_data = test_summarize::<IDMapEMData, Kmer8, _, _>(input_id_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "IDs: ['0', '1', '2', '3', '7', '8'], mapped IDs: ['1', '2', '3'], edge coverage: A: 1 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        print_ol: "IDs: ['0', '1', '2', '3', '7', '8'], mapped IDs: ['1', '2', '3'], edge coverage: A: 1, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string(),
        print_json: "\"ids\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"], \"mapped_ids\": [\"1\", \"2\", \"3\"], \"has_mapped_ids\": 1".to_string(),
        tags: None,
        mem: size_aligned(i_m_s + mi_m_s + em_m, i_m_h + mi_m_h, 8), // mapped ids empty
        sum: None,
        ids:  ids.clone(),
        p_value: None,
        fold_change: None,
        sample_count: None,
        edge_mults: edge_mults.clone(),
        quality: None,
        mapped_ids: mapped_ids.clone(),
        valid,
        summarizer: Summarizers::IDMapEM,
    };

    assert_eq!(test_data, compare_data);

    // IDMapEMQualityData

    let test_data = test_summarize::<IDMapEMQualityData, Kmer8, _, _>(input_id_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "IDs: ['0', '1', '2', '3', '7', '8'], mapped IDs: ['1', '2', '3'], quality: medium, edge coverage: A: 1 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        print_ol: "IDs: ['0', '1', '2', '3', '7', '8'], mapped IDs: ['1', '2', '3'], quality: medium, edge coverage: A: 1, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string(),
        print_json: "\"ids\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"], \"mapped_ids\": [\"1\", \"2\", \"3\"], \"has_mapped_ids\": 1, \"quality\": 2".to_string(),
        tags: None,
        mem: size_aligned(i_m_s + mi_m_s + em_m + q_m, i_m_h + mi_m_h, 8), // mapped ids empty
        sum: None,
        ids:  ids.clone(),
        p_value: None,
        fold_change: None,
        sample_count: None,
        edge_mults:  edge_mults.clone(),
        quality,
        mapped_ids:  mapped_ids.clone(),
        valid,
        summarizer: Summarizers::IDMapEMQuality,
    };

    assert_eq!(test_data, compare_data);

    // GroupCountData

    let test_data = test_summarize::<GroupCountData, Kmer8, _, _>(input_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "count 1: 4\ncount 2: 2".to_string(), 
        print_ol: "count 1: 4, count 2: 2".to_string(),
        print_json: "\"count1\": 4, \"count2\": 2".to_string(),
        tags: None,
        mem: 2*4,
        sum,
        ids: None,
        p_value: None,
        fold_change: None,
        sample_count: None,
        edge_mults: None,
        quality: None,
        mapped_ids: None,
        valid,
        summarizer: Summarizers::GroupCount,
    };

    assert_eq!(test_data, compare_data);

    // RelCountData

    let test_data = test_summarize::<RelCountData, Kmer8, _, _>(input_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "relative amount group 1: 66\ncount both: 6".to_string(), 
        print_ol: "relative amount group 1: 66, count both: 6".to_string(),
        print_json: "\"rel_count_1\": 66, \"sum\": 6".to_string(),
        tags: None,
        mem: 2*4,
        sum,
        ids: None,
        p_value: None,
        fold_change: None,
        sample_count: None,
        edge_mults: None,
        quality: None,
        mapped_ids: None,
        valid,
        summarizer: Summarizers::RelCount,
    };

    assert_eq!(test_data, compare_data);



    // test with different group frax settings

    // One
    config.set_group_frac(GroupFrac::One, 0.33);

    let test_data = test_summarize::<TagsData, Kmer8, _, _>(input_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "samples:\n0\n1\n2\n3\n7\n8\n".to_string(),
        print_ol: "samples: ['0', '1', '2', '3', '7', '8']".to_string(),
        print_json: "\"samples\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"]".to_string(),
        tags,
        mem: t_m,
        sum: None,
        ids: None,
        p_value: None,
        fold_change: None,
        sample_count,
        edge_mults: None,
        quality: None,
        mapped_ids: None,
        valid,
        summarizer: Summarizers::Tags,
    };

    assert_eq!(test_data, compare_data);

    // Both
    config.set_group_frac(GroupFrac::Both, 0.33);

    let test_data = test_summarize::<TagsData, Kmer8, _, _>(input_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "samples:\n0\n1\n2\n3\n7\n8\n".to_string(),
        print_ol: "samples: ['0', '1', '2', '3', '7', '8']".to_string(),
        print_json: "\"samples\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"]".to_string(),
        tags,
        mem: t_m,
        sum: None,
        ids: None,
        p_value: None,
        fold_change: None,
        sample_count,
        edge_mults: None,
        quality: None,
        mapped_ids: None,
        valid: false,
        summarizer: Summarizers::Tags,
    };

    assert_eq!(test_data, compare_data);

    // Both and with max p
    config.set_max_p(Some(0.05));

    let test_data = test_summarize::<TagsData, Kmer8, _, _>(input_tags.into_iter(), &config, &translator);
    let compare_data = SummaryTest {
        print: "samples:\n0\n1\n2\n3\n7\n8\n".to_string(),
        print_ol: "samples: ['0', '1', '2', '3', '7', '8']".to_string(),
        print_json: "\"samples\": [\"0\", \"1\", \"2\", \"3\", \"7\", \"8\"]".to_string(),
        tags,
        mem: t_m,
        sum: None,
        ids: None,
        p_value: None,
        fold_change: None,
        sample_count,
        edge_mults: None,
        quality: None,
        mapped_ids: None,
        valid: false,
        summarizer: Summarizers::Tags,
    };

    assert_eq!(test_data, compare_data);
}

/// add the alignment buffer to a structure
/// - `size_heap`: contents of boxes/vectors -> are stored separately and do not go into alignment calculation
fn size_aligned(size_stack: usize, size_heap: usize, align: usize) -> usize {
    let empty = size_stack % align;
    let buffer = match empty {
        0 => 0,
        _ => align - empty
    };

    buffer + size_heap + size_stack
}
