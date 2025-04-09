use std::mem;

use bimap::BiHashMap;
use debruijn::{kmer::Kmer8, summarizer::{self, GroupCountData, GroupFrac, RelCountData, SampleInfo, SummaryConfig, SummaryData, TagsCountsData, TagsCountsEMData, TagsCountsPData, TagsCountsPEMData, TagsCountsSumData, TagsSumData, M}, EdgeMult, Exts, Kmer, Tags};

fn test_summarize<SD: SummaryData<u8>, F, K: Kmer>(items: F, config: &SummaryConfig, tag_translator: &bimap::BiHashMap<String, u8> ) 
    -> (Option<usize>, Option<(Tags, u32)>, usize, Option<f32>, Option<f32>, Option<usize>, Option<EdgeMult>, bool, String, String)
where 
    F: Iterator<Item = (K, Exts, u8)>,
{
    let (valid, _, data) = SD::summarize(items, config);

    let count = data.count();
    let score = data.score();

    match count {
        Some(c) => assert_eq!(c, score as usize),
        None => assert_eq!(1., score)
    }
    
    let tags_sum = data.tags_sum();
    let sample_count = data.sample_count();

    if let Some(ts) = tags_sum {
        assert_eq!(ts.1, count.unwrap() as u32);
        assert_eq!(ts.0.to_u8_vec().len(), sample_count.unwrap());
    }

    assert_eq!(data.valid(config), valid);

    let em = data.edge_mults().cloned();

    (
        count,
        tags_sum,
        data.mem(),
        data.p_value(config),
        data.fold_change(config),
        data.sample_count(),
        em,
        valid,
        data.print(tag_translator, config),
        data.print_ol(tag_translator, config),
    )
}

#[test]
fn test_summary_data() {
    // input

    let marker0 = 0b0000000001111;
    let marker1 = 0b1111111110000;
    let count0 = 4;
    let count1 = 9;
    let sample_kmers = vec![2834, 2343, 12, 1234, 345345, 122, 234, 23455, 231, 2, 3564, 12344, 34555];
    let sample_info = SampleInfo::new(marker0, marker1, count0, count1, sample_kmers);
    let summary_config = SummaryConfig::new(1, None, GroupFrac::None, 0.33, sample_info.clone(), None, summarizer::StatTest::WelchsTTest);
    let mut tag_translator = BiHashMap::new();
    (0..15).for_each(|i| { tag_translator.insert(format!("{i}"), i as u8); } );

    let input = [
        (Kmer8::from_u64(12), Exts::new(1), 0u8),
        (Kmer8::from_u64(12), Exts::new(1), 1u8),
        (Kmer8::from_u64(12), Exts::new(1), 2u8),
        (Kmer8::from_u64(12), Exts::new(1), 3u8),
        (Kmer8::from_u64(12), Exts::new(1), 7u8),
        (Kmer8::from_u64(12), Exts::new(1), 8u8),           
    ];

    println!("kmer: {:?}", Kmer8::from_u64(12));

    let size_tags = mem::size_of::<M>();

    let count = Some(input.len());
    let tags_sum = Some((Tags::from_u8_vec(vec![0, 1, 2, 3, 7, 8]), input.len() as u32));
    let sample_count = Some(6);
    let p_value = Some(0.39023498);
    let fold_change = Some(5.4498405);
    let edge_mults = Some(EdgeMult::new_from([0, 0, 0, 0, 0, 0, 0, 6]));

    // test summarize: (Some(count), Some((tags, sum)), memory, Some(p_value), Some(fold_change), Some(sample_count), Some(edge_mults), valid)

    let data = test_summarize::<u32, _, _>(input.into_iter(), &summary_config, &tag_translator);
    assert_eq!(data, (count, None, 4, None, None, None, None, true, "count: 6".to_string(), "count: 6".to_string()));

    let data = test_summarize::<Vec<u8>, _, _>(input.into_iter(), &summary_config, &tag_translator);
    assert_eq!(data, (None, None, 30, None, None, sample_count, None, true, "samples: [0, 1, 2, 3, 7, 8]".to_string(), "samples: [0, 1, 2, 3, 7, 8]".to_string()));

    let data = test_summarize::<TagsSumData, _, _>(input.into_iter(), &summary_config, &tag_translator);
    assert_eq!(data, (count, tags_sum, size_tags * 2, None, None, sample_count, None, true, "samples:\n0\n1\n2\n3\n7\n8\nsum: 6".to_string(), "samples: ['0', '1', '2', '3', '7', '8'], sum: 6".to_string())); // mem: M + 4 + alignment buffer

    let data = test_summarize::<TagsCountsSumData, _, _>(input.into_iter(), &summary_config, &tag_translator);
    assert_eq!(data, (count, tags_sum, size_tags * 2 + 16 + 6*4, p_value, fold_change, sample_count, None, true, 
        "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(), 
        "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string()
    )); // mem: M + 2*8 + 4+ab + 6*4

    let data = test_summarize::<TagsCountsData, _, _>(input.into_iter(), &summary_config, &tag_translator);
    assert_eq!(data, (count, tags_sum, size_tags + 16 + 6*4, p_value, fold_change, sample_count, None, true,
        "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(), 
        "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string()
    )); 

    let data = test_summarize::<TagsCountsPData, _, _>(input.into_iter(), &summary_config, &tag_translator);
    assert_eq!(data, (count, tags_sum, size_tags * 2 + 16 + 6*4, p_value, fold_change, sample_count, None, true,
        "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(), 
        "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string()
    )); 

    let data = test_summarize::<TagsCountsEMData, _, _>(input.into_iter(), &summary_config, &tag_translator);
    assert_eq!(data, (count, tags_sum, size_tags + 16 + 6*4 + 4*8, p_value, fold_change, sample_count, edge_mults.clone(), true, 
        "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge multiplicities: \nA: 6 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge multiplicities: A: 6, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string()
    )); 

    let data = test_summarize::<TagsCountsPEMData, _, _>(input.into_iter(), &summary_config, &tag_translator);
    assert_eq!(data, (count, tags_sum, size_tags * 2 + 16 + 6*4 + 4*8, p_value, fold_change, sample_count, edge_mults.clone(), true,  
        "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge multiplicities: \nA: 6 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge multiplicities: A: 6, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string()
    )); 

    let data = test_summarize::<GroupCountData, _, _>(input.into_iter(), &summary_config, &tag_translator);
    assert_eq!(data, (count, None, 4 + 4, None, None, None, None, true, "count 1: 4\ncount 2: 2".to_string(), "count 1: 4, count 2: 2".to_string()));

    let data = test_summarize::<RelCountData, _, _>(input.into_iter(), &summary_config, &tag_translator);
    assert_eq!(data, (count, None, 4 + 4, None, None, None, None, true, "relative amount group 1: 66\ncount both: 6".to_string(), "relative amount group 1: 66, count both: 6".to_string()));

    let summary_config = SummaryConfig::new(1, None, GroupFrac::One, 0.33, sample_info.clone(), None, summarizer::StatTest::WelchsTTest);
    let data = test_summarize::<TagsCountsPEMData, _, _>(input.into_iter(), &summary_config, &tag_translator);
    assert_eq!(data, (count, tags_sum, size_tags * 2 + 16 + 6*4 + 4*8, p_value, fold_change, sample_count, edge_mults.clone(), true,  
        "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge multiplicities: \nA: 6 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge multiplicities: A: 6, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string()
    )); 

    let summary_config = SummaryConfig::new(1, None, GroupFrac::Both, 0.33, sample_info.clone(), None, summarizer::StatTest::WelchsTTest);
    let data = test_summarize::<TagsCountsPEMData, _, _>(input.into_iter(), &summary_config, &tag_translator);
    assert_eq!(data, (count, tags_sum, size_tags * 2 + 16 + 6*4 + 4*8, p_value, fold_change, sample_count, edge_mults.clone(), false,  
        "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge multiplicities: \nA: 6 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge multiplicities: A: 6, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string()
    )); 

    let summary_config = SummaryConfig::new(1, None, GroupFrac::One, 0.33, sample_info.clone(), Some(0.05), summarizer::StatTest::WelchsTTest);
    let data = test_summarize::<TagsCountsPEMData, _, _>(input.into_iter(), &summary_config, &tag_translator);
    assert_eq!(data, (count, tags_sum, size_tags * 2 + 16 + 6*4 + 4*8, p_value, fold_change, sample_count, edge_mults.clone(), false,  
        "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge multiplicities: \nA: 6 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge multiplicities: A: 6, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string()
    )); 
}