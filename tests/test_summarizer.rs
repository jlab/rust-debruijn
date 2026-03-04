
use bimap::BiHashMap;
use debruijn::{EdgeMult, Exts, Kmer, KmerDataItem, Tags, kmer::Kmer8, summarizer::{GroupCountData, ID, IDData, IDEMData, IDMapEMData, IDMapEMQualityData, IDSumData, IDTag, IDTagsCountsData, IDTagsCountsPEMData, RelCountData, SampleInfo, Summarizers, SummaryConfig, SummaryData, TagsCountsData, TagsCountsEMData, TagsCountsPData, TagsCountsPEMData, TagsCountsSumData, TagsData, TagsSumData, Translator}};

fn test_summarize<'a, SD: SummaryData<DI>, F, K: Kmer, DI>(items: F, config: &'a SummaryConfig, translator: &'a Translator ) 
    -> (Option<usize>, Option<Tags>, usize, Option<f32>, Option<f32>, Option<usize>, Option<Vec<ID>>, Option<EdgeMult>, bool, String, String, Summarizers,)
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

    let em = data.edge_mults().cloned();

    data.fix_edge_mults(Exts::new(0));
    if let Some(e) = data.edge_mults() {
        assert_eq!(e.edge_mults(), [0; 8]);
    }

    data.set_edge_mults(Some(EdgeMult::new_from([0, 0, 0, 0, 0, 0, 0, 1])));

    assert!(data.join_test(&data));

    let ids = data.ids().map(|ids| ids.to_vec());
    
    (
        sum,
        tags,
        data.mem(),
        data.p_value(config),
        data.fold_change(config),
        data.sample_count(),
        ids,
        em,
        valid,
        data.print(translator, config, None),
        data.print_ol(translator, config, None),
        SD::summarizer()
    )
}

// memory requirements of a SummaryData instance depends not only on the struct 
// but also on the feature

// neither features
#[cfg(all(not(feature= "id4b"), not(feature = "sample128")))]
//                       sum id         id-sum             vec       tags tags-sum  
const MEM: [usize; 15] = [4, 8+8 + 6*2, 8+8 + 6*2 + 4 + 4, 8+8+8 + 6, 8, 8 + 4 + 4, 
//  tags-counts-sum       tags-counts    tags-counts-p      tags-counts-em    
    8 + 8+8 + 6*4 + 4 + 4, 8 + 8+8 + 6*4, 8 + 8+8 + 6*4 + 4 + 4, 8 + 8+8 + 6*4 + 8*4, 
//  tags-counts-p-em             id-tags-counts             id-tags-counts-p-em
    8 + 8+8 + 6*4 + 4 + 8*4 + 4, 8+8 + 6*2 + 8 + 8+8 + 6*4, 8+8 + 6*2 + 8+8 + 6*4 + 8 + 4 + 8*4 + 4,
//  group-counts
    4 + 4, 4 + 4
];

// only sample128 feature
#[cfg(all(not(feature= "id4b"), feature = "sample128"))]
//                       sum id         id-sum             vec        tags tags-sum  
const MEM: [usize; 15] = [4, 8+8 + 6*2, 8+8 + 6*2 + 4 + 4, 8+8+8 + 6, 16, 16 + 4 + 12, 
//  tags-counts-sum          tags-counts     tags-counts-p           tags-counts-em    
    16 + 8+8 + 6*4 + 4 + 12, 16 + 8+8 + 6*4, 16 + 8+8 + 6*4 + 4 + 12, 16 + 8+8 + 6*4 + 8*4, 
//  tags-counts-p-em               id-tags-counts              id-tags-counts-p-em
    16 + 8+8 + 6*4 + 4 + 8*4 + 12, 8+8 + 6*2 + 16 + 8+8 + 6*4, 8+8 + 6*2 + 8+8 + 6*4 + 16 + 4 + 8*4 + 12,
//  group-counts
    4 + 4, 4 + 4
];

// only id4b feature
#[cfg(all(feature= "id4b", not(feature = "sample128")))]
//                       sum id         id-sum             vec       tags tags-sum  
const MEM: [usize; 15] = [4, 8+8 + 6*4, 8+8 + 6*4 + 4 + 4, 8+8+8 + 6, 8, 8 + 4 + 4, 
//  tags-counts-sum       tags-counts    tags-counts-p      tags-counts-em    
    8 + 8+8 + 6*4 + 4 + 4, 8 + 8+8 + 6*4, 8 + 8+8 + 6*4 + 4 + 4, 8 + 8+8 + 6*4 + 8*4, 
//  tags-counts-p-em             id-tags-counts             id-tags-counts-p-em
    8 + 8+8 + 6*4 + 4 + 8*4 + 4, 8+8 + 6*4 + 8 + 8+8 + 6*4, 8+8 + 6*4 + 8+8 + 6*4 + 8 + 4 + 8*4 + 4,
//  group-counts
    4 + 4, 4 + 4
];

// both features
#[cfg(all(feature= "id4b", feature = "sample128"))]
//                       sum id         id-sum             vec       tags tags-sum  
const MEM: [usize; 15] = [4, 8+8 + 6*4, 8+8 + 6*4 + 4 + 4, 8+8+8 + 6, 16, 16 + 4 + 12, 
//  tags-counts-sum       tags-counts    tags-counts-p      tags-counts-em    
    16 + 8+8 + 6*4 + 4 + 12, 16 + 8+8 + 6*4, 16 + 8+8 + 6*4 + 4 + 12, 16 + 8+8 + 6*4 + 8*4, 
//  tags-counts-p-em             id-tags-counts             id-tags-counts-p-em
    16 + 8+8 + 6*4 + 4 + 8*4 + 12, 8+8 + 6*4 + 16 + 8+8 + 6*4, 8+8 + 6*4 + 8+8 + 6*4 + 16 + 4 + 8*4 + 12,
//  group-counts
    4 + 4, 4 + 4
];


#[test]
fn test_summary_data() {
    // input

    let marker0 = 0b0000000001111;
    let marker1 = 0b1111111110000;
    let sample_kmers = vec![2834, 2343, 12, 1234, 345345, 122, 234, 23455, 231, 2, 3564, 12344, 34555];
    let sample_info = SampleInfo::new(marker0, marker1, sample_kmers);
    let summary_config = SummaryConfig::new(sample_info.clone());
    println!("summary config: {:?}", summary_config);
    let mut tag_translator = BiHashMap::new();
    (0..15).for_each(|i| { tag_translator.insert(format!("{i}"), i as u8); } );
    let mut id_translator = BiHashMap::new();
    (0..15).for_each(|i| { id_translator.insert(format!("{i}"), i as ID); } );
    let translator = Translator::new(id_translator, tag_translator);

    let input_tags = [
        KmerDataItem::new(Kmer8::from_u64(12), Exts::new(1), 0u8, Some(debruijn::BaseQuality::Medium)),
        KmerDataItem::new(Kmer8::from_u64(12), Exts::new(1), 1u8, Some(debruijn::BaseQuality::Medium)),
        KmerDataItem::new(Kmer8::from_u64(12), Exts::new(1), 2u8, Some(debruijn::BaseQuality::Medium)),
        KmerDataItem::new(Kmer8::from_u64(12), Exts::new(1), 3u8, Some(debruijn::BaseQuality::Medium)),
        KmerDataItem::new(Kmer8::from_u64(12), Exts::new(1), 7u8, Some(debruijn::BaseQuality::Medium)),
        KmerDataItem::new(Kmer8::from_u64(12), Exts::new(1), 8u8, Some(debruijn::BaseQuality::Medium)),           
    ];

    let input_ids = [
        KmerDataItem::new(Kmer8::from_u64(12), Exts::new(1), 0 as ID, Some(debruijn::BaseQuality::Medium)),
        KmerDataItem::new(Kmer8::from_u64(12), Exts::new(1), 1 as ID, Some(debruijn::BaseQuality::Medium)),
        KmerDataItem::new(Kmer8::from_u64(12), Exts::new(1), 2 as ID, Some(debruijn::BaseQuality::Medium)),
        KmerDataItem::new(Kmer8::from_u64(12), Exts::new(1), 3 as ID, Some(debruijn::BaseQuality::Medium)),
        KmerDataItem::new(Kmer8::from_u64(12), Exts::new(1), 7 as ID, Some(debruijn::BaseQuality::Medium)),
        KmerDataItem::new(Kmer8::from_u64(12), Exts::new(1), 8 as ID, Some(debruijn::BaseQuality::Medium)),   
    ];

    let input_id_tags = input_tags
        .iter()
        .zip(&input_ids)
        .map(|(tags_item, ids_item)| KmerDataItem::new(tags_item.kmer, tags_item.exts, IDTag::new(ids_item.data, tags_item.data), tags_item.quality))
        .collect::<Vec<_>>();

    println!("kmer: {:?}", Kmer8::from_u64(12));

    let count = Some(input_tags.len());
    let tags = Some(Tags::from_tag_vec(vec![0, 1, 2, 3, 7, 8]));
    let sample_count = Some(6);
    let p_value = Some(0.39023498);
    let fold_change = Some(5.4498405);
    let edge_mults = Some(EdgeMult::new_from([0, 0, 0, 0, 0, 0, 0, 6]));

    // test summarize: (Some(count), Some((tags, sum)), memory, Some(p_value), Some(fold_change), Some(sample_count), Some(edge_mults), valid, "print", "print ol")

    let data = test_summarize::<u32, _, _, _>(input_tags.into_iter(), &summary_config, &translator);
    assert_eq!(data, (count, None, MEM[0], None, None, None, None, None, true, 
        "sum: 6".to_string(), "sum: 6".to_string(), Summarizers::Sum));

    let data = test_summarize::<IDData, _, _, _>(input_ids.into_iter(), &summary_config, &translator);
    assert_eq!(data, (None, None, MEM[1], None, None, None, Some(vec![0, 1, 2, 3, 7, 8]), None, true, 
        "IDs: ['0', '1', '2', '3', '7', '8']".to_string(), "IDs: ['0', '1', '2', '3', '7', '8']".to_string(), Summarizers::ID));

    let data = test_summarize::<IDSumData, _, _, _>(input_ids.into_iter(), &summary_config, &translator);
    assert_eq!(data, (count, None, MEM[2], None, None, None, Some(vec![0, 1, 2, 3, 7, 8]), None, true, 
        "IDs: ['0', '1', '2', '3', '7', '8'], sum: 6".to_string(), "IDs: ['0', '1', '2', '3', '7', '8'], sum: 6".to_string(), Summarizers::IDSum));

    let data = test_summarize::<Vec<u8>, _, _, _>(input_tags.into_iter(), &summary_config, &translator);
    assert_eq!(data, (None, None, MEM[3], None, None, sample_count, None, None, true, 
        "samples: ['0', '1', '2', '3', '7', '8']".to_string(), "samples: ['0', '1', '2', '3', '7', '8']".to_string(), Summarizers::VecTags));

    let data = test_summarize::<TagsData, _, _, _>(input_tags.into_iter(), &summary_config, &translator);
    assert_eq!(data, (None, tags, MEM[4], None, None, sample_count, None, None, true, 
        "samples:\n0\n1\n2\n3\n7\n8\n".to_string(), "samples: ['0', '1', '2', '3', '7', '8']".to_string(), Summarizers::Tags));

    let data = test_summarize::<TagsSumData, _, _, _>(input_tags.into_iter(), &summary_config, &translator);
    assert_eq!(data, (count, tags, MEM[5], None, None, sample_count, None, None, true, 
        "samples:\n0\n1\n2\n3\n7\n8\nsum: 6".to_string(), "samples: ['0', '1', '2', '3', '7', '8'], sum: 6".to_string(), Summarizers::TagsSum)); // mem: M + 4 + alignment buffer

    let data = test_summarize::<TagsCountsSumData, _, _, _>(input_tags.into_iter(), &summary_config, &translator);
    assert_eq!(data, (count, tags, MEM[6], p_value, fold_change, sample_count, None, None, true, 
        "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(), 
        "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(),
        Summarizers::TagsCountsSum
    )); // mem: M + 2*8 + 4+ab + 6*4

    let data = test_summarize::<TagsCountsData, _, _, _>(input_tags.into_iter(), &summary_config, &translator);
    assert_eq!(data, (count, tags, MEM[7], p_value, fold_change, sample_count, None, None, true,
        "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(), 
        "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(),
        Summarizers::TagsCounts
    )); 

    let data = test_summarize::<TagsCountsPData, _, _, _>(input_tags.into_iter(), &summary_config, &translator);
    assert_eq!(data, (count, tags, MEM[8], p_value, fold_change, sample_count, None, None, true,
        "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(), 
        "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(),
        Summarizers::TagsCountsP
    )); 

    let data = test_summarize::<TagsCountsEMData, _, _, _>(input_tags.into_iter(), &summary_config, &translator);
    assert_eq!(data, (count, tags, MEM[9], p_value, fold_change, sample_count, None, edge_mults.clone(), true, 
        "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge coverage: \nA: 1 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge coverage: A: 1, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string(),
        Summarizers::TagsCountsEM
    )); 

    let data = test_summarize::<TagsCountsPEMData, _, _, _>(input_tags.into_iter(), &summary_config, &translator);
    assert_eq!(data, (count, tags, MEM[10], p_value, fold_change, sample_count, None, edge_mults.clone(), true,  
        "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge coverage: \nA: 1 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge coverage: A: 1, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string(),
        Summarizers::TagsCountsPEM
    )); 

    let data = test_summarize::<IDTagsCountsData, _, _, _>(input_id_tags.clone().into_iter(), &summary_config, &translator);
    assert_eq!(data, (count, tags, MEM[11], p_value, fold_change, sample_count, Some(vec![0, 1, 2, 3, 7, 8]), None, true,
        "IDs: ['0', '1', '2', '3', '7', '8'], samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(), 
        "IDs: ['0', '1', '2', '3', '7', '8'], samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405".to_string(),
        Summarizers::IDTagsCounts
    )); 

    // size: 8o16 + 6 * 4 + 8+8 + 6 * 2o4 + 8+8 + 4 + 8+8+8*1
    //       tag    counts        ids             p   ec
    //       
    //       b       b       b       b       b       b       b       b       b       b       b       b       b       b       b       b       b       b       b       b
    //       |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    //       |tag    |counts box     |counts                 |ids box        |ids        |p  |em box         |em     | 13 * 8 = 104
    //       |tag            |counts box     |counts                 |ids box        |ids        |p  |em box         |em     |buffer | 7 * 16 = 112
    //       |tag    |counts box     |counts                 |ids box        |ids                    |p  |em box         |em     |buf| 14 * 8 = 112
    //       |tag            |counts box     |counts                 |ids box        |ids                    |p  |em box         |em     |buf| 15 * 8 = 120

    let data = test_summarize::<IDTagsCountsPEMData, _, _, _>(input_id_tags.clone().into_iter(), &summary_config, &translator);
    assert_eq!(data, (count, tags, MEM[12], p_value, fold_change, sample_count, Some(vec![0, 1, 2, 3, 7, 8]), edge_mults.clone(), true,  
        "IDs: ['0', '1', '2', '3', '7', '8'], samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge coverage: \nA: 1 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        "IDs: ['0', '1', '2', '3', '7', '8'], samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge coverage: A: 1, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string(),
        Summarizers::IDTagsCountsPEM
    )); 

    let data = test_summarize::<IDEMData, _, _, _>(input_id_tags.clone().into_iter(), &summary_config, &translator);
    assert_eq!(data, (None, None, MEM[1] + 8*4, None, None, None, Some(vec![0, 1, 2, 3, 7, 8]), edge_mults.clone(), true, 
        "IDs: ['0', '1', '2', '3', '7', '8'], edge coverage: \nA: 1 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        "IDs: ['0', '1', '2', '3', '7', '8'], edge coverage: A: 1, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string(), 
        Summarizers::IDEM)
    );

    let data = test_summarize::<IDMapEMData, _, _, _>(input_id_tags.clone().into_iter(), &summary_config, &translator);
    assert_eq!(data, (None, None, MEM[1] + 8*4 + 2*8, None, None, None, Some(vec![0, 1, 2, 3, 7, 8]), edge_mults.clone(), true, 
        "IDs: ['0', '1', '2', '3', '7', '8'], mapped IDs: [], edge coverage: A: 1 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        "IDs: ['0', '1', '2', '3', '7', '8'], mapped IDs: [], edge coverage: A: 1, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string(), 
        Summarizers::IDMapEM)
    );

    let data = test_summarize::<IDMapEMQualityData, _, _, _>(input_id_tags.clone().into_iter(), &summary_config, &translator);
    assert_eq!(data, (None, None, MEM[1] + 8*4 + 2*8 + 8, None, None, None, Some(vec![0, 1, 2, 3, 7, 8]), edge_mults.clone(), true, 
        "IDs: ['0', '1', '2', '3', '7', '8'], mapped IDs: [], quality: medium, edge coverage: A: 1 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        "IDs: ['0', '1', '2', '3', '7', '8'], mapped IDs: [], quality: medium, edge coverage: A: 1, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string(), 
        Summarizers::IDMapEMQuality)
    );

    let data = test_summarize::<GroupCountData, _, _, _>(input_tags.into_iter(), &summary_config, &translator);
    assert_eq!(data, (count, None, MEM[13], None, None, None, None, None, true, 
        "count 1: 4\ncount 2: 2".to_string(), "count 1: 4, count 2: 2".to_string(), Summarizers::GroupCount));

    let data = test_summarize::<RelCountData, _, _, _>(input_tags.into_iter(), &summary_config, &translator);
    assert_eq!(data, (count, None, MEM[14], None, None, None, None, None, true, 
        "relative amount group 1: 66\ncount both: 6".to_string(), "relative amount group 1: 66, count both: 6".to_string(), Summarizers::RelCount));

    // test with different group frax settings
    let summary_config = SummaryConfig::new(sample_info.clone()).with_group_frac(debruijn::summarizer::GroupFrac::One, 0.33);
    let data = test_summarize::<TagsCountsPEMData, _, _, _>(input_tags.into_iter(), &summary_config, &translator);
    assert_eq!(data, (count, tags, MEM[10], p_value, fold_change, sample_count, None, edge_mults.clone(), true,  
        "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge coverage: \nA: 1 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge coverage: A: 1, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string(),
        Summarizers::TagsCountsPEM
    )); 

    let summary_config = SummaryConfig::new(sample_info.clone()).with_group_frac(debruijn::summarizer::GroupFrac::Both, 0.33);
    let data = test_summarize::<TagsCountsPEMData, _, _, _>(input_tags.into_iter(), &summary_config, &translator);
    assert_eq!(data, (count, tags, MEM[10], p_value, fold_change, sample_count, None, edge_mults.clone(), false,  
        "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge coverage: \nA: 1 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge coverage: A: 1, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string(), 
        Summarizers::TagsCountsPEM
    )); 

    let summary_config = SummaryConfig::new(sample_info.clone()).with_max_p(Some(0.05)).with_group_frac(debruijn::summarizer::GroupFrac::One, 0.33);
    let data = test_summarize::<TagsCountsPEMData, _, _, _>(input_tags.into_iter(), &summary_config, &translator);
    assert_eq!(data, (count, tags, MEM[10], p_value, fold_change, sample_count, None, edge_mults.clone(), false,  
        "samples              - counts\n0                    - 1\n1                    - 1\n2                    - 1\n3                    - 1\n7                    - 1\n8                    - 1\nsum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge coverage: \nA: 1 | 0\nC: 0 | 0\nG: 0 | 0\nT: 0 | 0\n".to_string(), 
        "samples: ['0', '1', '2', '3', '7', '8'], counts: [1, 1, 1, 1, 1, 1], sum: 6, p-value: 0.39023498, log2(fold change): 5.4498405, edge coverage: A: 1, C: 0, G: 0, T: 0 | A: 0, C: 0, G: 0, T: 0".to_string(),
        Summarizers::TagsCountsPEM
    )); 

}

/// add the alignment buffer to a structure
/// - `size_heap`: contents of boxes/vectors -> are stored separately and do not go into alignment calculation
fn size_aligned(size_stack: usize, size_heap: usize, align: usize) -> usize {
    let buffer = align - (size_stack % align);

    buffer + size_heap + size_stack
}
