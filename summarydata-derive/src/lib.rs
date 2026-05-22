/* pub fn add(left: u64, right: u64) -> u64 {
    left + right
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }
} */


use std::collections::HashSet;

use proc_macro::{self, TokenStream};
use quote::quote;
use syn::{Data, DeriveInput, parse_macro_input};

#[proc_macro_derive(MyTrait)]
pub fn derive_mytrait(input: TokenStream) -> TokenStream {
    let input: DeriveInput = parse_macro_input!(input);
    let ident = input.ident;
    let data = input.data;

    let output = match data {
        Data::Struct(s) => {
            let fields  = s.fields.iter().map(|f| f.ident.clone().unwrap()).collect::<Vec<_>>(); 
            let types = s.fields.iter().map(|f| f.ty.clone()).collect::<Vec<_>>(); 

            // conditional implementations
            // square_b
            let square_b = if let Some(_f) = fields.iter().filter(|f| **f == "b").next() {
                Some(quote! {
                    fn sq_b(&self) -> f32 {
                        self.b * self.b
                    }
                })
            } else {
                None
            };

            quote! {
                impl MyTrait for #ident {
                    #(
                        fn #fields(&self) -> Option<#types> {
                            Some(self.#fields)
                        }
                    )*

                    #square_b
                }
            }
        }
        _ => todo!()
    };
    output.into()
}


#[proc_macro_derive(SummaryData)]
pub fn derive(input: TokenStream) -> TokenStream {
    let input: DeriveInput = parse_macro_input!(input);
    let ident = input.ident;
    let data = input.data;

    let output = match data {
        Data::Struct(s) => {
            let fields  = s.fields.iter().map(|f| f.ident.clone().unwrap()).collect::<Vec<_>>(); 
            let types = s.fields.iter().map(|f| f.ty.clone()).collect::<Vec<_>>(); 

            // check presence of fields
            let has_sum = fields.iter().filter(|f| **f == "sum").next().is_some();
            let has_tags = fields.iter().filter(|f| **f == "tags").next().is_some();
            let has_tag_vec  = fields.iter().filter(|f| **f == "buf").next().is_some();
            let has_counts = fields.iter().filter(|f| **f == "counts").next().is_some();
            let has_p_value = fields.iter().filter(|f| **f == "p_value").next().is_some();
            let has_edge_mults = fields.iter().filter(|f| **f == "edge_mults").next().is_some();
            let has_ids = fields.iter().filter(|f| **f == "ids").next().is_some();
            let has_map_ids = fields.iter().filter(|f| **f == "map_ids").next().is_some();
            let has_edge_maps = fields.iter().filter(|f| **f == "edge_maps").next().is_some();
            let has_quality = fields.iter().filter(|f| **f == "quality").next().is_some();
            let has_groups = fields.iter().filter(|f| **f == "group1").next().is_some()
                && fields.iter().filter(|f| **f == "group2").next().is_some();
            let has_percent = fields.iter().filter(|f| **f == "percent").next().is_some();

            // conditional implementations
            // format (translate) IDs, mapped IDs and tags
            let id_format = if has_ids { Some(quote! {id_format(&self.ids, translator, id_group_translator)}) } else { None };
            let id_format_ol = id_format.as_ref().map(|f| quote! {string.push_str(&format!("IDs: {}, ", #f));});
            let id_format_json = id_format.as_ref().map(|f| quote! {string.push_str(&format!("\"ids\": {}, ", #f));});
            
            let map_id_format = if has_map_ids { Some(quote! {id_format(&self.map_ids, translator, id_group_translator)}) } else { None };
            let map_id_format_ol = map_id_format.as_ref().map(|f| quote! {string.push_str(&format!("mapped IDs (node): {}, ", #f));});
            let map_id_format_json = map_id_format.as_ref().map(|f| quote! {string.push_str(&format!("\"mapped_ids_nodes\": {}, \"has_mapped_ids\": {}, ", #f, !self.map_ids.is_empty() as usize));});
            
            let tag_format_ol = fields.iter().filter(|f| **f == "tags").next().map(|_|
                quote! {
                    let tags = if let Some(tag_translator) = translator.tag_translator() {
                        format!("samples: {:?}, ", self.tags.to_string_vec(tag_translator))
                    } else {
                        format!("samples: {:?}, ", self.tags.to_tag_vec())
                    };
                    string.push_str(&tags);
                }
            );
            let tag_format_json = fields.iter().filter(|f| **f == "tags").next().map(|_|
                quote! {
                    let tags = if let Some(tag_translator) = translator.tag_translator() {
                        format!("\"samples\": {:?}, ", self.tags.to_string_vec(tag_translator))
                    } else {
                        format!("\"samples\": {:?}, ", self.tags.to_tag_vec())
                    };
                    string.push_str(&tags);
                }
            );

            let sum_format = if has_sum { Some(quote! {self.sum}) } else if has_counts {Some(quote! {self.sum()})} else { None };
            let sum_format_ol = sum_format.as_ref().map(|s| Some(quote! { string.push_str(&format!("sum: {}, ", #s)); }));
            let sum_format_json = sum_format.as_ref().map(|s| Some(quote! { string.push_str(&format!("\"sum\": {}, ", #s)); }));

            let quality_format_json = if has_quality { Some(quote! { string.push_str(&format!("\"quality\": {}, ", self.quality as usize)); })} else { None };

            let stats_format = if has_tags & has_counts { Some(quote! { 
                let p = match self.p_value(config) {
                    Some(p) => format!("{}", p),
                    None => "".to_string()
                };

                let fc = match self.fold_change(config) {
                    Some(fc) => format!("{}", fc),
                    None => "".to_string()
                };
            }) } else { None };
            let stat_format_ol = stats_format.as_ref().map(|sf| quote! { 
                #sf
                string.push_str(&format!("p-value: {}, log2(fold change): {}, ", p, fc));
            });
            let stat_format_json = stats_format.as_ref().map(|sf| quote! { 
                #sf
                string.push_str(&format!("\"p_value\": {}, \"fold_change\": {}, ", p, fc)); 
            });

            // print in potentially multiple lines
            let print = {
                // edge_mults, edge_maps, ids, map_ids, tags are formatted separately
                // vec fields len and buf should not be included
                const NOPRINT_FIELDS: [&str; 9] = ["edge_mults", "edge_maps", "ids", "map_ids", "buf", "len", "tags", "sum", "p_value"];
                let print_fields = fields.iter().filter(|f| NOPRINT_FIELDS.iter().filter(|invf| f == invf).next().is_none()).collect::<Vec<_>>();

                quote! {
                    fn print(&self, translator: &Translator, config: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
                        let mut string = String::new();
                        #id_format_ol
                        #map_id_format_ol
                        #tag_format_ol
                        #(
                            string.push_str(&format!("{}: {:?}, ", stringify!(#print_fields), self.#print_fields));
                        )*
                        #sum_format_ol
                        #stat_format_ol
                        // remove last two characters ", "
                        string.pop();
                        string.pop();
                        string.replace("\"", "\'")
                    }
                }
            };
            // print in one line
            let print_ol = {
                // ids, map_ids, tags are formatted separately
                // vec fields len and buf should not be included
                const NOPRINT_FIELDS: [&str; 7] = ["ids", "map_ids", "buf", "len", "tags", "sum", "p_value"];
                let print_fields = fields.iter().filter(|f| NOPRINT_FIELDS.iter().filter(|invf| f == invf).next().is_none()).collect::<Vec<_>>();

                quote! {
                    fn print_ol(&self, translator: &Translator, config: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
                        let mut string = String::new();
                        #id_format_ol
                        #map_id_format_ol
                        #tag_format_ol
                        #(
                            string.push_str(&format!("{}: {:?}, ", stringify!(#print_fields), self.#print_fields));
                        )*
                        #sum_format_ol
                        #stat_format_ol
                        // remove last two characters ", "
                        string.pop();
                        string.pop();
                        string.replace("\"", "\'")
                    }
                }
            };
            // print for json
            let print_json = {
                // ids, map_ids, tags are formatted separately
                // vec fields len and buf should not be included
                const NOPRINT_FIELDS: [&str; 10] = ["ids", "map_ids", "buf", "len", "tags", "sum", "p_value", "edge_mults", "quality", "edge_maps"];
                let print_fields = fields.iter().filter(|f| NOPRINT_FIELDS.iter().filter(|invf| f == invf).next().is_none()).collect::<Vec<_>>();

                quote! {
                    fn print_json(&self, translator: &Translator, config: &SummaryConfig, id_group_translator: Option<&HashMap<ID, ID>>) -> String {
                        let mut string = String::new();
                        #id_format_json
                        #map_id_format_json
                        #tag_format_json
                        #(
                            string.push_str(&format!("\"{}\": {:?}, ", stringify!(#print_fields), self.#print_fields));
                        )*
                        #quality_format_json
                        #sum_format_json
                        #stat_format_json
                        // remove last two characters ", "
                        string.pop();
                        string.pop();
                        string
                    }
                }
            };

            // tags for in vec // TODO check if this works, if not, check for capacity field instead
            let tags_from_vec: Option<proc_macro2::TokenStream> = if ident == "Vec<Tag>" {
                Some(
                    quote! {
                        fn tags(&self) -> Option<Tags> { 
                            Some(Tags::from_tag_vec(self.clone()))
                        }
                    }
                )   
            } else { None };

            // sum for summarizers which only have counts or group sums (no regular sum)
            // have to have a separate Self::sum() impl
            let counts_sum = if has_sum { None } else if has_counts | has_groups {
                Some(quote! {
                        fn sum(&self) -> Option<u32> {
                            Some(self.sum())
                        }
                    }
                )
            } else { None };

            // include heap memory for
            let heap_vec_tag = if ident == "Vec<Tag>" { Some(quote! { + mem::size_of_val(&**self) }) } else { None }; // TODO does this work?
            let heap_ids = if has_ids { Some(quote! { + mem::size_of_val(&*self.ids) }) } else {None};
            let heap_map_ids = if has_map_ids { Some(quote! { + mem::size_of_val(&*self.map_ids) }) } else {None};
            let heap_counts = if has_counts { Some(quote! { + mem::size_of_val(&*self.counts) }) } else {None};
            let heap_edge_maps = if has_edge_maps { Some(quote! { + self.edge_maps.mem_heap() }) } else {None};

            // get or calculate the p-value
            let p_value = if has_p_value {
                // we have a p-value, use it directly or recalculate if the stat test has changed
                Some(quote! {
                    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {
                        if config.stat_test_changed {
                            Some(p_value(&self.tags.to_tag_vec(), &self.counts, config).unwrap())
                        } else {
                            Some(self.p_value)
                        } 
                    }
                })
            } else if has_tags & has_counts {
                // we have the data to calculate a p-value
                Some(quote! {
                    fn p_value(&self, config: &SummaryConfig) -> Option<f32> {
                        p_value(&self.tags.to_tag_vec(), &self.counts, config).ok()
                    }
                })
            } else { None };

            // calculate the fold change
            let fold_change = if has_tags & has_counts {
                // we have the data to calculate a p-value
                Some(quote! {
                    fn fold_change(&self, config: &SummaryConfig) -> Option<f32> {
                        Some(log2_fold_change(self.tags, &self.counts, &config.sample_info))
                    }
                })
            } else { None };

            // add sample count
            let sample_count = if has_tags {
                Some(quote! {
                    fn sample_count(&self) -> Option<usize> {
                        Some(self.tags.len())
                    }
                })
            } else if has_counts {
                Some(quote! {
                    fn sample_count(&self) -> Option<usize> {
                        Some(self.counts.len())
                    }
                })
            } else if has_tag_vec {
                Some(quote! {
                    fn sample_count(&self) -> Option<usize> {
                        Some(self.len())
                    }
                })
            }else { None };

            // add validity checks for each possible field, base starts with true in case no other checks
            // sample counts and sum are valid
            let valid_counts = if has_tags & has_sum {
                // we have tags and sum, use for valid_counts
                Some(quote! {&& valid_counts(self.tags, Some(self.sum), config)})
            } else if has_tags & has_counts {
                // we have tags and counts, use for valid_counts
                Some(quote! {&& valid_counts(self.tags, Some(self.sum()), config)})
            } else if has_tags {
                Some(quote! {&& valid_counts(self.tags, None, config)})
            } else if has_counts | has_sum {
                Some(quote! {&& self.sum().expect("missing sum") as usize >= config.min_kmer_obs})
            } else { None };
            
            // p value is valid
            let valid_p = if has_p_value {
                Some(quote! {&& valid_p(PInfo::PValue { p: self.p_value(config).expect("error getting p-values") }, config)})
            } else if has_counts & has_tags {
                // we have to calculate p // TODO hope this works
                Some(quote! {
                    && match config.max_p {
                        Some(p) => self.p_value(config).expect("error calculating p-value") <= p,
                        None => true,
                    }

                })
            } else { None };

            // quality is valid
            let valid_q = if has_quality { Some(quote! { && self.quality >= config.min_quality }) }  else { None };

            // edge mult methods
            let edge_mults = fields.iter().filter(|f| **f == "edge_mults").next().map(|_f| {
                quote! {
                    fn edge_mults(&self) -> Option<&EdgeMult> {
                        Some(&self.edge_mults)
                    }

                    fn set_edge_mults(&mut self, edge_mults: Option<EdgeMult>) {
                        self.edge_mults = edge_mults.expect("Error: no edge mults")
                    }
                }
            });

            // edge map methods
            let edge_maps = if has_edge_maps {
                Some(quote! {
                    fn mapped_edge_ids(&self) -> Option<&EdgeMap> {
                        Some(&self.edge_maps)
                    }

                    fn set_mapped_edge_ids(&mut self, mapped_edge_ids: Option<EdgeMap>) {
                        self.edge_maps = mapped_edge_ids.expect("Error: no mapped edge IDs")
                    }
                })
            } else { None };

            // fix edge data
            let fix_edge_data = if has_edge_mults | has_edge_maps {

                let edge_data_inner = if has_edge_maps & has_edge_mults {
                    quote!{
                        self.edge_mults.clean_edges(exts);
                        self.edge_maps.clean_edges(exts);
                    }
                } else if has_edge_mults {
                    quote! {self.edge_mults.clean_edges(exts);}
                } else {
                    quote! { self.edge_maps.clean_edges(exts); }
                };

                Some(quote! {
                    fn fix_edge_data(&mut self, exts: Exts) {
                        #edge_data_inner
                    }
                })
            } else { None };

            // ids
            let ids = if has_ids {
                Some(quote! {
                    fn ids(&self) -> Option<&[ID]> {
                        Some(&self.ids)
                    }
                })
            } else { None };

            // mapped ids
            let mapped_ids = if has_map_ids {
                Some(quote! {
                    fn mapped_ids(&self) -> Option<&[ID]> {
                        Some(&self.map_ids)
                    }

                    fn set_mapped_ids(&mut self, mapped_ids: Box<[ID]>) {
                        self.map_ids = mapped_ids
                    }
                })
            } else { None };

            // join test
            // edge_maps and edge_mults are ignored for join test
            let join_test = if has_edge_maps | has_edge_mults {
                let valid_fields = fields.iter().filter(|f| **f != "edge_maps" && **f != "edge_mults").collect::<Vec<_>>();

                Some(quote! {
                    fn join_test(&self, other: &Self) -> bool {
                        true
                        #(
                            && self.#valid_fields == other.#valid_fields
                        )*
                    }
                })
            } else { None };

            // summarize

            // Tag or IDTag?
            let summary_item = if has_ids { quote! {IDTag} } else {quote! {Tag}};
            // summarize items with or without quality
            let summary = if has_ids {
                quote! {let summary = summarize_tags_ids_edge_q(items, config);}
            } else {
                quote! {let summary = summarize_tags_edge_q(items, config);}
            };

            // is p-value valid?
            let summary_valid_p = if has_p_value {
                quote! {
                    // calculate p-value with chosen test
                    let p_value = p_value(&summary.tag_vec, &summary.tag_counts, config).unwrap();
                    let valid_p = valid_p(PInfo::PValue { p: p_value }, config);
                }
            } else {
                quote! {let valid_p = valid_p(PInfo::Calculate { tag_vec: &summary.tag_vec, tag_counts: &summary.tag_counts}, config);}
            };

            // is quality valid?
            let summary_valid_q = if has_quality {
                quote! {
                    let quality = summary.highest_quality.expect("missing quality score - required for summarizer");
                    let valid_q = quality >= config.min_quality;
                }
            } else {
                quote! {let valid_q = if let Some(q) = summary.highest_quality { q >= config.min_quality } else {true };}
            };

            // rename and initialize values so we can construct with just the field names
            let summarized_ids = if has_ids { Some(quote! {let ids: Box<[ID]> = summary.id_vec.into();}) } else { None };
            let summarized_edge_mults = if has_edge_mults { Some(quote! {let edge_mults = summary.edge_mults;}) } else { None };
            let summarized_edge_maps = if has_edge_maps { Some(quote! {let edge_maps = EdgeMap::default();}) } else { None };
            let summarized_map_ids = if has_map_ids { Some(quote! {let map_ids = Vec::new().into();}) } else { None };
            let summarized_sum = if has_sum { 
                Some(quote! {
                    let sum = match config.significant {
                        Some(digits) => round_digits(summary.sum, digits),
                        None => summary.sum  
                    };
                })
            } else { None };
            let summarized_groups = if has_groups | has_percent {
                Some(quote! {
                    let mut count1 = 0;
                    let mut count2 = 0;

                    for (count, tag) in counts.iter().zip(summary.tag_vec.clone()) {
                        let bin_tag = (2 as Marker).pow(tag as u32);
                        let group1 = ((config.sample_info.marker0 & bin_tag) > 0) as u32;
                        let group2 = ((config.sample_info.marker1 & bin_tag) > 0) as u32;

                        if (group1 + group2) != 1 { 
                            panic!(
                                "should not happen\n tag: {:#066b}\n m1:  {:#066b}\n m2:  {:#066b}\n g1:  {}\n g2:  {}", 
                                bin_tag, config.sample_info.marker0, config.sample_info.marker1, group1, group2
                            )
                        }
                        count1 += group1 * count;
                        count2 += group2 * count;
                    }

                    let (group1, group2) = match config.significant {
                        Some(digits) => (round_digits(count1, digits), round_digits(count2, digits)),
                        None => (count1, count2)
                    };
                })
            } else { None };
            let summarized_percent = fields.iter().filter(|f| **f == "percent").next().map(|_f| { 
                quote! {let percent = (group1 as f32 / (group1 + group2) as f32 * 100.) as u32;}
            });

            // base summarize method
            let summarize = quote! {
                fn summarize<K: Kmer, F: Iterator<Item = KmerDataItem<K, #summary_item>>>(items: F, config: &SummaryConfig) -> (bool, Exts, Self) {
                    #summary
                    #summary_valid_p
                    #summary_valid_q

                    let counts: Box<[u32]> = summary.tag_counts.into();

                    #summarized_ids
                    #summarized_edge_mults
                    #summarized_edge_maps
                    #summarized_map_ids
                    #summarized_sum
                    #summarized_groups
                    #summarized_percent

                    let tags = Tags::from_tag_vec(summary.tag_vec);
                    let valid  = valid_counts(tags, Some(summary.sum), config) && valid_p && valid_q;

                    (valid && valid_p, summary.all_exts, #ident { #(#fields, )* })
                }
            };

            // summarizer
            let summarizer = if ident == "u32" {
                quote! {Sum}
            } else if ident == "Vec<u32>" {
                quote! {VecTags}
            } else {
                quote! {#ident}
            };

            // exclude certain fields from writing getters for them
            // edge mults, edge maps, ids, map_ids -> all getters as reference
            // p_value -> more complicated getter
            // group1, group2, percent, counts -> do not get getters
            // buf, len -> fields of Vec<u32> -> also no getters
            const INVALID_FIELDS: [&str; 11] = ["edge_mults", "edge_maps", "ids", "map_ids", "group1", "group2", "percent", "counts", "buf", "len", "p_value"];

            let (getter_fields, getter_types) = fields.iter().zip(&types).filter(|(f, _t)| INVALID_FIELDS.iter().filter(|invf| f == invf).next().is_none()).collect::<(Vec<_>, Vec<_>)>();

            // final implementation for all
            let sd = quote! {
                impl SummaryData<#summary_item> for #ident {

                    #print
                    #print_ol
                    #print_json
                    
                    #(
                        fn #getter_fields(&self) -> Option<#getter_types> {
                            Some(self.#getter_fields)
                        }
                    )*

                    #tags_from_vec
                    #counts_sum
                    #sample_count

                    fn mem(&self) -> usize {
                        mem::size_of_val(self) #heap_vec_tag #heap_ids #heap_map_ids #heap_counts #heap_edge_maps
                    }

                    #p_value
                    #fold_change

                    fn valid(&self, config: &SummaryConfig) -> bool {
                        true #valid_counts #valid_p #valid_q
                    }

                    #edge_mults 
                    #edge_maps
                    #fix_edge_data
                    #ids
                    #mapped_ids
                    #join_test
                    #summarize

                    fn summarizer() -> Summarizers {
                        Summarizers::#summarizer
                    }

                }
            };
        println!("{}", sd);
        sd
        }
        _ => todo!()
    };
    output.into()
}