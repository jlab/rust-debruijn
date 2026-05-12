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
            let has_counts = fields.iter().filter(|f| **f == "counts").next().is_some();
            let has_p_value = fields.iter().filter(|f| **f == "p_value").next().is_some();
            let has_edge_mults = fields.iter().filter(|f| **f == "edge_mults").next().is_some();
            let has_ids = fields.iter().filter(|f| **f == "ids").next().is_some();
            let has_map_ids = fields.iter().filter(|f| **f == "map_ids").next().is_some();
            let has_edge_maps = fields.iter().filter(|f| **f == "edge_maps").next().is_some();
            let has_quality = fields.iter().filter(|f| **f == "quality").next().is_some();

            // conditional implementations
            // TODO
            // print in multiple lines
            // print in one line
            // print for json

            // tags for in vec // TODO check if this works, if not, check for capacity field instead
            let tags_from_vec = if ident == "Vec<Tag>" {
                Some(
                    quote! {
                        fn tags(&self) -> Option<Tags> { 
                            Some(Tags::from_tag_vec(self.clone()))
                        }
                    }
                )   
            } else { None };

            // sum for summarizers which only have counts (no sum)
            let counts_sum = if !has_sum {
                fields.iter().filter(|f| **f == "counts").next().map(|_f| {
                    quote! {
                        fn sum(&self) -> Option<u32> {
                            Some(self.sum())
                        }
                    }
                })
            } else { None };

            // include heap memory for
            // Vec<Tag>
            let heap_vec_tag = if ident == "Vec<Tag>" { 
                Some( quote! { + mem::size_of_val(&**self) } )   
            } else { None };
            // ids
            let heap_ids = fields.iter().filter(|f| **f == "ids").next().map(|_f| {
                quote! { + mem::size_of_val(&*self.ids) }
            });
            // map_ids
            let heap_map_ids = fields.iter().filter(|f| **f == "map_ids").next().map(|_f| {
                quote! { + mem::size_of_val(&*self.map_ids) }
            });
            // counts
            let heap_counts = fields.iter().filter(|f| **f == "counts").next().map(|_f| {
                quote! { + mem::size_of_val(&*self.counts) }
            });
            // edge_maps
            let heap_edge_maps = fields.iter().filter(|f| **f == "edge_maps").next().map(|_f| {
                quote! { + self.edge_maps.mem_heap() }
            });

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
            if has_tags {
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
            } else { None };

            // add validity checks for each possible field, base starts with true in case no other checks
            // sample counts and sum are valid
            let valid_counts = if has_tags & (has_counts | has_sum) {
                // we have tags and counts or sum, use for valid_counts
                Some(quote! {
                    && valid_counts(self.tags, Some(self.sum()), config)
                })
            } else if has_tags {
                Some(quote! {
                    && valid_counts(self.tags, None, config)
                })
            } else if has_counts | has_sum {
                Some(quote! {
                    && self.sum().expect("missing sum") as usize >= config.min_kmer_obs
                })
            } else { None };
            
            // p value is valid
            let valid_p = if has_p_value {
                Some(quote! {
                    && valid_p(PInfo::PValue { p: self.p_value(config).expect("error getting p-values") }, config)
                })
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
            let valid_q = fields.iter().filter(|f| **f == "edge_mults").next().map(|_f| {
            quote! { && self.quality >= config.min_quality }
            });

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
            let edge_maps = fields.iter().filter(|f| **f == "edge_maps").next().map(|_f| { 
                quote! {
                    fn mapped_edge_ids(&self) -> Option<&EdgeMap> {
                        Some(&self.edge_maps)
                    }

                    fn set_mapped_edge_ids(&mut self, mapped_edge_ids: Option<EdgeMap>) {
                        self.edge_maps = mapped_edge_ids.expect("Error: no mapped edge IDs")
                    }
                }
            });

            // fix edge data
            let fix_edge_data = if has_edge_mults & has_edge_maps {
                Some(quote! {
                    fn fix_edge_data(&mut self, exts: Exts) {
                        self.edge_mults.clean_edges(exts);
                        self.edge_maps.clean_edges(exts);
                    }
                })
            } else if has_edge_mults {
                Some(quote! {
                    fn fix_edge_data(&mut self, exts: Exts) {
                        self.edge_mults.clean_edges(exts);
                    }
                })
            } else if has_edge_maps {
                Some(quote! {
                    fn fix_edge_data(&mut self, exts: Exts) {
                        self.edge_maps.clean_edges(exts);
                    }
                })
            } else { None };

            // ids
            let ids = fields.iter().filter(|f| **f == "ids").next().map(|_f| { 
                quote! {
                    fn ids(&self) -> Option<&[ID]> {
                        Some(&self.ids)
                    }
                }
            });

            // mapped ids
            let mapped_ids = fields.iter().filter(|f| **f == "map_ids").next().map(|_f| { 
                quote! {
                    fn mapped_ids(&self) -> Option<&[ID]> {
                        Some(&self.map_ids)
                    }

                    fn set_mapped_ids(&mut self, mapped_ids: Box<[ID]>) {
                        self.map_ids = mapped_ids
                    }
                }
            });

            // TODO join test
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

            // TODO summarize

            // summarizer
            let summarizer = if ident == "u32" {
                quote! {Sum}
            } else if ident == "Vec<u32>" {
                quote! {VecTags}
            } else {
                quote! {#ident}
            };



            // TODO exclude certain fields from writing getters for them
            // edge mults, counts, edge maps, ids -> all getters as reference

            // final implementation for all
            quote! {
                impl SummaryData for #ident {
                    #(
                        fn #fields(&self) -> Option<#types> {
                            Some(self.#fields)
                        }
                    )*

                    #tags_from_vec
                    #counts_sum

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

                    fn summarizer() -> Summarizers {
                        Summarizers::#summarizer
                    }

                }
            }
        }
        _ => todo!()
    };
    output.into()
}