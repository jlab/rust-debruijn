use std::{collections::HashMap, fmt::Display, iter::Sum, marker::PhantomData};

use log::debug;

use crate::{graph::DebruijnGraph, summarizer::{SummaryConfig, SummaryData, ID, M}, Kmer};
use std::fmt::Debug;

/// mode for coloring nodes in dot files - check compatibility with node data
#[derive(Debug)]
pub enum ColorMode<'a> {
    /// only compatible with [`IDSumData`]
    IDGroups {id_group_ids: &'a HashMap<ID, ID>, n_id_groups: usize},
    /// only compatible with [`IDSumData`]
    IDS {n_ids: usize},
    /// compatible with all [`SummaryData`] containing `Tags`
    SampleGroups,
    /// compatible with all [`SummaryData`] containing `TagsCounts`
    FoldChange
}


/// contains the hues, the markers signifying which tag belongs to which group, 
/// the maximun kmer count and the average kmer count
#[derive(Clone, Debug, PartialEq)]
pub struct Colors<SD: SummaryData<DI>, DI> {
    // 2-bit encoded group associations of labels
    marker0: M,
    marker1: M,
    // factor (slope) for log2(fold change) to hue transformation
    log2_fc_factor: Option<f32>,
    // slope (m) and y intercept (b) for n obs to value transformation
    _log2_nobs_mb: Option<(f32, f32)>,
    // slope (m) and y intercept (b) for log10(p-value)) to saturation transformation
    log10_p_mb: Option<(f32, f32)>,
    // slope (m) and y intercept (b) for log10(edge multiplicity) to pen width transformation
    log10_em_mb: Option<(f32, f32)>,
    // phantom data who 
    phantom_data_sd: PhantomData<SD>,
    phantom_data_di: PhantomData<DI>,
}

// add `'b, 'a: 'b` in case of lifetime error
impl<SD: SummaryData<DI> + Debug, DI> Colors<SD, DI> {
    const HUE_RED: f32 = 0.;
    const HUE_YELLOW: f32 = 60. / 360.;
    const HUE_GREEN: f32 = 120. / 360.;
    const HUE_PURPLE: f32 = 289. / 360.;

    const SAT_MIN: f32 = 0.5;
    const SAT_MAX: f32 = 1.;
    const SAT_DEF: f32 = 1.;

    const VAL_MIN: f32 = 0.7;
    const VAL_MAX: f32 = 1.;
    const VAL_DEF: f32 = 1.;

    const FC_MAX: f32 = 5.;
    const FC_MIN: f32 = -5.;
    const P_MAX: f32 = -4.;

    const EDGE_WIDTH_MAX: f32 = 20.;
    const EDGE_WIDTH_MIN: f32 = 1.;
    const EDGE_WIDTH_DEF: f32 = 1.;


    /// Creates a new [`Colors<SD>`]. 
    pub fn new<K: Kmer>(graph: &DebruijnGraph<K, SD>, summary_config: &SummaryConfig) -> Self {
        
        // fold change factor: fold change of node will be multiplied by this factor to get hue
        // factor is the slope of a linear function
        let log2_fc_factor = match graph.get_node(0).data().p_value(summary_config) {
            Some(_) => {
                let (min_max, _, _) = get_min_max(
                &graph, 
                &|&graph| Box::new(graph
                        .iter_nodes()
                        .map(|node| node.data().fold_change(summary_config).expect("error getting fold change"))
                    )
                );
                debug!("log2_fc min max: {:?}", min_max);

                match min_max {
                    Some((min_fc, max_fc)) => {
                         // make symmetrical
                        let m_fc = if min_fc.abs() > max_fc.abs() {
                            min_fc.abs()
                        } else {
                            max_fc.abs()
                        };

                        // if too large, replace with max fold change
                        let val_fc = if m_fc > Self::FC_MAX { 
                            Self::FC_MAX
                        } else {
                            m_fc
                        };

                        // yellow should be where log2(fc) = 0
                        // division by zero should not happen bc max has to be larger than min
                        // and larger absolute of the two is used
                        Some(Self::HUE_YELLOW / val_fc)
                    },
                    None => None
                }
            },
            None => None
        };

        debug!("log2_fc factor: {:?}", log2_fc_factor);

        // number of observations (nobs)
        let (nobs, _, _) = get_min_max(
        &graph, 
        &|&graph| Box::new(graph
                .iter_nodes()
                .map(|node| node.data().score())
            )
        );
        debug!("nobs min max: {:?}", nobs);

        // calculate m and b for value = m * log10(nobs) + b
        let log2_nobs_mb = match nobs {
            Some((min_nobs, max_nobs)) => {
                let min = min_nobs.log2();
                let max = max_nobs.log2();
                
                let m = (Self::VAL_MAX - Self::VAL_MIN) / (max - min);
                let b = Self::VAL_MIN - m * min;

                Some((m, b))
            }
            None => None
        };

        debug!("log2_nobs m, b: {:?}", log2_nobs_mb);
        
        // calculate m and b for saturation = m * log10(p-value) + b
        let log10_p_mb = match graph.get_node(0).data().p_value(summary_config) {
            Some(_) => {
                // get min and max p-value
                let (min_max, _, _) = get_min_max(
                &graph, 
                &|&graph| Box::new(graph
                        .iter_nodes()
                        .map(|node| node.data().p_value(summary_config).expect("error getting p-value"))
                    )
                );                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
                debug!("p min max: {:?}", min_max);
                match min_max {
                    Some((min_p, max_p)) => {
                        let min = max_p.log10();
                        let max = min_p.log10();

                        // replace "max" of p with P_MAX if too small
                        let max = if max < Self::P_MAX { Self::P_MAX } else { max };
                        
                        // b should be 0 -> p = 1 -> log10(1) = 0 -> VAL_MIN=0
                        let m = (Self::SAT_MAX - Self::SAT_MIN) / (max - min);
                        let b = Self::SAT_MIN - m * min;
        
                        Some((m, b))
                    }
                    None => None
                }
            },
            None => None
        };

        debug!("log10_p m, b: {:?}", log10_p_mb);

        // calculate m and b for pen width = m * log10(edge mults) + b
        let log10_em_mb = match graph.get_node(0).data().edge_mults() {
            Some(_) => {
                //  get min and max
                let (min_max, _, _) = get_min_max(&graph, &|graph| Box::new(graph
                    .iter_nodes()
                    .flat_map(|node| node.data().edge_mults().expect("error getting edge mult").edge_mults())
                    .filter(|element| *element != 0)
                )); 
    
                debug!("edge mults min max: {:?}", min_max);
    
                min_max.map(|(min, max)|{
                    let min = (min as f32).log10();
                    let max = (max as f32).log10();
    
                    // b should be 1 -> em = 1 -> log10(1) = 0 -> EDGE_WITH_MIN=1
                    let m = (Self::EDGE_WIDTH_MAX - Self::EDGE_WIDTH_MIN) / (max - min);
                    let b = Self::EDGE_WIDTH_MIN - m * min;
    
                    (m, b)
                })
            }, 
            None => None
        };

        debug!("log10_em m, b: {:?}", log10_em_mb);

        let (marker0, marker1) = summary_config.get_markers();      

        Colors {
            marker0,
            marker1,
            log2_fc_factor,
            _log2_nobs_mb: log2_nobs_mb,
            log10_p_mb,
            log10_em_mb,
            phantom_data_sd: PhantomData,
            phantom_data_di: PhantomData
        }
    }


    /// get the color for a node in a HSV format
    pub fn node_color(&self, data: &SD, summary_config: &SummaryConfig, outline: bool, color_mode: ColorMode) -> String {

        // get hue
        let hue = match color_mode {
            ColorMode::FoldChange => {
                match self.log2_fc_factor {
                    // if fold change available calculate hue based on log2(fc)
                    Some(fc_factor) => {
                        match data.fold_change(summary_config).unwrap() {
                            Self::FC_MAX..=f32::INFINITY => Self::HUE_GREEN,
                            f32::NEG_INFINITY..=Self::FC_MIN => Self::HUE_RED,
                            fc => fc * fc_factor + Self::HUE_YELLOW
                        }
                    }, 
                    None => Self::HUE_PURPLE
                }
            },
            ColorMode::SampleGroups  => {
                match data.tags_sum() {
                    Some((tag, _)) => {
                        if tag.bit_and(self.marker0) & !tag.bit_and(self.marker1) {
                            // tags are only in marker0 group
                            Self::HUE_GREEN
                        } else if !tag.bit_and(self.marker0) & tag.bit_and(self.marker1) {
                            // tag is only in marker1 group
                            Self::HUE_RED
                        } else if tag.bit_and(self.marker0) & tag.bit_and(self.marker1) {
                            // tag is in both groups
                            Self::HUE_YELLOW
                        } else {
                            // tag is in neither group
                            // should only happen if the tags started with more than three distinct characters
                            // (overflow purple)
                            Self::HUE_PURPLE
                        }
                    },
                    None => Self::HUE_PURPLE
                }
            },
            ColorMode::IDGroups { id_group_ids, n_id_groups } => {
                match data.ids() {
                    Some(ids) => {
                        ids.iter().map(|id| *(id_group_ids.get(id).expect("id was not in HM")) as f32 / (n_id_groups * ids.len()) as f32).sum::<f32>()
                    },
                    None => Self::HUE_PURPLE
                }
            }
            ColorMode::IDS { n_ids } => {
                match data.ids() {
                    Some(ids) => {
                        ids.iter().map(|id| *id as f32 / (n_ids * ids.len()) as f32).sum::<f32>()
                    }
                    None => Self::HUE_PURPLE
                }
            }
        };

        // calculate saturation
        let saturation = match self.log10_p_mb {
            //Some((m, b)) => m * data.p_value(summary_config).expect("error getting p-value").log10() + b, // to broad
            Some((_m, _b)) => if data.p_value(summary_config).expect("error getting p-value") < 0.05 { Self::SAT_MAX } else { Self::SAT_MIN },
            None => Self::SAT_DEF
        };

        // set value as default value
        let value = Self::VAL_DEF;

        // adapt font color to value ( currently always black)
        let font_color= if value <= 0.5 { "white" } else { "black" };

        // set outline (eg if it is in a path)
        let prefix = if outline {
            "black, penwidth=7, fillcolor="
        } else {
            ""
        };

        // return formatted string for color, fillcolor, and fontcolor
        format!("color={prefix}\"{hue} {saturation} {value}\", fontcolor={font_color}")
    }

    /// get the edge width based on the edge multiplicity
    pub fn edge_width(&self, edge_mult: u32) -> f32 {
        match self.log10_em_mb {
            Some((m, b)) => (edge_mult as f32).log10() * m + b,
            None => Self::EDGE_WIDTH_DEF
        }
    }
    
}

/// get the minimum and maximum value of an iterator, and if they are three times as 
/// small/large as the next smallest/largest item, also return the latter
/// filters for -inf and inf with floats
pub fn get_min_max<I: Debug, F, II, N>(iter_struct: &I, iter_value: &F) -> (Option<(N, N)>, Option<N>, Option<N>)
where
    F: Fn(&I) -> Box<II>,
    II: Iterator<Item = N>,
    N: CFilter + Sum + PartialOrd + Copy + Display
{  
    let min_o = iter_value(iter_struct)
        .filter(|value| value.filter())
        .min_by(|a, b| a.partial_cmp(b).unwrap());

    // check if there are actually elements to process
    if let Some(min) = min_o {
        let max = iter_value(iter_struct)
            .filter(|value| value.filter())
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .expect("error: empty iterator");

        if max > min {
            let snd_max = iter_value(iter_struct)
                .filter(|value| *value != max && value.filter())
                .max_by(|a, b| a.partial_cmp(b).unwrap())
                .expect("error: empty iterator");
            let snd_min = iter_value(iter_struct)
                .filter(|value| *value != min && value.filter())
                .min_by(|a, b| a.partial_cmp(b).unwrap())
                .expect("error: empty iterator");
            
            // if max is OUTL times bigger than snd_max, it is an outlier
            const OUTL: f64 = 3.;
            let outlier_max = if (max.to_f64() / snd_max.to_f64()) > OUTL && snd_max > snd_min { Some(snd_max) } else { None };
            let outlier_min = if (min.to_f64() / snd_min.to_f64()) > OUTL && snd_max > snd_min { Some(snd_min) } else { None };

            return (Some((min, max)), outlier_min, outlier_max)
        }
    }   
    
    (None, None, None)
}

/// trait for filtering the values in [`get_min_max`]
pub trait CFilter {
    fn filter(self) -> bool;
    fn to_f64(self) -> f64;
}

impl CFilter for usize {
    fn filter(self) -> bool {
        self != usize::MAX
    }
    fn to_f64(self) -> f64 {
        self as f64
    }
}

impl CFilter for u32 {
    fn filter(self) -> bool {
        self != u32::MAX
    }
    fn to_f64(self) -> f64 {
        self as f64
    }
}

impl CFilter for f32 {
    fn filter(self) -> bool {
        self.is_finite()
    }
    fn to_f64(self) -> f64 {
        self as f64
    }
}

