use std::{fmt::Display, iter::Sum, marker::PhantomData};

use log::debug;

use crate::{graph::DebruijnGraph, summarizer::{SummaryConfig, SummaryData, M}, Kmer};
use std::fmt::Debug;



/// contains the hues, the markers signifying which tag belongs to which group, 
/// the maximun kmer count and the average kmer count
#[derive(Clone)]
pub struct Colors<SD: SummaryData<u8>> {
    marker0: M,
    marker1: M,
    log2_fc_factor: Option<f32>,
    _log2_nobs_mb: Option<(f32, f32)>,
    log10_p_mb: Option<(f32, f32)>,
    phantom_data_sd: PhantomData<SD>,
}

// add `'b, 'a: 'b` in case of lifetime error

impl<SD: SummaryData<u8> + Debug> Colors<SD> {
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

    /// Creates a new `Colors<SD>`. 
    /// 
    /// The markers binary encode the groups of the tags:
    /// `marker0 = 00010011` (shortened to `u8`) means the tags `[0, 1, 4]` are in group 1.
    pub fn new<K: Kmer>(graph: &DebruijnGraph<K, SD>, summary_config: &SummaryConfig) -> Self {
        
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
                // make it symmetrical
                match min_max {
                    Some((min_fc, max_fc)) => {
                        let m_fc = if min_fc.abs() > max_fc.abs() {
                            min_fc.abs()
                        } else {
                            max_fc.abs()
                        };
                        let val_fc = if m_fc > Self::FC_MAX { 
                            Self::FC_MAX
                        } else {
                            m_fc
                        };
                        Some(Self::HUE_YELLOW / val_fc)
                    },
                    None => None
                }
            },
            None => None
        };

        debug!("log2_fc factor: {:?}", log2_fc_factor);

        let (nobs, _, _) = get_min_max(
        &graph, 
        &|&graph| Box::new(graph
                .iter_nodes()
                .map(|node| node.data().score())
            )
        );
        debug!("nobs min max: {:?}", nobs);

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
        
        let log10_p_mb = match graph.get_node(0).data().p_value(summary_config) {
            Some(_) => {
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

        let (marker0, marker1) = summary_config.get_markers();      

        Colors {
            marker0,
            marker1,
            log2_fc_factor,
            _log2_nobs_mb: log2_nobs_mb,
            log10_p_mb,
            phantom_data_sd: PhantomData,
        }
    }


    /// get the color for a node in a HSV format
    pub fn color(&self, data: &SD, summary_config: &SummaryConfig, outline: bool) -> String {

        // get hue
        let hue = match self.log2_fc_factor {
            // if fold change available calculate hue based on log2(fc)
            Some(fc_factor) => {
                match data.fold_change(summary_config).unwrap() {
                    Self::FC_MAX..=f32::INFINITY => Self::HUE_GREEN,
                    f32::NEG_INFINITY..=Self::FC_MIN => Self::HUE_RED,
                    
                    fc => fc * fc_factor + Self::HUE_YELLOW
                }
            }, 
            // if not, try to use discrete colors based on tags
            None => {
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
                
            }
        };

        // get saturation
        /* let saturation = match self.log2_nobs_mb {
            Some((m, b)) => m * data.score().log2() + b,
            None => Self::SAT_DEF
        }; */
        let saturation = match self.log10_p_mb {
            //Some((m, b)) => m * data.p_value(summary_config).expect("error getting p-value").log10() + b,
            Some((_m, _b)) => if data.p_value(summary_config).expect("error getting p-value") < 0.05 { Self::SAT_MAX } else { Self::SAT_MIN },
            None => Self::SAT_DEF
        };

        // get value from p value significant or not
        let value = Self::VAL_DEF;

        // adapt font color to value
        let font_color= if value <= 0.5 { "\"white" } else { "\"black" };

        let prefix = if outline {
            "black, penwidth=7, fillcolor="
        } else {
            ""
        };

        // return formatted string for color, fillcolor, and fontcolor
        format!("color={prefix}\"{hue} {saturation} {value}\", fontcolor={font_color}")
    }
    
}

pub fn get_min_max<I: Debug, F, II, N>(iter_struct: &I, iter_value: &F) -> (Option<(N, N)>, Option<N>, Option<N>)
where
    F: Fn(&I) -> Box<II>,
    II: Iterator<Item = N>,
    N: CFilter + Sum + PartialOrd + Copy + Display
{  
    let min = iter_value(iter_struct)
        .filter(|value| value.filter())
        .min_by(|a, b| a.partial_cmp(b).unwrap())
        .expect("error: empty iterator");
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
        
        let outlier_max = if (max.to_f64() / snd_max.to_f64()) > 3. { Some(snd_max) } else { None };
        let outlier_min = if (min.to_f64() / snd_min.to_f64()) > 3. { Some(snd_min) } else { None };

        return (Some((min, max)), outlier_min, outlier_max)
    }

    (None, None, None)
}

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

impl CFilter for f32 {
    fn filter(self) -> bool {
        self.is_finite()
    }
    fn to_f64(self) -> f64 {
        self as f64
    }
}

