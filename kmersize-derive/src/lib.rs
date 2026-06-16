use proc_macro::TokenStream;
use quote::quote;
use syn::{Data, DeriveInput, parse_macro_input};

#[proc_macro_derive(KmerSize)]
pub fn derive(input: TokenStream) -> TokenStream {
    let input: DeriveInput = parse_macro_input!(input);
    let ident = input.ident;
    let data = input.data;

    let output = match data {
        Data::Struct(_s) => {
            /*
            impl KmerSize for K128 {
                #[inline(always)]
                fn K() -> usize {
                    128
                }
            }
             */

            let k_string = format!("{ident}");
            let mut k = k_string.chars();
            let _ = k.next(); // skip first char ('K')
            let k = k.as_str().parse::<usize>().unwrap();

            let out = quote!{
                impl KmerSize for #ident {
                    #[inline(always)]
                    fn K() -> usize {
                        #k
                    }
                }
            };

            println!("{}", out);
            out
        }
        _ => todo!()
    };
    output.into()
}
