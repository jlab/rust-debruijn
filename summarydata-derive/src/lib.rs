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
pub fn derive(input: TokenStream) -> TokenStream {
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