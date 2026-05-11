macro_rules! recognise_tree2 {
    ("birch") => { println!("is burch") };
    ("oak") => { println!("is oak") };
    ($($other:tt)*) => { println!("dunno") };
}

macro_rules! recognise_tree {
    (larch) => { println!("#1, the Larch.") };
    (redwood) => { println!("#2, the Mighty Redwood.") };
    (fir) => { println!("#3, the Fir.") };
    (chestnut) => { println!("#4, the Horse Chestnut.") };
    (pine) => { println!("#5, the Scots Pine.") };
    ($($other:tt)*) => { println!("I don't know; some kind of birch maybe?") };
}

#[cfg(test)]
mod tests1 {

    #[test]
    fn test_macros() {
        recognise_tree!(birch);
        recognise_tree2!(birch);
    }
}

use summarydata_derive::MyTrait;

trait MyTrait {
    fn a(&self) -> Option<u32> { None }
    fn b(&self) -> Option<f32> { None }
    fn sq_b(&self) -> f32 { 0. }
}

#[derive(MyTrait)]
pub struct MyStructOne {
    a: u32
}

#[derive(MyTrait)]
pub struct MyStructTwo {
    a: u32,
    b: f32
}


#[cfg(test)]
mod tests {
    use crate::macros::{MyStructOne, MyStructTwo, MyTrait};

    #[test]
    fn test_my_trait() {
        let struct1 = MyStructOne {a: 1};
        let struct2 = MyStructTwo {a: 2, b: 7.};

        println!("1a: {:?}, 1b:, {:?}, 2a: {:?}, 2b: {:?}", struct1.a(), struct1.b(), struct2.a(), struct2.b());

        println!("1sq b: {}, 2sq b: {}", struct1.sq_b(), struct2.sq_b());
    }
}




