use serde_derive::{Deserialize, Serialize};
use std::fmt::Debug;
use std::hash::Hash;
use std::usize;
use crate::dna_string::{DnaString, PackedDnaStringSet};
use crate::Exts;


/// Reads Structure: Option 1
/// storage: a constant no of u64s for each read, 2-bit encoded, space for 32 bases in each u64 (len n * 5 probably u64)
/// lens: actual length of read (rest of last u64 is 00) (len n u32) (could be u16)
/// max_len: maximum length of read, space which is reserved for each read - eg. 151 bases -> 5 u64s are reserved
/// exts: exts (probably empty) for each read (len n u8)
/// data: data (label/color) for each read (len n (u8))
/// 
/// would require methods to increads max_length -> add 0u64s after every read
#[derive(Ord, PartialOrd, Clone, PartialEq, Eq, Hash, Serialize, Deserialize, Debug)]
pub struct ReadsConst<D> {
    storage: Vec<u64>,
    lens: Vec<u32>,
    max_len: usize,
    exts: Vec<Exts>, 
    data: Vec<D>,
}

/// Packed Reads Structure: Option 2
/// seqs: PackedDnaStringSet containing all reads, packed in u64s, their positions and lengths
///     contains one Vec (len n * 4.2 avg u64) for reads as u64s , a Vec (len n u32) for lengths and a Vec (len n usize) for Positions
/// exts: exts (probably empty) for each read (len n u8)
/// data: data (label/color) for each read (len n (u8))
#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct ReadsPacked<D> {
    seqs: PackedDnaStringSet,
    exts: Vec<Exts>, 
    data: Vec<D>,
}

/// Unpacked Packed Reads Structure: Option 3
/// (basically same as above)
/// storage + lens + start: basically PackedDnaString
/// exts: exts (probably empty) for each read (len n u8)
/// data: data (label/color) for each read (len n (u8))
#[derive(Ord, PartialOrd, Clone, PartialEq, Eq, Hash, Serialize, Deserialize, Debug)]
pub struct ReadsUnPacked<D> {
    storage: Vec<u64>,
    lens: Vec<u32>,
    start: Vec<usize>,
    exts: Vec<Exts>, 
    data: Vec<D>,
}

/// Var Reads Structure: Option 4
/// storage: reads as u64s, also packed (len n * 4.2 probably u64)
/// ends: end position of each read in the u64 encoding (e.g. 120 would be somewhere in second u64)
///     length of each read can be calculated with distance to previous number, index by dividing by 32 and % (len n usize)
///     exceptions need to made for first read
/// exts: exts (probably empty) for each read (len n u8)
/// data: data (label/color) for each read (len n (u8))
#[derive(Ord, PartialOrd, Clone, PartialEq, Eq, Hash, Serialize, Deserialize, Debug)]
pub struct Reads<D> {
    storage: Vec<u64>,
    ends: Vec<usize>,
    exts: Vec<Exts>,
    data: Vec<D>,
    len: usize
}

impl<D: Clone + Copy> Reads<D> {

    pub fn new() -> Self {
        Reads {
            storage: Vec::new(),
            ends: Vec::new(),
            exts: Vec::new(),
            data: Vec::new(),
            len: 0
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn add_read(&mut self, read: DnaString, exts: Exts, data: D) {
        for base in read.iter() {
            self.push_base(base);
        }
        self.ends.push(self.len());
        self.data.push(data);
        self.exts.push(exts);
    }

    fn push_base(&mut self, base: u8) {
        let bit = (self.len() % 32) * 2;
        //println!("bit: {}, base: {}", bit, base);
        if bit != 0 {
            match self.storage.pop() {
                Some(last) => {
                    let last = last + ((base as u64) << (64 - bit - 2));
                    //println!("{:#066b}", last);
                    self.storage.push(last);
                },
                None => panic!("tried to push base to empty vector (?)")
            }
        } else {
            self.storage.push((base as u64) << 62);
        }
        self.len += 1; 
    }

    fn addr(&self, i: &usize) -> (usize, usize) {
        (i / 32, (i % 32 ) * 2)
    }

    /// get the `i`th read in a `Reads`
    pub fn get_read(&self, i: usize) -> (DnaString, Exts, D) {
        let mut sequence = DnaString::new();
        let end = self.ends[i];
        //let start = if i != 0 { self.ends[i-1] } else { 0 };
        let start = match i {
            0 => 0,
            1.. => self.ends[i-1]
        };

        for b in start..end {
            let (block, bit) = self.addr(&b);
            let base = ((self.storage[block] >> (62 - bit)) & 3u64) as u8;
            sequence.push(base);
        }

        (sequence, self.exts[i], self.data[i])
    }

}


#[cfg(test)]
mod tests {
    use crate::{dna_string::DnaString, Exts};
    use super::Reads;

    #[test]
    fn test_add() {

        let fastq = vec![
            (DnaString::from_acgt_bytes(str::as_bytes("ACGATCGT")), Exts::empty(), 6u8),
            (DnaString::from_acgt_bytes(str::as_bytes("GGGGGG")), Exts::empty(), 5u8),
            (DnaString::from_acgt_bytes(str::as_bytes("TTGGTT")), Exts::empty(), 7u8),
            (DnaString::from_acgt_bytes(str::as_bytes("ACCAC")), Exts::empty(), 8u8),
            (DnaString::from_acgt_bytes(str::as_bytes("TCCCT")), Exts::empty(), 9u8),
            (DnaString::from_acgt_bytes(str::as_bytes("ACCAC")), Exts::empty(), 8u8),
            (DnaString::from_acgt_bytes(str::as_bytes("TCCCT")), Exts::empty(), 9u8),

        ];

        let mut reads = Reads::new();
        for (read, ext, data) in fastq.clone() {
            reads.add_read(read, ext, data);
        }

        /* for no in reads.storage.iter() {
            println!("{:#b}", no)
        } */

         assert_eq!(reads.storage, vec![1791212948343256433, 5140577499666710528]);

        for i in 0..fastq.len() {
            //println!("read {}: {:?}", i, reads.get_read(i))
            assert_eq!(fastq[i], reads.get_read(i))
        }
    }
}
