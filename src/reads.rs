use std::ops::Range;
use itertools::Itertools;
use serde_derive::{Deserialize, Serialize};
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::{mem, str};
use crate::dna_string::DnaString;
use crate::{base_to_bits, base_to_bits_checked, Exts, Vmer};

#[derive(Debug, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Clone, Hash, Copy)]
pub enum Stranded {
    Forward,
    Reverse,
    Unstranded
}

/// Store many DNA sequences together with an Exts and data each compactly packed together
/// 
/// #### fields:
/// 
/// * `storage`: `Vec` with 2-bit encoded DNA bases of all sequences
/// * `ends`:  `Vec` with the ends (exclusive) of the separate sequences in the `Reads`
/// * `exts`: `Vec` with one Exts for each sequence
/// * `data`: `Vec` with data for each sequence
/// * `len`: length of all sequences together
/// * `stranded`: [`Stranded`] conveying the strandedness and direction of the reads
#[derive(Ord, PartialOrd, Clone, PartialEq, Eq, Hash, Serialize, Deserialize, Debug)]
pub struct Reads<D> {
    storage: Vec<u64>,
    ends: Vec<usize>,
    exts: Vec<Exts>,
    data: Vec<D>,
    len: usize,
    stranded: Stranded
}

impl<D: Clone + Copy> Reads<D> {

    /// Returns a new `Reads`
    pub fn new(stranded: Stranded) -> Self {
        Reads {
            storage: Vec::new(),
            ends: Vec::new(),
            exts: Vec::new(),
            data: Vec::new(),
            len: 0,
            stranded
        }
    }

    #[inline(always)]
    pub fn stranded(&self) -> Stranded {
        self.stranded
    }

    #[inline(always)]
    /// Returns the number of reads stored
    pub fn n_reads(&self) -> usize {
        self.ends.len()
    }

    pub fn mem(&self) -> usize {
        mem::size_of_val(self) + size_of_val(&*self.storage) + size_of_val(&*self.data) + size_of_val(&*self.ends) + size_of_val(&*self.exts)
    }

    /// Adds a new read to the `Reads`
    // maybe push_base until u64 is full and then do extend like in DnaString::extend ? with accellerated mode
    pub fn add_read<V: Vmer>(&mut self, read: V, exts: Exts, data: D) {
        for base in read.iter() {
            self.push_base(base);
        }
        self.ends.push(self.len);
        self.data.push(data);
        self.exts.push(exts);
       
    }

    /// Transforms a `[(vmer, exts, data)]` into a `Reads` - watch for memory usage
    // TODO test if memory efficient
    pub fn from_vmer_vec<V: Vmer, S: IntoIterator<Item=(V, Exts, D)>>(vec_iter: S, stranded: Stranded) -> Self {
        let mut reads = Reads::new(stranded);
        for (vmer, exts, data) in vec_iter {
            for base in vmer.iter() {
                reads.push_base(base);
            }
            reads.ends.push(reads.len);
            reads.data.push(data);
            reads.exts.push(exts);
        }

        reads.storage.shrink_to_fit();
        reads.data.shrink_to_fit();
        reads.exts.shrink_to_fit();
        reads.ends.shrink_to_fit();
        
        reads
    }


    /// add ASCII encoded bases to the Reads
    /// 
    /// will transform all ascii characters outside of ACGTacgt into A
    /// see also: [`Reads::add_from_bytes_checked`]
    pub fn add_from_bytes(&mut self, bytes: &[u8], exts: Exts, data: D) {
        
        // fill the last incomplete u64 block
        let missing = 32 - (self.len % 32);
        if missing != 0 {
            if  missing > bytes.len() {
                let fill = bytes.iter().map(|c| base_to_bits(*c));
                self.extend(fill);
                self.exts.push(exts);
                self.data.push(data);
                self.ends.push(self.len);
                return;
            } else {
                let fill = bytes[0..missing].iter().map(|c| base_to_bits(*c));
                self.extend(fill);
            }
        }
        
        // Accelerated avx2 mode. Should run on most machines made since 2013.
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                for chunk in bytes[missing..bytes.len()].chunks(32) {
                    if chunk.len() == 32 {
                        let (conv_chunk, _) = unsafe { crate::bitops_avx2::convert_bases(chunk) };
                        let packed = unsafe { crate::bitops_avx2::pack_32_bases(conv_chunk) };
                        self.storage.push(packed);
                        self.len += 32;
                    } else {
                        let b = chunk.iter().map(|c| base_to_bits(*c));
                        self.extend(b);
                    }
                }

                self.exts.push(exts);
                self.data.push(data);
                self.ends.push(self.len);
                return;
            }
        }

        let b = bytes.iter().map(|c| base_to_bits(*c));
        self.extend(b);

        self.exts.push(exts);
        self.data.push(data);
        self.ends.push(self.len);
        
    }

    /// add ASCII encoded bases to the Reads
    /// 
    /// will return `false` if the bytes contained characters outside of `ACGTacgt`, otherwise return true and add the bases
    /// see also: [`Reads::add_from_bytes`]
    pub fn add_from_bytes_checked(&mut self, bytes: &[u8], exts: Exts, data: D) -> bool {

        let (_, corrects): (Vec<u8>, Vec<bool>) = bytes.iter().map(|c| base_to_bits_checked(*c)).collect();
        if corrects.iter().contains(&false) { return false }

        
        // fill the last incomplete u64 block
        let missing = 32 - (self.len % 32);
        if missing != 0 {
            if  missing > bytes.len() {
                let fill = bytes.iter().map(|c| base_to_bits(*c));
                self.extend(fill);
                self.exts.push(exts);
                self.data.push(data);
                self.ends.push(self.len);
                return true;
            } else {
                let fill = bytes[0..missing].iter().map(|c| base_to_bits(*c));
                self.extend(fill);
            }
        }
        
        // Accelerated avx2 mode. Should run on most machines made since 2013.
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                for chunk in bytes[missing..bytes.len()].chunks(32) {
                    if chunk.len() == 32 {
                        let (conv_chunk, _) = unsafe { crate::bitops_avx2::convert_bases(chunk) };
                        let packed = unsafe { crate::bitops_avx2::pack_32_bases(conv_chunk) };
                        self.storage.push(packed);
                        self.len += 32;
                    } else {
                        let b = chunk.iter().map(|c| base_to_bits(*c));
                        self.extend(b);
                    }
                }

                self.exts.push(exts);
                self.data.push(data);
                self.ends.push(self.len);
                return true;
            }
        }

        let b = bytes.iter().map(|c| base_to_bits(*c));
        self.extend(b);

        self.exts.push(exts);
        self.data.push(data);
        self.ends.push(self.len);

        true         
    }


    /// Add new 2-bit encoded base to the `Reads`
    fn push_base(&mut self, base: u8) {
        let bit = (self.len % 32) * 2;
        if bit != 0 {
            match self.storage.pop() {
                Some(last) => {
                    let last = last + ((base as u64) << (64 - bit - 2));
                    self.storage.push(last);
                },
                None => panic!("tried to push base to empty vector (?)")
            }
        } else {
            self.storage.push((base as u64) << 62);
        }
        self.len += 1; 
    }


    /// extend the reads' storage by 2-bit encoded bases
    fn extend(&mut self, mut bytes: impl Iterator<Item = u8>) {
        // fill the last incomplete u64 block
        while self.len % 32 != 0 {
            match bytes.next() {
                Some(b) => self.push_base(b),
                None => return,
            }
        }

        let mut bytes = bytes.peekable();

        // chunk the remaining items into groups of at most 32 and handle them together
        while bytes.peek().is_some() {
            let mut val: u64 = 0;
            let mut offset = 62;
            let mut n_added = 0;

            for _ in 0..32 {
                if let Some(b) = bytes.next() {
                    assert!(b < 4);
                    val |= (b as u64) << offset;
                    offset -= 2;
                    n_added += 1;
                } else {
                    break;
                }
            }

            self.storage.push(val);
            self.len += n_added;
        }
    }

    #[inline(always)]
    fn addr(&self, i: &usize) -> (usize, usize) {
        (i / 32, (i % 32 ) * 2)
    }

    /// get the `i`th read in a `Reads`
    pub fn get_read(&self, i: usize) -> Option<(DnaString, Exts, D)> {
        if i >= self.n_reads() { return None }

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

        Some((sequence, self.exts[i], self.data[i]))
    }


    /// shrink the vectors' capacity to fit the length
    /// 
    /// use sparsely
    pub fn shrink_to_fit(&mut self)  {
        self.storage.shrink_to_fit();
        self.data.shrink_to_fit();
        self.exts.shrink_to_fit();
        self.ends.shrink_to_fit();
    }



    /// Iterate over the reads as (DnaString, Exts, D).
    pub fn iter(&self) -> ReadsIter<'_, D> {
        ReadsIter {
            reads: self,
            i: 0,
            end: self.n_reads(),
            length: self.n_reads(),
        }
    }

    /// Iterate over a range start reads as (DnaString, Exts, D).
    pub fn partial_iter(&self, range: Range<usize>) -> ReadsIter<'_, D> {
        assert!(range.end <= self.n_reads());
        assert!(range.start < self.n_reads());
        assert!(range.start < range.end);
        ReadsIter {
            reads: self,
            i: range.start,
            end: range.end,
            length: (range.end - range.start)
        }
    }

    pub fn print_bin(&self) {
        for element in self.storage.iter() {
            print!("{:#066b} - ", element)
        }
        println!()
    }

}

impl<D: Clone + Copy> Default for Reads<D> {
    fn default() -> Self {
        Self::new(Stranded::Unstranded)
    }
}




/// Iterator over values of a DnaStringoded sequence (values will be unpacked into bytes).
pub struct ReadsIter<'a, D> {
    reads: &'a Reads<D>,
    i: usize,
    end: usize,
    length: usize,
}

impl<D: Clone + Copy> Iterator for ReadsIter<'_, D> {
    type Item = (DnaString, Exts, D);

    fn next(&mut self) -> Option<Self::Item> {
        if (self.i < self.reads.n_reads()) && (self.i < self.end) {
            let value = self.reads.get_read(self.i);
            self.i += 1;
            value
        } else {
            None
        }
    }
}

impl<D: Copy> ExactSizeIterator for ReadsIter<'_, D> {
    fn len(&self) -> usize {
        self.length
    }
}

impl<D: Clone + Copy + Debug> Display for Reads<D> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let vec: Vec<(DnaString, Exts, D)> = self.iter().collect();
        write!(f, "{:?}", vec)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ReadsPaired<D> {
    Empty,
    Unpaired { reads: Reads<D> },
    Paired { r1: Reads<D>, r2: Reads<D> },
    Combined {r1: Reads<D>, r2: Reads<D>, unpaired: Reads<D>}
}

impl<D: Clone + Copy> ReadsPaired<D> {
    pub fn iterable(&self) -> Vec<&Reads<D>> {
        match self {
            Self::Empty => vec![],
            Self::Unpaired { reads  } => vec![reads],
            Self::Paired { r1, r2 } => vec![r1, r2],
            Self::Combined { r1, r2, unpaired } => vec![r1, r2, unpaired],
        }
    }

    pub fn n_reads(&self) -> usize {
        match self {
            Self::Empty => 0,
            Self::Unpaired { reads  } => reads.n_reads(),
            Self::Paired { r1, r2 } => r1.n_reads() + r2.n_reads(),
            Self::Combined { r1, r2, unpaired } => r1.n_reads() + r2.n_reads() + unpaired.n_reads(),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::time;

    use itertools::enumerate;
    use rand::random;

    use crate::{dna_string::DnaString, reads::Stranded, Exts};
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


        let mut reads = Reads::new(Stranded::Unstranded);
        for (read, ext, data) in fastq.clone() {
            reads.add_read(read, ext, data);
        }

        println!("reads: {:#?}", reads);


        /* for no in reads.storage.iter() {
            println!("{:#b}", no)
        } */

         assert_eq!(reads.storage, vec![1791212948343256433, 5140577499666710528]);

        for (i, _) in fastq.iter().enumerate() {
            //println!("read {}: {:?}", i, reads.get_read(i))
            assert_eq!(fastq[i], reads.get_read(i).unwrap())
        }

        for (seq, _, _) in reads.iter() {
            println!("{:?}, {}", seq, seq.len())
        }
        println!();

        for read in reads.partial_iter(5..7) {
            println!("{:?}", read)
        }

        println!("memory usage: {}", reads.mem())
    }

    #[test]
    fn test_get_read() {
        let mut reads = Reads::new(Stranded::Unstranded);
        //reads.add_read(DnaString::from_acgt_bytes("AGCTAGCTAGC".as_bytes()), Exts::empty(), 67u8);
        reads.add_from_bytes("ACGATCGNATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACTGATCGATCG".as_bytes(), Exts::empty(), 67u8);
        let read = reads.get_read(0);
        println!("{:?}", read);
        println!("{:#066b}", reads.storage[0]);
        println!("{:?}", reads)
    }

    #[test]
    fn test_add_from_bytes() {
        let dna = [
            "AAGCGGAGATTATTCACGAGCATCGCGTAC".as_bytes(),
            "GATCGATGCATGCTAGA".as_bytes(),
            "ACGTAAAAAAAAAATTATATAACGTACGTAAAAAAAAAATTATATAACGTAACGTAAAAAAAAAAATTATAATAACGT".as_bytes(),
            "AGCTAGCTAGCTGACTGAGCGACTGA".as_bytes(),
            "AGCTAGCTAGCTGACTGAGCGACTGACGGATC".as_bytes(),
            "TTTTTTTTTTTTTTTTTTTTTTTT".as_bytes(),
            "ACGATCGAATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACTGATCGATCG".as_bytes(),
            "ACGATCGATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACTGATCGATCGAAGGGCAGTTAGGCCGTAAGCGCGAT".as_bytes(),
        ];

        let mut reads: Reads<u8> = Reads::new(Stranded::Unstranded);
        for seq in dna {
            reads.add_from_bytes(seq, Exts::empty(), random());

        }

        for (i, seq) in enumerate(reads.iter()) {
            let sequence = DnaString::from_acgt_bytes(dna[i]);
            assert_eq!(seq.0, sequence);
        }
    }

    #[test]
    fn test_add_from_bytes_checked() {
        let dna = [
            "AAGCGGAGATTATTCACGAGCATCGCGTAC".as_bytes(),
            "GATCGATGCATGCTAGA".as_bytes(),
            "ACGTAAAAAAAAAATTATATAACGTACGTAAAAAAAAAANTTATATAACGTAACGTAAAAAAAAAAATTATAATAACGT".as_bytes(),
            "AGCTAGCTAGCTGACNGAGCGACTGA".as_bytes(),
            "AGCTAGCTAGCTGACTGAGCGACTGACGGATC".as_bytes(),
            "TTTTTTTTTTTTTTTTTTTTTTTT".as_bytes(),
            "ACGATCGAATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACNNNTGATCGATCG".as_bytes(),
            "ACGATCGATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACTGATCGATCGAAGGGCAGTTAGGCCGTAAGCGCGAT".as_bytes(),
            "A".as_bytes(),
            "AAAAN".as_bytes(),
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN".as_bytes(),
        ];

        let mut reads: Reads<u8> = Reads::new(Stranded::Unstranded);
        let mut corrects = Vec::new();
        for seq in dna {
            corrects.push(reads.add_from_bytes_checked(seq, Exts::empty(), random()));
        }

        let mut read_counter = 0;

        for (i, correct) in enumerate(corrects) {
            let sequence = DnaString::from_acgt_bytes_checked(dna[i]);
            match correct {
                true => {
                    let read = reads.get_read(read_counter).unwrap();
                    assert_eq!(sequence.unwrap(), read.0);
                    read_counter += 1;

                },
                false => assert!(sequence.is_err()),
            }
        }
    }

    #[test]
    fn test_speed_from_bytes() {
        let dnas = [
            "AAGCGGAGATTATTCACGAGCATCGCGTAC".as_bytes(),
            "GATCGATGCATGCTAGA".as_bytes(),
            "ACGTAAAAAAAAAATTATATAACGTACGTAAAAAAAAAANTTATATAACGTAACGTAAAAAAAAAAATTATAATAACGT".as_bytes(),
            "AGCTAGCTAGCTGACNGAGCGACTGA".as_bytes(),
            "AGCTAGCTAGCTGACTGAGCGACTGACGGATC".as_bytes(),
            "TTTTTTTTTTTTTTTTTTTTTTTT".as_bytes(),
            "ACGATCGAATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACNNNTGATCGATCG".as_bytes(),
            "ACGATCGATGCTAGCTGATCGGCGACGATCGATGCTAGCTGATCGTAGCTGACTGATCGATCGAAGGGCAGTTAGGCCGTAAGCGCGAT".as_bytes(),
            "A".as_bytes(),
            "AAAAN".as_bytes(),
            "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN".as_bytes(),
        ];

        const REPS: usize = 500000;
        /*
        test with 5 ml:
            through DnaString: 88.09323 s 
            direct to Read: 44.41913 s
         */


        let ds_start= time::Instant::now();
        let mut reads: Reads<u8> = Reads::new(Stranded::Unstranded);
        for _i in 0..REPS {
            for dna in dnas {
                reads.add_read(DnaString::from_acgt_bytes(dna), Exts::empty(), random());
            }
        }
        let ds_finish = ds_start.elapsed();

        let r_start= time::Instant::now();
        let mut reads: Reads<u8> = Reads::new(Stranded::Unstranded);
        for _i in 0..REPS {
            for dna in dnas {
                reads.add_from_bytes(dna, Exts::empty(), random());
            }
        }
        let r_finish = r_start.elapsed();

        println!("through DnaString: {} s \n direct to Read: {} s", ds_finish.as_secs_f32(), r_finish.as_secs_f32())


    }
}
