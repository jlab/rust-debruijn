use std::io::BufRead;
use std::vec::IntoIter;

const BUF: usize = 64*1024;

pub struct FastqReader<R: BufRead> {
    buf_reader: R,
}

impl<R: BufRead> FastqReader<R> {
    pub fn new(inner: R) -> Self {
        FastqReader {
            buf_reader: inner
        }
    }

    pub fn sequences(mut self) -> Option<FastqSequenceIterator<R>> {
        // read into buffer
        let mut new_buffer = [0; BUF];
        let new_bytes = self.buf_reader.read(&mut new_buffer).expect("error reading bytes from fastq file");

        // return none if there is nothing in the file
        if new_bytes == 0 {return None}

        // turn buffer into iterator
        let iter_buffer = new_buffer[..new_bytes].to_vec().into_iter();

        Some(FastqSequenceIterator {
            fastq_reader: self,
            iter_buffer,
        })
    }
}

pub struct FastqSequenceIterator<R: BufRead> {
    fastq_reader: FastqReader<R>,
    iter_buffer: IntoIter<u8>,
}

impl<R: BufRead> Iterator for FastqSequenceIterator<R> {
    type Item = Vec<u8>;

    fn next(&mut self) -> Option<Self::Item> {
                
        loop {

            let mut is_rest = false;

            let mut rest= Vec::new();
            
            // TODO move while loops to function, is repetetive
            // return read seq, take iterator and rest as &mut
            let mut name = Vec::new();
            while !is_rest {
                match self.iter_buffer.next() {
                    Some(10) => { rest.push(10); break },
                    None => is_rest = true,
                    Some(byte) => { rest.push(byte); name.push(byte) }, 
                }
            }

            let mut sequence = Vec::new();
            while !is_rest {
                match self.iter_buffer.next() {
                    Some(10) => { rest.push(10); break },
                    None => is_rest = true,
                    Some(byte) => { rest.push(byte); sequence.push(byte) }, 
                }
            }
            
            let mut plus = Vec::new();
            while !is_rest {
                match self.iter_buffer.next() {
                    Some(10) => { rest.push(10); break },
                    None => is_rest = true,
                    Some(byte) => { rest.push(byte); plus.push(byte) }, 
                }
            }
            
            let mut quality = Vec::new();
            while !is_rest {
                match self.iter_buffer.next() {
                    Some(10) => { rest.push(10); break },
                    None => is_rest = true,
                    Some(byte) => { rest.push(byte); quality.push(byte) }, 
                }
            }


            if !is_rest {
                return Some(sequence);
            }
            // regular buffer is empty and needs to be filled, process needs to be repeated until a full record with a sequence is found
            else {
                // read new bytes and append them to the rest-buffer
                // return none if no new bytes were read
                let mut vec_buffer = rest;
                
                let mut new_buffer = [0; BUF];
                let new_bytes = self.fastq_reader.buf_reader.read(&mut new_buffer).expect("error reading bytes from fastq file");

                // TODO implement Error if there is still something in 
                if new_bytes == 0 {return None}

                vec_buffer.append(&mut new_buffer[..new_bytes].to_vec());
                self.iter_buffer = vec_buffer.into_iter();
            }

        }

    }


}

#[cfg(test)]
mod tests{
    use std::{fs::File, str};

    use super::FastqReader;

    #[test]
    fn test_read_sequence() {

        let path = "../dbg/data/test_100.fastq";
        let reader = std::io::BufReader::new(File::open(path).expect("error opening fastq file"));

        let fastq_reader = FastqReader::new(reader);
        
        match fastq_reader.sequences() {
            Some(seq_iter) => {
                for seq in seq_iter {
                    println!("sequence: {}", str::from_utf8(&seq).expect("?"))
                }
            },
            None => panic!("hä")
        }
    }
}

