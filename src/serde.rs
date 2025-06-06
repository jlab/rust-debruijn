use std::fmt::Debug;
use std::hash::Hash;
use std::path::Path;
use std::{fs::File, io::{BufReader, BufWriter}};

use boomphf::hashmap::BoomHashMap2;
use serde::{de::DeserializeOwned, Deserialize, Serialize};

use crate::graph::DebruijnGraph;
use crate::reads::{ReadData, ReadDatas};
use crate::summarizer::{Summarizers, SummaryConfig, SummaryData, Translator};
use crate::Kmer;
use crate::{reads::ReadsPaired, Exts};

/// serialize a [`ReadsPaired`] together with its hashed tags and IDs
#[derive(Debug, Serialize, Deserialize, PartialEq, Eq, Clone)]
pub struct SerReads<DI> {
    reads: ReadsPaired<DI>,
    read_datas: ReadDatas,
    translator: Translator
}

impl<DI> SerReads<DI> {
    /// make a new [`SerReads`]
    pub fn new(reads: ReadsPaired<DI>, translator: Translator) -> SerReads<DI> 
    where DI: ReadData
    {
        SerReads {
            reads,
            read_datas: DI::read_datas(),
            translator
        }
    }

    /// serialize a [`SerReads`]
    pub fn serialize<P: AsRef<Path>>(&self, path: P) 
    where DI: Serialize
    {
        let file = File::create(path).expect("error creating file for serialized reads");
        let writer = BufWriter::new(file);
        bincode::serialize_into(writer, &self).expect("error serializing reads");
    }

    /// deserialize a [`SerReads`]
    pub fn deserialize_from<P: AsRef<Path> + Debug>(path: P) -> SerReads<DI> 
    where DI: DeserializeOwned
    {
        let file = File::open(&path).expect("error opening file with serialized reads");
        let reader = BufReader::new(file);

        match bincode::deserialize_from(reader) {
            Ok(ser_reads) => ser_reads,
            Err(err) => panic!("Error deserializing cached reads: {}\n Make sure the file was cached with a compatible version and parameters. \nFile: {:?}", err, path)
        }
    }

    /// get a reference of the underlying [`ReadsPaired<DI>`]
    pub fn reads(&self) -> &ReadsPaired<DI> {
        &self.reads
    }

    /// get a reference of the underlying translator (hashed tags & IDs)
    pub fn translator(&self) -> &Translator {
        &self.translator
    }

    /// get the parameters (read data kind)
    pub fn parameters(&self) -> ReadDatas {
        self.read_datas
    }

    /// get the underlying read, hashed tags and hashed labels
    pub fn dissolve(self) -> (ReadsPaired<DI>, Translator) {
        (self.reads, self.translator)
    }
}


/// serialize a [`SerKmers`] together with its hashed tags and IDs and a config
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct SerKmers<K: Hash, SD> {
    kmers: BoomHashMap2<K, Exts, SD>,
    k: usize,
    summarizer: Summarizers,
    translator: Translator,
    config: SummaryConfig
}

impl<K: Kmer, SD> SerKmers<K, SD> {
    /// make a new [`SerKmers`]
    pub fn new<DI>(kmers: BoomHashMap2<K, Exts, SD>, translator: Translator, config: SummaryConfig) -> SerKmers<K, SD> 
    where SD: SummaryData<DI>
    {
        SerKmers {
            kmers,
            summarizer: SD::summarizer(),
            k: K::k(),
            translator,
            config
        }
    }

    /// serialize a [`SerKmers`]
    pub fn serialize<P: AsRef<Path>>(&self, path: P) 
    where 
        K: Serialize,
        SD: Serialize
    {
        let file = File::create(path).expect("error creating file for serialized k-mers");
        let writer = BufWriter::new(file);
        bincode::serialize_into(writer, &self).expect("error serializing k-mers");
    }

    /// deserialize a [`SerKmers`]
    pub fn deserialize_from<P: AsRef<Path> + Debug>(path: P) -> SerKmers<K, SD> 
    where
        K: DeserializeOwned,
        SD: DeserializeOwned
    {
        let file = File::open(&path).expect("error opening file with serialized k-mers");
        let reader = BufReader::new(file);

        match bincode::deserialize_from(reader) {
            Ok(ser_reads) => ser_reads,
            Err(err) => panic!("Error deserializing cached k-mers: {}\n Make sure the file was cached with a compatible version and parameters. \nFile: {:?}", err, path)
        }
    }

    /// get a reference of the underlying [`BoomHashMap2<K, Exts, SD>`]
    pub fn kmers(&self) -> &BoomHashMap2<K, Exts, SD> {
        &self.kmers
    }

    /// get a reference of the underlying translator (hashed tags & IDs)
    pub fn translator(&self) -> &Translator {
        &self.translator
    }

    /// get a reference of the underlying config
    pub fn config(&self) -> &SummaryConfig {
        &self.config
    }

    /// get the set parameters (k, summarizer) of the `SerKmers`
    pub fn parameters(&self) -> (usize, Summarizers) {
        (self.k, self.summarizer)
    }

    /// get the underlying read, hashed tags, hashed labels and config
    pub fn dissolve(self) -> (BoomHashMap2<K, Exts, SD>, Translator, SummaryConfig) {
        (self.kmers, self.translator, self.config)
    }
}

/// serialize a [`SerGraph`] together with its hashed tags and IDs and a config
#[derive(Debug, Serialize, Deserialize)]
pub struct SerGraph<K: Hash, SD> {
    graph: DebruijnGraph<K, SD>,
    summarizer: Summarizers,
    k: usize,
    translator: Translator,
    config: SummaryConfig
}

impl<K: Kmer, SD> SerGraph<K, SD> {
    /// make a new [`SerGraph`]
    pub fn new<DI>(graph: DebruijnGraph<K, SD>, translator: Translator, config: SummaryConfig) -> SerGraph<K, SD> 
    where SD: SummaryData<DI>
    {
        SerGraph {
            graph,
            summarizer: SD::summarizer(),
            k: K::k(),
            translator,
            config
        }
    }

    /// serialize a [`SerGraph`]
    pub fn serialize<P: AsRef<Path>>(&self, path: P) 
    where 
        K: Serialize,
        SD: Serialize
    {
        let file = File::create(path).expect("error creating file for a serialized graph");
        let writer = BufWriter::new(file);
        bincode::serialize_into(writer, &self).expect("error serializing graph");
    }

    /// deserialize a [`SerGraph`]
    pub fn deserialize_from<P: AsRef<Path> + Debug>(path: P) -> SerGraph<K, SD>  
    where 
        K: DeserializeOwned,
        SD: DeserializeOwned
    {
        let file = File::open(&path).expect("error opening file with a serialized graph");
        let reader = BufReader::new(file);

        match bincode::deserialize_from(reader) {
            Ok(ser_reads) => ser_reads,
            Err(err) => panic!("Error deserializing cached graph: {}\n Make sure the file was cached with a compatible version and parameters. \n File: {:?}", err, path)
        }
    }

    /// get a reference of the underlying [`DebruijnGraph<K, SD>`]
    pub fn graph(&self) -> &DebruijnGraph<K, SD> {
        &self.graph
    }

    /// get a reference of the underlying translator (hashed tags and IDs)
    pub fn translator(&self) -> &Translator {
        &self.translator
    }
    /// get a reference of the underlying [`SummaryConfig`]
    pub fn config(&self) -> &SummaryConfig {
        &self.config
    }

    /// get the set parameters (k, summarizer) of the `SerGraph`
    pub fn parameters(&self) -> (usize, Summarizers) {
        (self.k, self.summarizer)
    }

    /// get the underlying graph, hashed tags, hashed labels and config
    pub fn dissolve(self) -> (DebruijnGraph<K, SD>, Translator, SummaryConfig) {
        (self.graph, self.translator, self.config)
    }
}

#[cfg(test)]
mod test {
    use std::fs::remove_file;

    use crate::{kmer::Kmer16, reads::ReadDatas, serde::{SerGraph, SerKmers}, summarizer::{IDSumData, Summarizers, ID}};

    use super::SerReads;

    // test files: cargo run -- -c ../marbel_datasets/sim_reads_200.csv -s id-sum --checkpoint -o ../rust-debruijn/test_data/test_graph_ids

    #[test]
    fn test_ser_reads() {
        let ser_reads: SerReads<ID> = SerReads::deserialize_from("test_data/test_graph_ids.reads.dbg");

        let cloned_ser_reads = ser_reads.clone();

        let all = ser_reads.dissolve();
        assert_eq!(cloned_ser_reads.reads(), &all.0);
        assert_eq!(cloned_ser_reads.translator(), &all.1);
        assert_eq!(cloned_ser_reads.parameters(), ReadDatas::ID);

        let ser_reads = SerReads::new(all.0, all.1);
        let ser_path = "test_data/new_ser_reads";
        ser_reads.serialize(ser_path);

        let new_ser_reads: SerReads<ID> = SerReads::deserialize_from(ser_path);
        assert_eq!(ser_reads, new_ser_reads);
        remove_file(ser_path).unwrap();
    }

    #[test]
    fn test_ser_kmers() {
        let ser_kmers: SerKmers<Kmer16, IDSumData> = SerKmers::deserialize_from("test_data/test_graph_ids.kmers.dbg");

        let cloned_ser_kmers = ser_kmers.clone();

        let all = ser_kmers.dissolve();
        assert_eq!(cloned_ser_kmers.kmers().len(), all.0.len());
        assert_eq!(cloned_ser_kmers.translator(), &all.1);
        assert_eq!(cloned_ser_kmers.config(), &all.2);
        assert_eq!(cloned_ser_kmers.parameters(), (16, Summarizers::IDSum));

        let ser_kmers = SerKmers::new(all.0, all.1, all.2);
        let ser_path = "test_data/new_ser_kmers";
        ser_kmers.serialize(ser_path);

        let new_ser_kmers: SerKmers<Kmer16, IDSumData> = SerKmers::deserialize_from(ser_path);
        assert_eq!(ser_kmers.kmers().len(), new_ser_kmers.kmers().len());
        assert_eq!(ser_kmers.translator(), new_ser_kmers.translator());
        assert_eq!(ser_kmers.config(), new_ser_kmers.config());

        remove_file(ser_path).unwrap();
    }

    #[test]
    fn test_ser_graph() {
        let ser_graph: SerGraph<Kmer16, IDSumData> = SerGraph::deserialize_from("test_data/test_graph_ids.graph.dbg");

        let all = ser_graph.dissolve();

        let ser_graph = SerGraph::new(all.0, all.1, all.2);
        let ser_path = "test_data/new_ser_graph";
        ser_graph.serialize(ser_path);

        let new_ser_kmers: SerGraph<Kmer16, IDSumData> = SerGraph::deserialize_from(ser_path);
        assert_eq!(ser_graph.graph().len(), new_ser_kmers.graph().len());
        assert_eq!(ser_graph.translator(), new_ser_kmers.translator());
        assert_eq!(ser_graph.config(), new_ser_kmers.config());
        assert_eq!(ser_graph.parameters(), (16, Summarizers::IDSum));

        remove_file(ser_path).unwrap();
    }
}

