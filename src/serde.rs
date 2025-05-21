use std::hash::Hash;
use std::path::Path;
use std::{fs::File, io::{BufReader, BufWriter}};

use bimap::BiMap;
use boomphf::hashmap::BoomHashMap2;
use serde::{de::DeserializeOwned, Deserialize, Serialize};

use crate::graph::DebruijnGraph;
use crate::summarizer::SummaryConfig;
use crate::{reads::ReadsPaired, summarizer::ID, Exts};

/// serialize a [`ReadsPaired`] together with its hashed tags and IDs
#[derive(Debug, Serialize, Deserialize, PartialEq, Eq, Clone)]
pub struct SerReads<DI> {
    reads: ReadsPaired<DI>,
    hashed_tags: BiMap<String, u8>,
    hashed_ids: BiMap<String, ID>
}

impl<DI: Serialize + DeserializeOwned> SerReads<DI> {
    /// make a new [`SerReads`]
    pub fn new(reads: ReadsPaired<DI>, hashed_tags: BiMap<String, u8>, hashed_ids: BiMap<String, ID>) -> SerReads<DI> {
        SerReads {
            reads,
            hashed_tags,
            hashed_ids
        }
    }

    /// serialize a [`SerReads`]
    pub fn serialize<P: AsRef<Path>>(&self, path: P) {
        let file = File::create(path).expect("error creating file for serialized reads");
        let writer = BufWriter::new(file);
        bincode::serialize_into(writer, &self).expect("error serializing reads");
    }

    /// deserialize a [`SerReads`]
    pub fn deserialize_from<P: AsRef<Path>>(path: P) -> SerReads<DI> {
        let file = File::open(path).expect("error opening file with serialized reads");
        let reader = BufReader::new(file);

        match bincode::deserialize_from(reader) {
            Ok(ser_reads) => ser_reads,
            Err(err) => panic!("Error deserializing cached reads: {}\n Make sure the file was cached with a compatible version and parameters", err)
        }
    }

    /// get a reference of the underlying [`ReadsPaired<DI>`]
    pub fn reads(&self) -> &ReadsPaired<DI> {
        &self.reads
    }

    /// get a reference of the underlying hashed tags
    pub fn hashed_tags(&self) -> &BiMap<String, u8> {
        &self.hashed_tags
    }

    /// get a reference of the underlying hashed ids
    pub fn hashed_ids(&self) -> &BiMap<String, ID> {
        &self.hashed_ids
    }

    /// get the underlying read, hashed tags and hashed labels
    pub fn dissolve(self) -> (ReadsPaired<DI>, BiMap<String, u8>, BiMap<String, ID>) {
        (self.reads, self.hashed_tags, self.hashed_ids)
    }
}


/// serialize a [`SerKmers`] together with its hashed tags and IDs and a config
#[derive(Debug, Serialize, Deserialize)]
pub struct SerKmers<K: Hash, SD> {
    kmers: BoomHashMap2<K, Exts, SD>,
    hashed_tags: BiMap<String, u8>,
    hashed_ids: BiMap<String, ID>,
    config: SummaryConfig
}

impl<K, SD> SerKmers<K, SD> 
where 
    K: Serialize + DeserializeOwned + Hash,
    SD: Serialize + DeserializeOwned
{
    /// make a new [`SerKmers`]
    pub fn new(kmers: BoomHashMap2<K, Exts, SD>, hashed_tags: BiMap<String, u8>, hashed_ids: BiMap<String, ID>, config: SummaryConfig) -> SerKmers<K, SD> {
        SerKmers {
            kmers,
            hashed_tags,
            hashed_ids,
            config
        }
    }

    /// serialize a [`SerKmers`]
    pub fn serialize<P: AsRef<Path>>(&self, path: P) {
        let file = File::create(path).expect("error creating file for serialized k-mers");
        let writer = BufWriter::new(file);
        bincode::serialize_into(writer, &self).expect("error serializing k-mers");
    }

    /// deserialize a [`SerKmers`]
    pub fn deserialize_from<P: AsRef<Path>>(path: P) -> SerKmers<K, SD>  {
        let file = File::open(path).expect("error opening file with serialized k-mers");
        let reader = BufReader::new(file);

        match bincode::deserialize_from(reader) {
            Ok(ser_reads) => ser_reads,
            Err(err) => panic!("Error deserializing cached k-mers: {}\n Make sure the file was cached with a compatible version and parameters", err)
        }
    }

    /// get a reference of the underlying [`BoomHashMap2<K, Exts, SD>`]
    pub fn kmers(&self) -> &BoomHashMap2<K, Exts, SD> {
        &self.kmers
    }

    /// get a reference of the underlying hashed tags
    pub fn hashed_tags(&self) -> &BiMap<String, u8> {
        &self.hashed_tags
    }

    /// get a reference of the underlying hashed ids
    pub fn hashed_ids(&self) -> &BiMap<String, ID> {
        &self.hashed_ids
    }

    pub fn config(&self) -> &SummaryConfig {
        &self.config
    }

    /// get the underlying read, hashed tags, hashed labels and config
    pub fn dissolve(self) -> (BoomHashMap2<K, Exts, SD>, BiMap<String, u8>, BiMap<String, ID>, SummaryConfig) {
        (self.kmers, self.hashed_tags, self.hashed_ids, self.config)
    }
}

/// serialize a [`SerGraph`] together with its hashed tags and IDs and a config
#[derive(Debug, Serialize, Deserialize)]
pub struct SerGraph<K: Hash, SD> {
    graph: DebruijnGraph<K, SD>,
    hashed_tags: BiMap<String, u8>,
    hashed_ids: BiMap<String, ID>,
    config: SummaryConfig
}

impl<K, SD> SerGraph<K, SD> 
where 
    K: Serialize + DeserializeOwned + Hash,
    SD: Serialize + DeserializeOwned
{
    /// make a new [`SerGraph`]
    pub fn new(graph: DebruijnGraph<K, SD>, hashed_tags: BiMap<String, u8>, hashed_ids: BiMap<String, ID>, config: SummaryConfig) -> SerGraph<K, SD> {
        SerGraph {
            graph,
            hashed_tags,
            hashed_ids,
            config
        }
    }

    /// serialize a [`SerGraph`]
    pub fn serialize<P: AsRef<Path>>(&self, path: P) {
        let file = File::create(path).expect("error creating file for a serialized graph");
        let writer = BufWriter::new(file);
        bincode::serialize_into(writer, &self).expect("error serializing graph");
    }

    /// deserialize a [`SerGraph`]
    pub fn deserialize_from<P: AsRef<Path>>(path: P) -> SerGraph<K, SD>  {
        let file = File::open(path).expect("error opening file with a serialized graph");
        let reader = BufReader::new(file);

        match bincode::deserialize_from(reader) {
            Ok(ser_reads) => ser_reads,
            Err(err) => panic!("Error deserializing cached graph: {}\n Make sure the file was cached with a compatible version and parameters", err)
        }
    }

    /// get a reference of the underlying [`DebruijnGraph<K, SD>`]
    pub fn graph(&self) -> &DebruijnGraph<K, SD> {
        &self.graph
    }

    /// get a reference of the underlying hashed tags
    pub fn hashed_tags(&self) -> &BiMap<String, u8> {
        &self.hashed_tags
    }

    /// get a reference of the underlying hashed ids
    pub fn hashed_ids(&self) -> &BiMap<String, ID> {
        &self.hashed_ids
    }

    /// get a reference of the underlying [`SummaryConfig`]
    pub fn config(&self) -> &SummaryConfig {
        &self.config
    }

    /// get the underlying graph, hashed tags, hashed labels and config
    pub fn dissolve(self) -> (DebruijnGraph<K, SD>, BiMap<String, u8>, BiMap<String, ID>, SummaryConfig) {
        (self.graph, self.hashed_tags, self.hashed_ids, self.config)
    }
}

#[cfg(test)]
mod test {
    use std::fs::remove_file;

    use crate::summarizer::ID;

    use super::SerReads;

    #[test]
    fn test_ser_reads() {
        let ser_reads: SerReads<ID> = SerReads::deserialize_from("test_data/dlfksdkl");

        let cloned_ser_reads = ser_reads.clone();

        let all = ser_reads.dissolve();
        assert_eq!(cloned_ser_reads.reads(), &all.0);
        assert_eq!(cloned_ser_reads.hashed_tags(), &all.1);
        assert_eq!(cloned_ser_reads.hashed_ids(), &all.2);

        let ser_reads = SerReads::new(all.0, all.1, all.2);
        let ser_path = "test_data/new_ser_reads.reads.dbg";
        ser_reads.serialize(ser_path);

        let new_ser_reads: SerReads<ID> = SerReads::deserialize_from(ser_path);
        assert_eq!(ser_reads, new_ser_reads);
        remove_file(ser_path).unwrap();
    }
}

