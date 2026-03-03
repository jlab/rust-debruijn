#!/bin/bash

# update dbg first

# dbg dir
DBG=$1

cd $DBG

cargo run -- -c ../marbel_datasets/sim_reads_100.csv -s sum --stranded -o ../rust-debruijn/test_data/marbel_100_sum --checkpoint -k 22
cargo run --features sample128 -- -c ../marbel_datasets/sim_reads_100.csv -s sum --stranded -o ../rust-debruijn/test_data/marbel_100_sum_128 --checkpoint -k 22
cargo run --features low-ks -- -c ../marbel_datasets/sim_reads_200.csv -s id-sum --checkpoint -o ../rust-debruijn/test_data/test_graph_ids -k 16
cargo run --features id4b --features low-ks -- -c ../marbel_datasets/sim_reads_200.csv -s id-sum --checkpoint -o ../rust-debruijn/test_data/test_graph_ids-4b -k16
cargo run --features low-ks -- -i data/test_400.fastq.gz -o ../rust-debruijn/test_data/400 -s tags-counts-sum -k 16 -t t
cargo run --features low-ks -- -c data/test2.csv -o ../rust-debruijn/test_data/sided -s tags-counts-p-em -k 16
cargo run --features sample128 --features low-ks -- -c data/test2.csv -o ../rust-debruijn/test_data/sided-128 -s tags-counts-p-em -k 16


rm ../rust-debruijn/test_data/marbel_100_sum_128.config.yaml
rm ../rust-debruijn/test_data/marbel_100_sum_128.graph.dbg
rm ../rust-debruijn/test_data/marbel_100_sum_128.reads.dbg
rm ../rust-debruijn/test_data/marbel_100_sum.config.yaml
rm ../rust-debruijn/test_data/marbel_100_sum.graph.dbg
rm ../rust-debruijn/test_data/marbel_100_sum.reads.dbg
rm ../rust-debruijn/test_data/test_graph_ids-4b.config.yaml
rm ../rust-debruijn/test_data/test_graph_ids.config.yaml
rm ../rust-debruijn/test_data/sided-128.config.yaml
rm ../rust-debruijn/test_data/sided.config.yaml
rm ../rust-debruijn/test_data/400.config.yaml
