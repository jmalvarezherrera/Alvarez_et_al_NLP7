#The code provide in these folders allows to analyze raw (fastq.qz) data for RNA-seq, ChIP-seq and DamID-seq.
rna_seq/ folder contains scripts for Preprocessing, Alignment, Get Read Counts and calculate differential gene expression with DESeq2.
We used the following versions of the programs:
fastx_toolkit (v0.0.14)
Tophat2 (v2.1.1)
genomicFeatures (v1.34.1)

chip_seq_and_damid_seq/aligments_peak_calling/ folder contains scripts for Preprocessing, Alignment, Peak calling for ChIP-seq and DamID-seq data.
We used the following versions of the programs:
fastx_toolkit (v0.0.14)
Bowtie2 (v2.3.2)
MACS2 (version 2.1.1)

chip_seq_and_damid_seq/chip_seq_profiles/ folder contains instructions to generate a heatmap and a profile of ChIP-seq data.
We used the following versions of the program:
DeepTools (v3.0.1) 

chip_seq_and_damid_seq/damid_profiles/ folder contains instructions to generate a profile of DamID-seq data.
We used the following versions of the program:
DeepTools (v3.0.1)
