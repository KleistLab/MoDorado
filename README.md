# MoDorado (off-label modification detection by repurposing Dorado models)
[![DOI](https://img.shields.io/badge/DOI-10.1093%2Fnar%2Fgkaf795-blue)](https://doi.org/10.1093/nar/gkaf795)

MoDorado is a light-weight algorithm that detects modification by the off-label use of pre-trained modification-specific models in nanopore direct RNA sequencing (**SQK-RNA004**). 

As of August 2025, the current version of the Dorado RNA model is rna004_130bps_sup@v5.2.0, supporting 8 modifications over A/U/C/G (inosine_m6A_2OmeA, pseU_2OmeU, m5C_2OmeC, 2OmeG). 

MoDorado currently supports analysis of tRNA modifications, but the pipeline can be adapted to other types of RNAs (see details below).

## Installation
To install MoDorado, you need Python >= 3.10 and ideally a dedicated virtual environment via conda/mamba to avoid conflicts.
```
conda create --name modo_env python=3.10
conda activate modo_env
```

To install, run
```
git clone https://github.com/KleistLab/MoDorado.git
cd MoDorado
pip install .
```

## Basecalling and modification calling with Dorado
To run [Dorado](https://github.com/nanoporetech/dorado), we select the **sup** model and all its modification models (Dorado should automatically download the latest version). 
```
dorado basecaller sup,inosine_m6A_2OmeA,m5C_2OmeC,pseU_2OmeU,2OmeG pod5s/ > basecalls.bam
```
`pod5s` is the directory of the raw pod5 files. 

## tRNA-specific analysis
Here, we focus on tRNA-specific analysis. For other types of RNAs, please refer to the section "General RNA analysis". 
### tRNA Alignment 
Popular aligners usually expect a fastq file as input, so we first need to convert `.sam` output from Dorado with samtools using 
```
samtools fastq -T "*" basecalls.sam > basecalls.fastq
```
For tRNA alignment, we recommend [Parasail](https://github.com/jeffdaily/parasail), based on evaluation performed in [Sun et al 2023](https://doi.org/10.1093/nar/gkad826). 
```
parasail_aligner -a sw_trace_striped_sse41_128_16 -M 2 -X 1 -c 10 -x -d  -O SAMH -t 6 -b 1000 -f tests/data/reference.fasta -q basecalls.fastq -g basecalls_parasail.sam
```
### tRNA read filtering (Parasail specific)
As Parasail performs all-versus-all pairwise alignment between every read and every reference, the alignments need to be filtered for the best candidate hit(s). Additionally, we need to specify the threshold for alignment score (AS) for filtering alignments, and the alignment start/length to obtain full length reads (in our dataset, there is a 23-nt adapter at the 5', thus `align_start` is set at 25).
```
modorado filter_parasail -i basecalls_parasail.sam -d basecalls.sam -o basecalls_parasail_filtered_fulllen.sam --AS 50 --align_start 25 --align_len 80
```
### tRNA parsing Dorado modification predictions 
This step is the same for tRNA or other types of RNA. If an alignment file with secondary/supplementary alignments is provided, only the primary alignment will be considered. The output is a tsv file with with all the dorado prediction scores aligned to a particular reference position, by each modification type.
```
modorado parse_dorado -r tests/data/reference.fasta -a basecalls_parasail_filtered_fulllen.sam -o dorado_preds.tsv
```
The alignment file can be either in sam or bam format.
### tRNA compare samples with KL divergence 
To be updated...

## General RNA analysis
The following has been tested on mRNA or viral transcriptomes.
### General RNA - Alignment with minimap2 or others
The current pipeline should work with all popular aligners for RNA reads. They usually require the reads to be in fastq format, so it would be necessary to convert the basecalled sam/bam file by 
```
samtools fastq -T "*" basecalls.sam > basecalls.fastq
```
The `-T "*"` option makes sure that all the optional fields in sam/bam are copied over to the fastq, which can then be aligned. For example, it has been tested with minimap2 with the following generic parameters (the `-y` option below copies the optional fields again from fastq to the output sam).
```
minimap2 -y -x map-ont -a --secondary=no reference.fasta basecalls.fastq > basecalls_minimap.sam 
```
### General RNA - parsing Dorado predictions
This step is exactly the same as the one for tRNA above. 
```
modorado parse_dorado -r reference.fasta -a basecalls_minimap.sam -o dorado_preds.tsv
```
### General RNA -  compare samples with KL divergence 
To be updated...
<!--
## 3. Distribution comparison with KL Divergence 
With Dorado results parsed, we can now compare two samples at each position of the tRNAs using the KL Divergence. To do this, we run 
```
modorado compare -p trna2mods.pckl -r reference.fasta -a tests/data/20241031_data_shifted_mods.xlsx -s sample1,sample2 --cov 100 -o kl_symmetric_mincov100.tsv
```
Here, the samples are again listed as a string separated by commas (when more than 2 samples are added, the wildtype sample is assumed to be the first in the list). The `100` at the end is a minimum coverage threshold for each tRNA.

This will generate a kl_symmetric_mincov100.tsv file, which contains the KL Divergence for each position of each tRNA. We provided a small example input to in the `data` folder (which results in 0 KLs when the minimum coverage is not fulfilled).
--->

## Miscellaneous: making signal plots
For making signal plots, include the options `--emit-moves --emit-sam` for Dorado basecalling.
### Signal extraction by subsampling reads from pod5 files
First, we subsample reads from the pod5 files by running
```
modorado extract_signal --sample FH017,FH028 -a tests/data/FH017_parasail_reference_filtered_fulllen.sam tests/data/FH028_parasail_reference_filtered_fulllen.sam --ref tests/data/reference.fasta --pod5_dir tests/data/ --subsample 200 -o tests/output/signals_FH017,FH028_200.pckl

```
Here, we need to specify the samples (comma separated, as many as needed, but the pod5 files should be listed in the same order as the sample list), and the location of their pod5 files. The subsample parameter is the number of subsampled reads.
This generates a pickle object in the output folder, containing the extracted signals for subsequent analysis or plotting.
### Plotting the signals of two samples (example Fig.5B in paper)
We can make signal plots comparing two samples by running the following. The example given is for Fig.5B of the manuscript on bioRxiv. 
```
modorado plot --sample1 FH028 --sample2 FH017 --signals tests/output/signals_FH017,FH028_200.pckl --trna tRNA-Cys-GCA-1-1 --pos 58 --kmer 11 --annotation tests/data/SI_table1.xlsx -o tests/output/FH028_FH017_Cys-GCA-1_11mer.svg
```
Here, we need to specify the signal file from the first step, the name of the tRNA, the position to be plotted and the length of the kmer centering the position. 
This should generate the desired plot in the output folder.
![plot](tests/data/FH028_FH017_Cys-GCA-1_11mer.svg)
<!--
## Quick start 
Here we show a quick toy example with two small samples in the tests folder, starting from parsing dorado model predictions (Step 2). The preprocessing steps require the original dorado basecalls, which exceed github's file size limits.
```
mkdir tests/output

# To compute KL divergence between two samples
modorado parse_dorado -r tests/data/reference.fasta -s FH028,FH017 -a tests/data/FH028_parasail_reference_filtered_fulllen.sam tests/data/FH017_parasail_reference_filtered_fulllen.sam -o tests/output/trna2mods.pckl
modorado compare -p tests/output/trna2mods.pckl -r tests/data/reference.fasta -a tests/data/20241031_data_shifted_mods.xlsx -s FH028,FH017 --cov 100 -o tests/output/kl_symmetric_mincov100_test.tsv

# To plot signals between two samples at a certain position  
modorado extract_signal --sample FH017,FH028 -a tests/data/FH017_parasail_reference_filtered_fulllen.sam tests/data/FH028_parasail_reference_filtered_fulllen.sam --ref tests/data/reference.fasta --pod5_dir tests/data/ --subsample 200 -o tests/output/signals_FH017,FH028_200.pckl
modorado plot --sample1 FH028 --sample2 FH017 --signals tests/output/signals_FH017,FH028_200.pckl --trna tRNA-Cys-GCA-1-1 --pos 58 --kmer 11 --annotation tests/data/SI_table1.xlsx -o tests/output/FH028_FH017_Cys-GCA-1_11mer.svg 
```
--->
