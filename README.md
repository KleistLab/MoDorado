# MoDorado (off-label modification detection by repurposing Dorado models)
MoDorado is a light-weight algorithm that detects modification by off-label use of pre-trained modification-specific models in nanopore direct RNA sequencing (**SQK-RNA004**). The main features include the following:

- Directly utilizing the basecalling output from Dorado
- Supporting **A/U/C** bases (corresponding to the currently availalbe m6A/Ψ/m5C models for RNA004)
- **Off-label** predictions for modifications beyond the existing models by comparing prediction score distributions 

The current version of MoDorado supports analysis of tRNA modifications, with future extensions planned for mRNA and other types of RNAs.

## 1. Preprocessing 
### 1.1 Basecalling and modification calling with Dorado models
First, we perform basecalling and modification calling [Dorado](https://github.com/nanoporetech/dorado). Currently, the latest basecalling model for RNA004 is `rna004_130bps_sup@v5.1.0`, with the four models for m6A/Ψ/m5C/inosine selected. Additionally, for further processing (e.g. signal visualisation), the options `--emit-moves --emit-sam` should be used.


### 1.2 Alignment and read filtering (currently tRNA-specific)
To run [Parasail](https://github.com/jeffdaily/parasail), we first need to convert `.sam` output from Dorado by `samtools fastq -T "*" basecalls.sam > basecalls.fastq`.  
```
parasail_aligner -a sw_trace_striped_sse41_128_16 -M 2 -X 1 -c 10 -x -d  -O SAMH -t 6 -b 1000 -f data/reference.fasta -q basecalls.fastq -g basecalls_parasail.sam
```
As Parasail performs all-versus-all pairwise alignment between basecalled reads and each reference, the alignments need to be filtered for the best candidate hit(s). Additionally, we need to specify the threshold for alignment score (AS) for filtering alignments, and the alignment start/length for full length reads.
```
python filter_parasail.py -i basecalls_parasail.sam -o basecalls_parasail_filtered_fulllen.sam -d basecalls.sam --AS 50 --align_start 25 --align_len 80
```

## 2. Parsing Dorado model predictions
Dorado stores the modification information in `MM` and `ML` tags (for detailed description see [the SAM documentation](https://samtools.github.io/hts-specs/SAMv1.pdf)). To parse these into a data structure from each tRNA for each sequencing sample (e.g. FH017, FH028, etc), we run
```
python parse_dorado.py data/reference.fasta basecalls_parasail_filtered_fulllen.sam FH017,FH028
```
Here, all sequencing samples can be written as a long string separated by a comma, i.e. `sample1,sample2,sample3,...,samplen`. 

## 3. Distribution comparison with KL Divergence 
With Dorado results parsed, we can now compare two samples at each position of the tRNAs using the KL Divergence. To do this, we run 
```
python dist_compare.py output/trna2mods.pckl data/reference.fasta data/20241031_data_shifted_mods.xlsx FH028,FH017 100
```
Here, the samples are again listed as a string separated by commas (when more than 2 samples are added, the wildtype sample is assumed to be the first in the list). The `100` at the end is a minimum coverage threshold for each tRNA.

This will generate a kl_symmetric_mincov100.tsv file, which contains the KL Divergence for each position of each tRNA. We provided a small example input to in the `data` folder (which results in 0 KLs when the minimum coverage is not fulfilled).

## Add-on plotting functionalities
To reproduce the signal plots in the manuscript, we can do the following two steps. 

### Signal extraction by subsampling reads from pod5 files
First, we subsample reads from the pod5 files by running
```
python signal_tools.py extract --samples FH017,FH028 --ref data/reference.fasta --pod5_dir pod5/FH017.pod5,pod5/FH028.pod5  --subsample 200
```
Here, we need to specify the samples (comma separated, as many as needed, but the pod5 files should be listed in the same order as the sample list), and the location of their pod5 files. The subsample parameter is the number of subsampled reads.
This generates a pickle object in the output folder, containing the extracted signals for subsequent analysis or plotting.

### Plotting the signals of two samples (example Fig.5B in paper)
We can make signal plots comparing two samples by running the following. The example given is for Fig.5B of the manuscript on bioRxiv. 
```
python signal_tools.py plot --sample1 FH028 --sample2 FH017 --signals output/signals_FH017,FH028_200.pckl --trna tRNA-Cys-GCA-1-1 --pos 58 --kmer 11
```
Here, we need to specify the signal file from the first step, the name of the tRNA, the position to be plotted and the length of the kmer centering the position. 
This should generate the desired plot in the output folder.
![plot](./output/FH028_FH017_Cys-GCA-1_11mer.svg)

