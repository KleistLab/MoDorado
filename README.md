# MoDorado (off-label modification detection with Dorado models)
MoDorado is a light-weight algorithm that detects modification by off-label use of pre-trained modification-specific models in nanopore direct RNA sequencing. The main features include the following:

- Directly utilizing the basecalling output from Dorado
- Supporting **A/U/C** bases (corresponding to the currently availalbe m6A/Ψ/m5C models for RNA004)
- **Off-label** predictions for modifications beyond the existing models by comparing prediction score distributions 

The current version of MoDorado supports analysis of tRNA modifications, with future extensions planned for mRNA and other types of RNAs.

## Step 0: Basecalling with Dorado
To run MoDorado, it is necessary to have used the latest version Dorado for basecalling and modification calling. Currently for RNA004 data, this is the `rna004_130bps_sup@v5.1.0` model for basecalling, with the four models for m6A/Ψ/m5C/inosine selected. 

Additionally, for visualisation and easy processing, the options `--emit-moves --emit-sam` should be used.


## Step 1: Alignment with Parasail and read filtering
To run [Parasail](https://github.com/jeffdaily/parasail), we first need to convert `.sam` output from Dorado by `samtools fastq -T "*" basecalls.sam > basecalls.fastq`.  
```
parasail_aligner -a sw_trace_striped_sse41_128_16 -M 2 -X 1 -c 10 -x -d  -O SAMH -t 6 -b 1000 -f data/reference.fasta -q basecalls.fastq -g basecalls_parasail_reference.sam
```
As Parasail performs all-versus-all pairwise alignment between basecalled reads and each reference, we need to filter the alignments by finding the best hits
```
python filter_parasail.py basecalls_parasail_reference.sam
```
And filter full-length reads by
```
python filter_fulllen.py basecalls_parasail_reference_filtered.sam
```

## Step 2: Parse Dorado model predictions
Dorado stores the modification information in `MM` and `ML` tags (for detailed description see [the SAM documentation](https://samtools.github.io/hts-specs/SAMv1.pdf)). To parse these into a data structure from each tRNA and their nucleotides for each sequencing sample (e.g. FH017, FH028, etc), we run
```
python parse_dorado.py data/reference.fasta alignment FH017,FH028
```
Here, all sequencing samples can be written as a long string separated by a comma, i.e. `sample1,sample2,sample3,...,samplen`. 

## Step 3: Distribution comparison with KL Divergence 
With Dorado results parsed, we can now compare two samples at each position of the tRNAs using the KL Divergence. To do this, we run 
```
python dist_compare.py output/trna2mods.pckl data/reference.fasta data/20241031_data_shifted_mods.xlsx FH028,FH017 100
```
Here, the samples are again listed as a string separated by commas (when more than 2 samples are added, the wildtype sample is assumed to be the first in the list). The `100` at the end is a minimum coverage threshold for each tRNA.

This will generate a kl_symmetric_mincov100.tsv file, which contains the KL Divergence for each position of each tRNA. We provided a small example input to in the `data` folder (which results in 0 KLs when the minimum coverage is not fulfilled).
