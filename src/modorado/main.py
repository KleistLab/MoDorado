import argparse
from .signal_tools import extract_signal, plot_signal
from .filter_parasail import filter
from .parse_dorado import parse_dorado
from .dist_compare import KL
from .plot_trna import plot_trna


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    # filter parasail alignments based on alignment score and/or full length reads
    para_parser = subparsers.add_parser("filter_parasail", help="Filter parasail alignments based on alignment score (AS) and/or full length reads")
    para_parser.add_argument("-i", "--input", required=True, type=str, help="Input parasail alignment samfile")
    para_parser.add_argument("-o", "--output", required=True, type=str, help="Output filtered parasail alignment samfile")
    para_parser.add_argument("-d", "--dorado", required=True, type=str, help="Raw Dorado basecall output for addtional read information, either in sam or bam format")
    para_parser.add_argument("--AS", type=int, required=True, help="Minimum threshold of alignment score (AS) to filter parasail alignments")
    para_parser.add_argument("--align_start", type=int, required=True, help="Maximum threshold of alignment start position to filter parasail alignments")
    para_parser.add_argument("--align_len", type=int, required=True, help="Minimum threshold of alignment length to filter parasail alignments")
    para_parser.set_defaults(func=filter)

    dorado_parser = subparsers.add_parser("parse_dorado", help="Parse dorado MM and ML fields")
    dorado_parser.add_argument("-r", "--reference", required=True, type=str, help="The reference fasta file")
    dorado_parser.add_argument("-o", "--output", required=True, type=str, help="Output tsv with dorado prediction scores")
    dorado_parser.add_argument("-a", "--alignment", required=True, type=str, help="The input filtered alignment samfiles")
    # dorado_parser.add_argument("-a", "--alignment", nargs='+', required=True, type=str, help="The input filtered alignment samfiles")
    dorado_parser.set_defaults(func=parse_dorado)

    compare_parser = subparsers.add_parser("compare", help="Compute KL divergence between two samples")
    compare_parser.add_argument("-o", "--output", required=True, type=str, help="The output table of KL divergence over all positions between the two samples")
    compare_parser.add_argument("--sample1", required=True, type=str, help="Name of sample1")
    compare_parser.add_argument("--sample2", required=True, type=str, help="Name of sample2")
    compare_parser.add_argument("--file1", required=True, type=str, help="File location of sample1")
    compare_parser.add_argument("--file2", required=True, type=str, help="File location of sample1")
    compare_parser.add_argument("--mincov", required=False, default = 200, type=int, help="The minimum coverage threshold for computing KL")
    compare_parser.set_defaults(func=KL)           

    plot_trna_parser = subparsers.add_parser("plot_trna", help="Plot KL divergence values of tRNAs in a heatmap")
    plot_trna_parser.add_argument("-o", "--output", required=True, type=str, help="The heatmap visualisation of KL divergence over all tRNAs")
    plot_trna_parser.add_argument("--clip_5", required=True, type=int, help="The number of nt to be clipped at the 5' of each tRNA, so that the lengths match the multiple sequence alignment")
    plot_trna_parser.add_argument("--clip_3", required=True, type=int, help="The number of nt to be clipped at the 3' of each tRNA, so that the lengths match the multiple sequence alignment")
    plot_trna_parser.add_argument("-r", "--reference", required=True, type=str, help="The reference fasta file")
    plot_trna_parser.add_argument("--msa", required=True, type=str, help="The multiple sequence alignment file")
    plot_trna_parser.add_argument("--kl", required=True, type=str, help="The KL divergence output table")
    plot_trna_parser.add_argument("--mod", required=True, choices = ["m6A", "inosine", "m5C", "pseU", "Am", "Cm", "Gm", "Um"], help="The modification type in Dorado, must be one of the following: m6A, inosine, m5C, pseU, Am, Cm, Gm, Um")
    plot_trna_parser.set_defaults(func=plot_trna)               

    # To extract signals from a list of samples 
    extract_parser = subparsers.add_parser("extract_signal", help="Extract signals from a list of samples")
    extract_parser.add_argument("--samples", nargs='+', type=str, required=True, help="A list of samples separated by comma, e.g. FH017,FH018,FH019")
    extract_parser.add_argument("-a", "--alignment", nargs='+', type=str, required=True, help="The alignment files corresponding to the samples")
    extract_parser.add_argument("--ref", type=str, required=True, help="The reference fasta file")
    extract_parser.add_argument("--pod5s",  nargs='+', type=str, required=True, help="The locations of the pod5 files separated by comma, e.g. FH017.pod5,FH028.pod5")
    extract_parser.add_argument("--subsample", type=str, required=True, help="The (maximum) number of subsamples per tRNA (which is not necessarily reached by low coverage samples)")
    extract_parser.add_argument("-o", "--output", required=True, type=str, help="Output python pickle object storing signal data of subsampled tRNA reads")
    extract_parser.set_defaults(func=extract_signal)

    # To plot signals centering a particular position between two samples 
    plot_parser = subparsers.add_parser("plot", help="Plot signals centering a particular position between two samples")
    plot_parser.add_argument("-s1", "--sample1", type=str, required=True, help="The wildtype sample to be plotted")
    plot_parser.add_argument("-s2", "--sample2", type=str, required=True, help="The mutant sample to be plotted")
    plot_parser.add_argument("--signals", type=str, required=True, help="The signal pickle file")
    plot_parser.add_argument("--trna", type=str, required=True, help="The tRNA to be plotted, e.g. tRNA-Cys-GCA-1-1")
    plot_parser.add_argument("--pos", type=int, required=True, help="The position to be plotted")
    plot_parser.add_argument("--kmer", type=int, required=True, help="The kmer length to be plotted")
    plot_parser.add_argument("--annotation", type=str, required=True, help="The annotation file with nucleotide and their positions")
    plot_parser.add_argument("-o", "--output", type=str, required=True, help="The output plot location")
    plot_parser.add_argument("--ymax", type=int, required=False, default = 110, help="Optional max for the y-axis")
    plot_parser.add_argument("--offset", type=float, required=False, help="If set, offset is no longer computed from data")
    plot_parser.add_argument("--norm", type=float, required=False, default = 1, help="If set, normalise the signal") # to be improved
    plot_parser.set_defaults(func=plot_signal)

    args = parser.parse_args()
    args.func(args)
    
if __name__ == "__main__":
    main()