import argparse
from .signal_tools import extract_signal, plot_signal
from .filter_parasail import filter
from .parse_dorado import parse_dorado
from .dist_compare import KL


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    # filter parasail alignments based on alignment score and/or full length reads
    para_parser = subparsers.add_parser("filter_parasail", help="Filter parasail alignments based on alignment score (AS) and/or full length reads")
    para_parser.add_argument("-i", "--input", required=True, type=str, help="Input parasail alignment samfile")
    para_parser.add_argument("-o", "--output", required=True, type=str, help="Output filtered parasail alignment samfile")
    para_parser.add_argument("-d", "--dorado", required=True, type=str, help="Output filtered parasail alignment samfile")
    para_parser.add_argument("--AS", type=int, default=50, help="Minimum threshold of alignment score (AS) to filter parasail alignments")
    para_parser.add_argument("--align_start", type=int, default=25, help="Maximum threshold of alignment start position to filter parasail alignments")
    para_parser.add_argument("--align_len", type=int, default=80, help="Minimum threshold of alignment length to filter parasail alignments")
    para_parser.set_defaults(func=filter)

    dorado_parser = subparsers.add_parser("parse_dorado", help="Parse dorado MM and ML fields")
    dorado_parser.add_argument("-r", "--reference", required=True, type=str, help="The reference fasta file")
    dorado_parser.add_argument("-o", "--output", required=True, type=str, help="Output pickle object")
    dorado_parser.add_argument("-a", "--alignment", nargs='+', required=True, type=str, help="The input filtered alignment samfiles")
    dorado_parser.add_argument("-s", "--samples", required=True, type=str, help="The list of samples separated by comma, e.g. FH017,FH028,...")
    dorado_parser.set_defaults(func=parse_dorado)

    compare_parser = subparsers.add_parser("compare", help="Compute KL divergence between two samples")
    compare_parser.add_argument("-p", "--pickle", required=True, type=str, help="The input pickle object")
    compare_parser.add_argument("-o", "--output", required=True, type=str, help="The output table of KL divergence over all positions between the two samples")
    compare_parser.add_argument("-r", "--reference", required=True, type=str, help="The input reference fasta file")
    compare_parser.add_argument("-a", "--annotation", required=True, type=str, help="The input annotation file")
    compare_parser.add_argument("-s", "--samples", required=True, type=str, help="The list of samples, when more than 2 samples are added, the wildtype sample is assumed to be the first in the list")
    compare_parser.add_argument("--cov", required=False, default = 200, type=int, help="The minimum coverage threshold for computing KL")
    compare_parser.set_defaults(func=KL)                    

    # To extract signals from a list of samples 
    extract_parser = subparsers.add_parser("extract_signal", help="Extract signals from a list of samples")
    extract_parser.add_argument("--samples", type=str, required=True, help="A list of samples separated by comma, e.g. FH017,FH018,FH019")
    extract_parser.add_argument("-a", "--alignments", type=str, required=True, help="The alignment files corresponding to the samples")
    extract_parser.add_argument("--ref", type=str, required=True, help="The reference fasta file")
    extract_parser.add_argument("--pod5_dir", type=str, required=True, help="The locations of the pod5 files separated by comma, e.g. FH017.pod5,FH028.pod5")
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
    plot_parser.set_defaults(func=plot_signal)

    args = parser.parse_args()
    args.func(args)
    
if __name__ == "__main__":
    main()