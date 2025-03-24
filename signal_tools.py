import pysam
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np 
import pod5
import argparse
import pickle

def get_signal(dim, x, signal_calibrated, signal_cpts, stride):
    signal = np.zeros(dim)
    seq2ref = dict(x.get_aligned_pairs())
    
    for base_i in range(1, x.query_length): # skipping last base for now
        if seq2ref[base_i] != None: # the base is not an insertion
            ref_pos = seq2ref[base_i]
            chunk = signal_calibrated[signal_cpts[-base_i-1] : signal_cpts[-base_i]]
            step_size = len(chunk) // stride
            for step_i in range(stride):
                subchunk = chunk[step_i*step_size : (step_i+1)*step_size]
                signal[ref_pos * stride + step_i] = np.median(subchunk)
    return signal

def extract_signal(args):
    reference = pysam.FastaFile(args.ref)
    samples = args.samples.split(",")
    pod5s = args.pod5_dir.split(",")
    n_sub = int(args.subsample)

    ref2sigs = {}
    for i in range(len(samples)):
        sample = samples[i]
        pod5_file = pod5s[i]
        print(sample)
        samfile = pysam.AlignmentFile("data/" + sample + "_parasail_reference_filtered_fulllen.sam", "r")
        pod5file = pod5.DatasetReader(pod5_file) 

        ref2sigs[sample] = {}
        for ref in reference.references:
            ref2sigs[sample][ref] = []

        counter, iter = 0, samfile.fetch()
        for x in iter:
            if x.reference_name[0] != "S" and len(ref2sigs[sample][x.reference_name]) < n_sub: 
                counter += 1
                stride, move, ts = x.get_tag("mv")[0], np.array(x.get_tag("mv")[1:]), x.get_tag('ts')
                move_cpts = np.where(move == 1)[0]
                signal_cpts = ts + move_cpts * stride

                read_record, ref_seq = pod5file.get_read(x.query_name), reference.fetch(x.reference_name)
                
                if read_record != None: # no read splitting, original read_id exists
                    signal_calibrated = read_record.calibrate_signal_array(read_record.signal)
                    ref_signal = get_signal(stride * len(ref_seq), x, signal_calibrated, signal_cpts, stride)
                    
                
                else: # read has been split, new id applies
                    parent_id = x.get_tag("pi")
                    parent_record = pod5file.get_read(parent_id)
                    signal_calibrated = parent_record.calibrate_signal_array(parent_record.signal)
                    signal_cpts = signal_cpts + x.get_tag("sp")
                    ref_signal = get_signal(stride * len(ref_seq), x, signal_calibrated, signal_cpts, stride)
                
                matches = set(list(zip(*x.get_aligned_pairs(matches_only=True)))[1]) # returns the matched positions on the reference
                ref2sigs[sample][x.reference_name].append((ref_signal, matches))

        samfile.close()

    with open('output/signals_' + args.samples + '_' + args.subsample + '.pckl', 'wb') as handle:
        pickle.dump(ref2sigs, handle, protocol=pickle.HIGHEST_PROTOCOL)

def get_signal_median(sample, trna, ref2sigs):
    signal_list = np.zeros((len(ref2sigs[sample][trna]), len(ref2sigs[sample][trna][0][0])))
    for i in range(len(ref2sigs[sample][trna])):
        signal = ref2sigs[sample][trna][i][0]
        signal_list[i, :] = signal

    signal_list_masked = np.ma.masked_where(signal_list == 0, signal_list)
    signal_median = np.ma.median(signal_list_masked, axis=0).filled(0)
    return signal_median

def plot_signal(args, stride = 6):
    wt, mt = args.sample1, args.sample2
    print("WT: ", wt, "MT: ", mt)
    with open(args.signals, 'rb') as handle:
        ref2sigs = pickle.load(handle) 
    trna = args.trna
    kmer_len = args.kmer
    pos_of_interest = args.pos
    df_annotation = pd.read_excel("data/SI_table1.xlsx")

    df_trna = df_annotation[df_annotation["tRNA"] == trna[5:-2]]
    df_trna_cleared = df_trna.dropna(subset=["nucleotide"], ignore_index=True)
    row_pos = df_trna_cleared[df_trna_cleared["position"] == pos_of_interest].index[0]

    kmer_seq = ""
    for row in range(row_pos-kmer_len//2, row_pos+kmer_len//2 + 1):
        kmer_seq = kmer_seq + df_trna_cleared.iloc[row]["nucleotide"]
    print(kmer_seq)

    diff_median = get_signal_median(mt, trna, ref2sigs) - get_signal_median(wt, trna, ref2sigs)
    offset = np.median(diff_median)
    print("offset = ", offset)
    
    sns.set(rc = {'figure.figsize': (10, 6)})   
    sns.set_context("paper", font_scale = 2)
    sns.set_style("ticks")

    kmer_i = 23 + row_pos # 23 is the 5' adapter length in the reference
    length = kmer_len * stride
    x, y, hue = [], [], []
    for sample in [wt, mt]:
        for i in range(len(ref2sigs[sample][trna])):
            signal_chunk = ref2sigs[sample][trna][i][0][(kmer_i - kmer_len//2 - 2) * stride : (kmer_i + kmer_len//2 - 1) * stride]
            if 0 not in signal_chunk:
                if sample == wt:
                    y = y + list(signal_chunk + offset)
                else:
                    y = y + list(signal_chunk)       
                x = x + list(np.arange(0, length)) # signal positions to be plotted on the x-axis
                hue = hue + [sample] * length
    sns.lineplot(x = x, y = y, hue = hue, errorbar="sd", palette = ["#008080","#FF7F50"], lw=2, alpha=0.7)
    plt.vlines(x = np.arange(stride, length, stride), ymin = 0, ymax = 145, ls="--", color = "grey", lw=0.5)
    for base_i in range(len(kmer_seq)):
            plt.text(base_i * 6 +1, 125, kmer_seq[base_i])

    plt.ylabel("Signal intensity")
    plt.xticks([])
    plt.ylim(30, 110)
    plt.yticks(np.arange(30, 110+1, 20))
    sns.despine()
    plt.savefig("output/" + wt + "_" + mt + "_" + trna[5:-2] + "_" + str(kmer_len) + "mer.svg", format="svg", bbox_inches='tight')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    # To extract signals from a list of samples 
    extract_parser = subparsers.add_parser("extract", help="Extract signals from a list of samples")
    extract_parser.add_argument("--samples", type=str, required=True, help="A list of samples separated by comma, e.g. FH017,FH018,FH019")
    extract_parser.add_argument("--ref", type=str, required=True, help="The reference fasta file")
    extract_parser.add_argument("--pod5_dir", type=str, required=True, help="The locations of the pod5 files separated by comma, e.g. FH017.pod5,FH028.pod5")
    extract_parser.add_argument("--subsample", type=str, required=True, help="The (maximum) number of subsamples per tRNA (which is not necessarily reached by low coverage samples)")
    extract_parser.set_defaults(func=extract_signal)

    # To plot signals centering a particular position between two samples 
    plot_parser = subparsers.add_parser("plot", help="Plot signals centering a particular position between two samples")
    plot_parser.add_argument("--sample1", type=str, required=True, help="The wildtype sample to be plotted")
    plot_parser.add_argument("--sample2", type=str, required=True, help="The mutant sample to be plotted")
    plot_parser.add_argument("--signals", type=str, required=True, help="The signal pickle file")
    plot_parser.add_argument("--trna", type=str, required=True, help="The tRNA to be plotted, e.g. tRNA-Cys-GCA-1-1")
    plot_parser.add_argument("--pos", type=int, required=True, help="The position to be plotted")
    plot_parser.add_argument("--kmer", type=int, required=True, help="The kmer length to be plotted")
    plot_parser.set_defaults(func=plot_signal)

    args = parser.parse_args()
    args.func(args)