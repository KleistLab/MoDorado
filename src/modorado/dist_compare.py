import pickle
import pysam
import pandas as pd 
import matplotlib.pyplot as plt
from scipy.stats import entropy
import argparse

def KL(args):
    with open(args.pickle, 'rb') as pfile:
        trna2mods = pickle.load(pfile)
    reference = pysam.FastaFile(args.reference)
    df_annotation = pd.read_excel(args.annotation)
    samples = args.samples.split(",") #["FH028", "FH017", "FH018", "FH019", "FH020", "FH029", "FH031", "elp6"]

    offset, min_cov = 23, args.cov #200
    output_name = args.output

    f = open(output_name, 'w')
    print("tRNA", "position", "nucleotide", sep="\t", end = '\t', file = f)
    for i in range(1, len(samples)):
        print(samples[i], end = '\t', file = f)
    print(file = f)

    for trna in reference.references:
        ref_pt = 0
        for index, row in df_annotation.iterrows():
            if row["tRNA"] == trna[5:-2]:
                if not pd.isna(row["nucleotide"]):
                    freqs = {}
                    print(trna[5:-2], row["position"], row["nucleotide"], sep="\t", end = "\t", file = f)
                    for sample_i in range(len(samples)):
                        G_flag = True # check if the position is a G
                        for mod in ['T+17802.', 'A+a.', 'C+m.']: #['A+17596.']
                            if ref_pt+offset in trna2mods[samples[sample_i]][trna][mod]:
                                G_flag = False
                                preds = trna2mods[samples[sample_i]][trna][mod][ref_pt+offset]
                                hist_data = plt.hist(preds, range = (11, 256), bins = 10)
                                freqs[sample_i] = hist_data[0] + 1 # pseudo-count for frequency 
                        if G_flag == True:
                            freqs[sample_i] = []
                        if sample_i != 0:
                            if sum(freqs[0]) > min_cov and sum(freqs[sample_i]) > min_cov: # check if min_cov is passed
                                kl = entropy(freqs[0], freqs[sample_i]) + entropy(freqs[sample_i], freqs[0])
                                # kl = distance.jensenshannon(freqs[0], freqs[sample_i])
                                print(kl, end="\t", file = f)
                            else:
                                kl = 0
                                print(kl, end="\t", file = f)
                        plt.close()
                    print(file = f)
                    ref_pt += 1
    f.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pickle", required=True, type=str, help="The input pickle object")
    parser.add_argument("-o", "--output", required=True, type=str, help="The output table of KL divergence over all positions between the two samples")
    parser.add_argument("-r", "--reference", required=True, type=str, help="The input reference fasta file")
    parser.add_argument("-a", "--annotation", required=True, type=str, help="The input annotation file")
    parser.add_argument("-s", "--samples", required=True, type=str, help="The list of samples, when more than 2 samples are added, the wildtype sample is assumed to be the first in the list")
    parser.add_argument("--cov", required=False, default = 200, type=int, help="The minimum coverage threshold for computing KL")
                        
    args = parser.parse_args()
    KL(args)