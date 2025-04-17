import pickle
import pysam
import pandas as pd 
import matplotlib.pyplot as plt
from scipy.stats import entropy
import sys

with open(sys.argv[1], 'rb') as pfile:
    trna2mods = pickle.load(pfile)
reference = pysam.FastaFile(sys.argv[2])
df_annotation = pd.read_excel(sys.argv[3])
samples = sys.argv[4].split(",") #["FH028", "FH017", "FH018", "FH019", "FH020", "FH029", "FH031", "elp6"]

offset, min_cov = 23, int(sys.argv[5]) #200
output_name = "./output/kl_symmetric_mincov" + sys.argv[5] + ".tsv"

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