from modorado.utils import Sample, Reference, Position, Mod
import matplotlib.pyplot as plt
from scipy.stats import entropy
import argparse

def KL(args):
    sample1 = Sample(args.sample1, args.file1)
    sample2 = Sample(args.sample2, args.file2)
    min_cov = args.mincov
    outfile = args.output

    for sample in [sample1, sample2]:
        with open(sample.file) as f:
            next(f)
            for line in f:
                ref_str, mod_str, pos_str, scores_str = line.strip().split()
                ref = sample.references.setdefault(ref_str, Reference(ref_str))
                pos = ref.positions.setdefault(pos_str, Position(pos_str))
                mod = pos.mods.setdefault(mod_str, Mod(mod_str))
                mod.scores = list(map(int, scores_str.split(",")))

    with open(outfile, "w") as f:
        f.write(f"reference\tposition\tmod\tKL_divergence\n")
        for ref_name in sample1.references:
            ref = sample1.references[ref_name]
            for pos_name in ref.positions:
                pos = sample1.references[ref_name].positions[pos_name]
                for mod_name in pos.mods:
                    mod = sample1.references[ref_name].positions[pos_name].mods[mod_name]
                    scores1 = mod.scores
                    if ref_name in sample2.references and \
                        pos_name in sample2.references[ref_name].positions and \
                        mod_name in sample2.references[ref_name].positions[pos_name].mods:

                        scores2 = sample2.references[ref_name].positions[pos_name].mods[mod_name].scores
                        if len(scores1) > min_cov and len(scores2) > min_cov:
                            preds = [n for n in scores1 if n > 12] # filter out scores below 12 as it might distort visualisation
                            hist_data1 = plt.hist(preds, range = (11, 256), bins = 10)[0] + 1 
                            plt.close()

                            preds = [n for n in scores2 if n > 12] # 
                            hist_data2 = plt.hist(preds, range = (11, 256), bins = 10)[0] + 1 
                            plt.close()

                            kl = entropy(hist_data1, hist_data2) + entropy(hist_data2, hist_data1)

                            f.write(f"{ref_name}\t{pos_name}\t{mod_name}\t{kl}\n")

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