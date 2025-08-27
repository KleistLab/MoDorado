from modorado.utils import ReferenceFile
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def plot_trna(args):

    clip_5, clip_3 = args.clip_5, args.clip_3
    rf = ReferenceFile(args.reference)
    rf_msa = ReferenceFile(args.msa)
    kl_file = args.kl
    mod2symbol = {"m6A":"A+a.", "inosine":"A+17596.", "m5C":"C+m.", "pseU":"T+17802.", 
                  "Am":"A+69426.", "Cm":"C+19228.", "Gm":"G+19229.", "Um":"T+19227."}
    mod_to_plot = mod2symbol[args.mod]

    for trna in rf_msa.getReferences():
        rf_seq = rf.reference.fetch(trna)[clip_5:-clip_3]
        rf_msa_seq = rf_msa.reference.fetch(trna)
        if rf_seq != rf_msa_seq.replace("-", ""):
            raise Exception(f"{trna} reference and msa sequence do not match. Adjust clip_5 and clip_3? \n")

    trna2index = {trna: index for index, trna in enumerate(rf_msa.getReferences())}
    score_matrix = np.zeros((len(trna2index), len(rf_msa_seq)))
    print(f"The score matrix is of dimension {score_matrix.shape}.")

    n_row = 0
    trna2gaps = {trna: {} for trna in trna2index.keys()}
    for trna in trna2index:
        rf_msa_seq = rf_msa.reference.fetch(trna)
        del_indices = [i for i, c in enumerate(rf_msa_seq) if c == "-"]
        score_matrix[n_row, del_indices] = np.nan
        n_row += 1

        cum_gaps = 0
        for del_index in del_indices:
            cum_gaps += 1
            trna2gaps[trna][del_index - cum_gaps] = cum_gaps

    with open(kl_file) as f:
        next(f)
        for line in f:
            items = line.strip().split()
            trna, pos, mod, kl = items[0], int(items[1]), items[2], float(items[3])
            if trna not in trna2index:
                continue 

            if mod != mod_to_plot:# or trna != "tRNA-Arg-ACG-1-1":
                continue

            rf_seq = rf.reference.fetch(trna)[clip_5:-clip_3]
            rf_msa_seq = rf_msa.reference.fetch(trna)
            adjusted_pos = pos - clip_5
            
            offset = 0
            for pos_gap in trna2gaps[trna]:
                if pos_gap < adjusted_pos:
                    offset = trna2gaps[trna][pos_gap]
            adjusted_pos += offset
            
            if adjusted_pos >= 0 and adjusted_pos < len(rf_msa_seq):
                score_matrix[trna2index[trna], adjusted_pos] = kl

    sns.set(rc = {'figure.figsize': (20, 12)})   
    sns.set_context("paper", font_scale = 1)
    sns.set_style("ticks")

    colors = ["#FAFAD9", "#FD8D3C", "#800080", "#000000"]
    custom_cmap = mcolors.LinearSegmentedColormap.from_list(
        "custom", colors
    )

    # Set NaN values to grey
    cmap = custom_cmap.copy()
    cmap.set_bad(color="lightgrey")

    ax = sns.heatmap(score_matrix, 
                cmap=cmap, #"YlOrRd", 
                cbar=True, 
                cbar_kws={"shrink": 0.5},  
                yticklabels=[key[5:-2] for key in trna2index.keys()],
                xticklabels=False,
                edgecolors = 'white',
                linewidths=1,
                mask=np.isnan(score_matrix))

    cbar = ax.collections[0].colorbar

    cbar.ax.xaxis.set_label_position("top")   # put label on top
    cbar.ax.set_xlabel("KL_divergence", labelpad=10)  # overwrite default, adjust spacing

    num_cols = score_matrix.shape[1]
    ax.set_xticks(np.arange(4, num_cols, 5))          # positions (0-based, so 4=column 5)
    ax.set_xticklabels(np.arange(5, num_cols+1, 5))
    plt.xlabel("Adjusted Position")
    plt.ylabel("tRNA")
    plt.savefig(args.output, bbox_inches='tight')
    print(f"Plot saved in {args.output}.")

if __name__ == "__main__":
    plot_trna()