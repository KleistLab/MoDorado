import pickle
import pysam
import sys 

reference = pysam.FastaFile(sys.argv[1])
align_dir = sys.argv[2]

#samples = ["FH017", "FH018", "FH019", "FH020", "FH028", "FH029", "FH031", "elp6", "SY001", "SY002", "SY003", "SY004"]
samples = sys.argv[3].split(",")

trna2mods, trna2cov = {}, {}
for sample in samples:
    print(sample)
    samfile = pysam.AlignmentFile(align_dir + "/" + sample + "_parasail_reference_filtered_fulllen.sam", "r")
    trna2mods[sample], trna2cov[sample] = dict(), dict()
    for ref in reference.references:
        trna2mods[sample][ref] = {}
        trna2cov[sample][ref] = 0
        ref_seq = reference.fetch(ref)
        base2pos = {
            "T+17802.": [pos for pos, char in enumerate(ref_seq) if char == "T"],
            "A+17596.": [pos for pos, char in enumerate(ref_seq) if char == "A"],
            "A+a.": [pos for pos, char in enumerate(ref_seq) if char == "A"],
            "C+m.": [pos for pos, char in enumerate(ref_seq) if char == "C"]
        }

        for mod in base2pos:
            trna2mods[sample][ref][mod] = {}
            for pos in base2pos[mod]:
                trna2mods[sample][ref][mod][pos] = []

    count = 0
    iter = samfile.fetch()
    for x in iter:
        trna2cov[sample][x.reference_name] += 1

        count += 1
        if count % 500000 == 0:
            print(count)
            
        index_dict, ML_dict = {}, {}
        for mod in x.get_tag("MM").split(";")[:-1]:
            items = mod.split(",")
            mod, mm = items[0], list(map(int, items[1:]))
            index, loc = [], -1
            for skip in mm:
                loc += skip + 1
                index.append(loc)
            index_dict[mod] = index 

        chunk = 0
        for mod in index_dict:
            ML_dict[mod] = x.get_tag("ML")[chunk : chunk + len(index_dict[mod])].tolist()
            chunk = chunk + len(index_dict[mod])

        seq2ref = dict(x.get_aligned_pairs()) # seq pos first, ref pos second
        base2pos_seq = {
            "T": [pos for pos, char in enumerate(x.seq) if char == "T"],
            "A": [pos for pos, char in enumerate(x.seq) if char == "A"],
            "C": [pos for pos, char in enumerate(x.seq) if char == "C"]
        }
        
        for mod in index_dict:
            for i, idx in enumerate(index_dict[mod]):
                base = mod[0] 
                seq_pos = base2pos_seq[base][idx]
                ref_pos = seq2ref[seq_pos]
                if ref_pos is not None:
                    ref_base = reference.fetch(x.reference_name, start = ref_pos, end=ref_pos+1)
                    if ref_base == base: # the ref base might be a mismatch, though quite rarely
                        trna2mods[sample][x.reference_name][mod][ref_pos].append(ML_dict[mod][i])

    samfile.close()


with open("./output/trna2mods.pckl", 'wb') as pfile:
    pickle.dump(trna2mods, pfile, protocol=pickle.HIGHEST_PROTOCOL)