import pysam
import argparse

class ReferenceFile:
    def __init__(self, filename: str):
        self.reference = pysam.FastaFile(filename)
        self.alphabet = ["A", "C", "G", "T"]
        
        self.base2pos = {ref: {base: [] for base in self.alphabet} for ref in self.reference.references}
        for ref in self.base2pos:
            ref_seq = self.reference.fetch(ref)
            for pos, char in enumerate(ref_seq):
                self.base2pos[ref][char.upper()].append(pos)
        
    def getReferences(self):
        return self.reference.references

class AlignedBasecalledFile:
    def __init__(self, filename: str, rf: ReferenceFile):
        self.filename = filename
        self.reference = rf
        self.alignment = pysam.AlignmentFile(filename)
        self.mods = {}
        # only read in the first line to get the mod types 
        iter = self.alignment.fetch()
        x = next(iter) 
        for mm in x.get_tag("MM").split(";")[:-1]:
            self.mods[mm.split(",")[0]] = None
        self.ref2mod = {ref: {mod: {} for mod in self.mods} for ref in self.reference.getReferences()}
    
    def parseMM(self):
        iter = self.alignment.fetch()
        for x in iter:
            # filter out non-primary alignments
            if x.is_supplementary or x.is_secondary or x.is_reverse:
                continue

            # to parse the MM field and find the (read) positions with a mod prediction 
            loc_dict = {}
            for mod in x.get_tag("MM").split(";")[:-1]:
                items = mod.split(",")
                mod, mm = items[0], list(map(int, items[1:]))
                locs, loc = [], -1
                for skip in mm:
                    loc += skip + 1 # for info on the MM tag, see Samfile Optional Fields Specification 
                    locs.append(loc)
                loc_dict[mod] = locs 

            # to find the corresponding ML scores based on the dimension of MM
            ML_dict, chunk = {}, 0
            for mod in self.mods:
                ML_dict[mod] = x.get_tag("ML")[chunk : chunk + len(loc_dict[mod])].tolist()
                chunk = chunk + len(loc_dict[mod])

            # to filter MM/ML based on alignment and only include correctly basecalled positions 
            seq2ref = dict(x.get_aligned_pairs(matches_only = True)) # seq pos first, ref pos second; filter out indels
            base2pos_seq = {base: [] for base in self.reference.alphabet}
            for pos, char in enumerate(x.query_sequence):
                base2pos_seq[char].append(pos)
            
            for mod in self.mods:
                for i, idx in enumerate(loc_dict[mod]):
                    base = mod[0]
                    seq_pos = base2pos_seq[base][idx]
                    if seq_pos in seq2ref:
                        ref_pos = seq2ref[seq_pos]
                        ref_base = self.reference.reference.fetch(x.reference_name, start = ref_pos, end=ref_pos+1)
                        if ref_base == base: # the ref base might be a mismatch, though quite rarely
                            if ref_pos not in self.ref2mod[x.reference_name][mod]:
                                self.ref2mod[x.reference_name][mod][ref_pos] = []
                            self.ref2mod[x.reference_name][mod][ref_pos].append(ML_dict[mod][i])

def parse_dorado(args):
    # samples = args.samples.split(",") if "," in args.samples else sample
    rf = ReferenceFile(args.reference)
    af = AlignedBasecalledFile(args.alignment, rf=rf)
    af.parseMM()
    with open(args.output, "w") as f:
        f.write(f"reference\tmodification\tposition\tscores\n")
        for ref in af.reference.getReferences():
            for mod in af.mods:
                for pos in sorted(af.ref2mod[ref][mod]): # sort the position numericallyq
                    preds = af.ref2mod[ref][mod][pos]
                    preds_str = ",".join(map(str, preds))
                    f.write(f"{ref}\t{mod}\t{pos}\t{preds_str}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--reference", required=True, type=str, help="The reference fasta file")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output pickle object")
    parser.add_argument("-a", "--alignment", nargs='+', required=True, type=str, help="The input filtered alignment samfiles")
    parser.add_argument("-s", "--samples", required=True, type=str, help="The list of samples separated by comma, e.g. FH017,FH028,...")
    
    args = parser.parse_args()
    parse_dorado(args)