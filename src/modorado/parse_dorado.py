import argparse
from modorado.utils import ReferenceFile, AlignedBasecalledFile

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