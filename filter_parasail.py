import pysam
import argparse
from copy import deepcopy

def filter(input_file, output_file, dorado_file, as_thres, align_start, align_len):
    samfile = pysam.AlignmentFile(input_file, "r")
    outfile = pysam.AlignmentFile(output_file, 'w', template=samfile)
    movefile = pysam.AlignmentFile(dorado_file, "r", check_sq=False)
    
    iter = samfile.fetch()
    read2best = {}
    for x in iter:
        if x.get_tag("AS") < as_thres:
            continue

        if x.query_name not in read2best:
            read2best[x.query_name] = []
            read2best[x.query_name].append(deepcopy(x))
            continue
            
        if x.get_tag("AS") > read2best[x.query_name][0].get_tag("AS"):
            read2best[x.query_name] = []
            read2best[x.query_name].append(deepcopy(x))
            continue
        
        if x.get_tag("AS") == read2best[x.query_name][0].get_tag("AS"):
            read2best[x.query_name].append(deepcopy(x))
    samfile.close()

    iter = movefile.fetch()
    for x in iter:
        if x.query_name in read2best:
            for i in range(len(read2best[x.query_name])):
                read2best[x.query_name][i].set_tag("ts", x.get_tag("ts"))
                read2best[x.query_name][i].set_tag("mv", x.get_tag("mv"))
                if x.has_tag("MM"):
                    read2best[x.query_name][i].set_tag("MM", x.get_tag("MM"))
                    read2best[x.query_name][i].set_tag("ML", x.get_tag("ML"))

                if x.has_tag("pi"):
                    read2best[x.query_name][i].set_tag("pi", x.get_tag("pi"))
                    read2best[x.query_name][i].set_tag("sp", x.get_tag("sp"))
                
                if read2best[x.query_name][i].reference_start < align_start and read2best[x.query_name][i].reference_length > align_len: # Filter full length reads
                    outfile.write(read2best[x.query_name][i])

    movefile.close()
    outfile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=str, help="Input parasail alignment samfile")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output filtered parasail alignment samfile")
    parser.add_argument("-d", "--dorado", required=True, type=str, help="Output filtered parasail alignment samfile")
    parser.add_argument("--AS", type=int, default=50, help="Minimum threshold of alignment score (AS) to filter parasail alignments")
    parser.add_argument("--align_start", type=int, default=25, help="Maximum threshold of alignment start position to filter parasail alignments")
    parser.add_argument("--align_len", type=int, default=80, help="Minimum threshold of alignment length to filter parasail alignments")

    args = parser.parse_args()
    filter(args.input, args.output, args.dorado, args.AS, args.align_start, args.align_len)