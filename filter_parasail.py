import pysam
from copy import deepcopy
import sys

input_file = sys.argv[1]
output_file = input_file.split(".sam")[0] + "_filtered.sam"
dorado_file = input_file.split("_")[0] + ".sam"

samfile = pysam.AlignmentFile(input_file, "r")
outfile = pysam.AlignmentFile(output_file, 'w', template=samfile)
movefile = pysam.AlignmentFile(dorado_file, "r", check_sq=False)

iter = samfile.fetch()
read2best = {}
as_thres = 50

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
            
            outfile.write(read2best[x.query_name][i])

movefile.close()
outfile.close()
