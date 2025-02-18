import pysam
import sys

input_file = sys.argv[1]
output_file = input_file.split(".sam")[0] + "_fulllen.sam"

samfile = pysam.AlignmentFile(input_file, "r")
outfile = pysam.AlignmentFile(output_file, 'w', template=samfile)

iter = samfile.fetch()

for x in iter:
    if x.reference_start < 25 and x.reference_length > 80:
        outfile.write(x)

samfile.close()
outfile.close()
