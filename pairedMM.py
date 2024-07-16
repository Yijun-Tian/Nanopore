#!/app/eb/software/Python/3.8.2-GCCcore-9.3.0/bin/python
## grabnode
## module load Python/3.8.2-GCCcore-9.3.0
## source ~/venv/bin/activate
## make sure pysam version == 0.22.0

import pysam
import sys
import argparse
parser=argparse.ArgumentParser(description='''This script is used to extract paired modification tags (MM:Z:C+h;C+m ML:B:C) in a BAM file.\n
                                              Example usage: ./pairedMM.py R9P3.Haptag.mod.bam -o R9P3.pysam.out ''')

parser.add_argument('bam', help = 'A BAM file with MM and ML tags')
parser.add_argument('-o','--output', help = 'An output file with contig [reference_id], coordinate [reference_start] and modified base [modified_bases] annotation ')

args = parser.parse_args()
with open(args.output, "w") as outfile:
  samfile = pysam.AlignmentFile(args.bam, "rb")
  for read in samfile:
      outfile.write(str(read.reference_id) + "\t" + str(read.reference_start) + "\t" + str(read.modified_bases) + "\n")

  samfile.close()
outfile.close()





