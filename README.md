
# Fine Mapping Regulatory Variants by Characterizing Native CpG Methylation with Nanopore Long Read Sequencing 
Some PERL scripts I developed to parse Nanopore DNA sequenicng data.


**_MMMLparse.pl_**: Generate in-silico bisulfite conversion alignments from MM/ML BAM. The perl script is aware of indel and works for alignment report MM/ML tags on the top strand of SEQ, i.e. with C+m in MM:Z:tag.  
  Need to set the **base modification threshold** empirically, say ML>=204, roughly 80% confident that the CpG site is modified.  
  Example usage:  
```
threshold=204
perl MMMLparse.pl example.bam output_prefix $threshold
samtools sort output_prefix.sam -o output_prefix.sorted.bam
```

Check a few short read to make sure the conversion works properly. In below alignment, first 2 CpG sites were converted to TG since the ML tag values (24, 4) are too low. 
```
samtools view -F 20 example.bam | awk '{if (length($10)<=500) print $10,$NF}' | grep CG

```
TTGTACTT **__CG__** TTCAGTTA **__CG__** TATTGCTTTCTACCACACACATGCTCTTCTGTTTCCTTTTGTTCAACAGATTTCACTGGCCCATT **__CG__** CAGAAAAATGGTAACAACCTGTTAGCTGTTTTCATCAATTTATGTGATGTATTGTGTATTAATTACTAGTATTCTCTATGTCATATTATTATTATTAGAGATGATGGAAGAGGAAGGGCATTGTATATTAATTATCAGTATATTTCATATATATATATGTATGTATGTGTATTTGTTTGAGACAAGGTCTTGCTTTGTTGCCCAACTGGAGTGCAGCAGCACTGCAGCCTTGAGCTCCTGGCCTTGAACTCCTGGCCTCAAGTGATCCTCCCACCTTGGACTCCCAAAGTGCTAGGATTACAGG **__CG__** TGAGTCACCATGCCCAGCT ML:B:C,24,4,255,255
```
samtools view -F 20 output_prefix.sam | awk '{if (length($10)<=500) print $10,$NF}' | grep [CT]G
```
TTGTACTT **__TG__** TTCAGTTA **__TG__** TATTGCTTTCTACCACACACATGCTCTTCTGTTTCCTTTTGTTCAACAGATTTCACTGGCCCATT **_CG_** CAGAAAAATGGTAACAACCTGTTAGCTGTTTTCATCAATTTATGTGATGTATTGTGTATTAATTACTAGTATTCTCTATGTCATATTATTATTATTAGAGATGATGGAAGAGGAAGGGCATTGTATATTAATTATCAGTATATTTCATATATATATATGTATGTATGTGTATTTGTTTGAGACAAGGTCTTGCTTTGTTGCCCAACTGGAGTGCAGCAGCACTGCAGCCTTGAGCTCCTGGCCTTGAACTCCTGGCCTCAAGTGATCCTCCCACCTTGGACTCCCAAAGTGCTAGGATTACAGG **_CG_** TGAGTCACCATGCCCAGCT ML:B:C,24,4,255,255







