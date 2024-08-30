# Identify Regulatory Variation by Long-Read Nanopore Sequencing 
# Fine Mapping Regulatory Variants by Characterizing Native CpG Methylation with Nanopore Long Read Sequencing 
Some PERL scripts I developed to parse Nanopore DNA sequenicng data.
![Picture1](https://github.com/user-attachments/assets/338a8778-32f5-424f-95b2-aa01ea7a1fe2)

Generate in-silico bisulfite conversion alignments from MM/ML BAM. The perl script is aware of indel and works for alignment report MM/ML tags on the top strand of SEQ, i.e. with C+m in MM:Z:tag. example usage:  
```
perl MMMLparse.pl example.bam output_prefix
samtools sort output_prefix.sam -o output_prefix.sorted.bam
```

