
# Fine Mapping Regulatory Variants by Nanopore Long Read Sequencing 
Some PERL scripts to parse Nanopore DNA sequenicng data.

![Picture1](https://github.com/Yijun-Tian/Nanopore/blob/main/NanoporeMethylationUse.png)

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
TTGTACTT **__CG__** TTCAGTTA **__CG__** TATTGCTTTCTACCACACACATGCTCTTCTGTTTCCTTTTGTTCAACAGATTTCACTGGCCCATT **__CG__** TTA...AGG **__CG__** TGAGTCACCATGCCCAGCT ML:B:C,24,4,255,255
```
samtools view -F 20 output_prefix.sam | awk '{if (length($10)<=500) print $10,$NF}' | grep [CT]G
```
TTGTACTT **__TG__** TTCAGTTA **__TG__** TATTGCTTTCTACCACACACATGCTCTTCTGTTTCCTTTTGTTCAACAGATTTCACTGGCCCATT **_CG_** TAA...AGG **_CG_** TGAGTCACCATGCCCAGCT ML:B:C,24,4,255,255



**_getHaplo_SE_cgOnly.pl_**: Generate methylation pattern from in-silico converted nanopore base calling BAM files.  
This script is adapted from MONOD2 script "getHaplo_PE_cgOnly.pl"
Need to generate a CpG position file from the FASTA file, and then run sebam2cghap.sh to call the PERL script (getHaplo_SE_cgOnly.pl)
```
seqkit locate -P -p CG hg38.align.fa | awk '{if (NR>1) print $1":W\t"$5"\t"$7}' > hg38.allcpgs.txt
gzip hg38.allcpgs.txt
sebam2cghap.sh hg38.allcpgs.txt.gz output_prefix.sorted.bam > output_prefix.cgPE.hapinfo.txt
```


**_cghap2mhbs.sh_**: Define methylation haplotype block (MHB) with the methylation pattern file (output_prefix.cgPE.hapinfo.txt).
This script merges methylation pattern file into high-depth regions and call MHB with minimum LD R2 cutoff.


```
bedtools genomecov -bg -split -ibam output_prefix.sorted.bam | awk -v min=10 '$4>=min { print $1"\t"$2"\t"$3}' | bedtools merge -d 10 -i - | awk '$3-$2>80 {print $1"\t"$2"\t"$3"\t"$3-$2+1}' > output_prefix.RD10_80up.genomecov.bed
cghap2mhbs.sh output_prefix.cgPE.hapinfo.txt output_prefix.RD10_80up.genomecov.bed 0.5 dbSnp153.txt output_prefix.r2p5
```





