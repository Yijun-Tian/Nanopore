#!/bin/bash
# before run
# ml load BEDTools/2.30.0-GCC-11.2.0
# ml load Python/3.8.6-GCCcore-10.2.0
# ml load SAMtools/1.17-GCC-12.2.0
# Check for the correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input.vcf> <input.bam> <output.dmr.txt>"
    exit 1
fi

# Assign arguments to variables
input_vcf="$1"
input_bam="$2"
output="$3"

# Temporary file to store processed variants
temp_variants_file=$(mktemp)
touch $output
# Read the VCF line by line
while IFS= read -r line; do
    # Skip comment lines
    if [[ $line =~ ^# ]]; then
        continue
    fi

    # Extract fields from the VCF line
    IFS=$'\t' read -r chr pos id ref alt <<< "$line"

    # Adjust position by subtracting 1
    adjusted_pos=$((pos - 1))

    # Get the alternate alleles into an array
    IFS=',' read -r -a alt_array <<< "$alt"

    # Create a variant string in the required format
    variant_string="$id:$ref:${alt_array[*]}:$chr:$adjusted_pos"

    # Run the Python script for the current variant
    python label.py -i "$input_bam" -o $id.bam -v "$variant_string"
    samtools index $id.bam

    # modkit call methylation for the labeled bam ($id.bam)
    ~/dist/modkit pileup -t 2 --partition-tag VA --partition-tag RG --ref ../hg38.alignmt.fa --only-tabs --cpg $id.bam --combine-strands $id.bg/


    cd $id.bg/
    for j in $(ls *.sorted.bed | rev | cut -c 12- | rev | sort | uniq );
        do
	  awk '{FS=OFS="\t"} $4~/m/ && $5>=5 {print $1,$2,$3,$11}' $j.sorted.bed | bedtools sort -i - > $j.bg
          if [ -s $j.bg ]; then
             echo 'goodboy'
          else
             rm -f $j.bg
          fi
        done
    	grp1=$(ls ${id}@${ref}_*.bg | paste -sd,)
    	n1=$(echo $grp1 | sed 's/,/\n/g' | wc -l)
    	grp2=$(ls ${id}@${alt_array[0]}_*.bg | paste -sd,)
    	n2=$(echo $grp2 | sed 's/,/\n/g' | wc -l)
    if [ "${n1}" -ge 3 ] && [ "${n2}" -ge 3 ]; then
        perl /share/lab_wangl/Yijun_Tian/WGS_Prostate/metilene_v0.2-8/metilene_input.pl --in1 $grp1 --in2 $grp2 --out $id.ref.alt.5mCG.input --h1 ref --h2 alt
        /share/lab_wangl/Yijun_Tian/WGS_Prostate/metilene_v0.2-8/metilene_linux64 -t 4 -a ref -b alt -d 10 -M 300 -m 10 $id.ref.alt.5mCG.input > $id.ref.alt.5mCG.output
        awk -v va=$variant_string -v N1=$n1 -v N2=$n2 '{FS=OFS="\t"}{if ($7<=0.05 || $8<=0.05) print $0,va,N1,N2}' $id.ref.alt.5mCG.output >> $output
        mv $id.ref.alt.5mCG.input ../inputs/
     else
        echo "this is a mission impossible"
     fi
     cd /share/lab_wangl/Yijun_Tian/WGS_Prostate/labelAlign/
     rm -rf $id.bg/
     rm $id.bam
     rm $id.bam.bai

done < "$input_vcf"

# Clean up
rm "$temp_variants_file"

echo "Finished processing VCF"
