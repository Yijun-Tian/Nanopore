import pysam
import argparse
import json

def parse_variants(variant_str):
    """Parse a string of variants into a list of tuples."""
    variants = []
    for variant in variant_str.split(","):
        variant_data = variant.split(":")
        if len(variant_data) < 5:
            raise ValueError("Each variant must include id, ref, alt1, alt2, chromosome, and start coordinate.")
        variant_id = variant_data[0]
        ref_allele = variant_data[1]
        alt_alleles = variant_data[2:-2]  # Capture all alt alleles until chromosome
        chromosome = variant_data[-2]
        start_coordinate = int(variant_data[-1])
        variants.append((variant_id, ref_allele, alt_alleles, chromosome, start_coordinate))
    return variants

def main(input_bam, output_bam, variants):
    # Open the BAM file
    bam_in = pysam.AlignmentFile(input_bam, "rb")
    bam_out = pysam.AlignmentFile(output_bam, "wb", header=bam_in.header)

    # Iterate through each variant
    for variant_id, ref_allele, alt_alleles, chromosome, start in variants:
        stop = start + 1  # Assuming variants are point mutations

        # Iterate through each read in the BAM file that overlaps with the variant position
        for read in bam_in.fetch(chromosome, start, stop):
            if read.query_sequence is None:  # Skip if read has no sequence
                continue
            
            # Get the positions the read covers
            read_positions = read.get_reference_positions(full_length=True)

            # Check if the variant position is within the read's positions
            if start in read_positions:
                read_index = read_positions.index(start)
                read_base = read.query_sequence[read_index]

                # Add tag information
                tag_info = f"{variant_id}@{read_base}"

                # Check if the read base matches the reference or any alternate alleles
                matched = False
                
                if read_base == ref_allele:
                    matched = True
                else:
                    # Check against all alternate alleles
                    for alt in alt_alleles:
                        if read_base == alt:
                            matched = True
                            break
                
                # Record as is for unexpected alleles if no match
                if not matched:
                    tag_info = f"{variant_id}@{read_base}"  # Still capture the unexpected allele

                # If the read already has a VA tag, append the new tag information
                if read.has_tag("VA"):
                    existing_tag = read.get_tag("VA")
                    new_tag = f"{existing_tag};{tag_info}"  # Use semicolon to separate entries
                    read.set_tag("VA", new_tag)
                else:
                    read.set_tag("VA", tag_info)  # Set the VA tag for the first time

             # Only write the read if it has a VA tag
            if read.has_tag("VA"):
                bam_out.write(read)

    # Close the BAM and VCF files
    bam_in.close()
    bam_out.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate BAM reads with variant allele tags.")
    parser.add_argument("-i", "--input_bam", required=True, help="Input BAM file")
    parser.add_argument("-o", "--output_bam", required=True, help="Output BAM file")
    parser.add_argument("-v", "--variants", required=True,
                        help="Variants in the format: id:ref:alt1:alt2:chromosome:start, ...")
    
    args = parser.parse_args()
    
    # Parse the variants string into a list of tuples
    variant_list = parse_variants(args.variants)

    # Call the main function with the arguments
    main(args.input_bam, args.output_bam, variant_list)
