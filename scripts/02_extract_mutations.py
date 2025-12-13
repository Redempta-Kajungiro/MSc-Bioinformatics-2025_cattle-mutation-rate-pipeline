from Bio import AlignIO

def extract_mutations(maf_file, output_file):
    mutations = []

    # Parse the MAF file
    with open(maf_file) as handle:
        alignments = AlignIO.parse(handle, "maf")
        
        # Iterate over each alignment block
        for alignment in alignments:
            # Get the length of the alignment
            alignment_length = alignment.get_alignment_length()
            
            # Get the start position of each sequence
            species1_start = alignment[0].annotations['start']
            species2_start = alignment[2].annotations['start']
                
            # Iterate over each position in the alignment block
            for column in range(alignment_length):
                # Get the nucleotides at the current position for both species and ancestral sequence
                species1_nucleotide = alignment[0, column]
                ancestral_nucleotide = alignment[1, column]
                species2_nucleotide = alignment[2, column]
                
                # Convert nucleotides to uppercase for comparison
                species1_nucleotide_upper = species1_nucleotide.upper()
                ancestral_nucleotide_upper = ancestral_nucleotide.upper()
                species2_nucleotide_upper = species2_nucleotide.upper()
                
                # Extract chromosome information from sequence IDs
                species1_chromosome = alignment[0].id.split('.')[1]
                species2_chromosome = alignment[2].id.split('.')[1]
                
                # Check if all nucleotides are valid (not gaps)
                if species1_nucleotide != '-' and species2_nucleotide != '-' and ancestral_nucleotide != '-':
                    # Calculate the position within the alignment block including the start position
                    species1_position = species1_start + column
                    species2_position = species2_start + column
                    
                    # Check for mutation between species 1 and the ancestral sequence
                    if species1_nucleotide_upper != ancestral_nucleotide_upper:
                        mutations.append((alignment[0].id, species1_chromosome, species1_position, ancestral_nucleotide, species1_nucleotide))
                    
                    # Check for mutation between species 2 and the ancestral sequence
                    if species2_nucleotide_upper != ancestral_nucleotide_upper:
                        mutations.append((alignment[2].id, species2_chromosome, species2_position, ancestral_nucleotide, species2_nucleotide))
    
    # Write mutations to the output file
    with open(output_file, 'w') as out_handle:
        out_handle.write("Sequence ID\tChromosome\tPosition\tAncestral\tMutation\tOriginal Nucleotide\n")
        for mutation in mutations:
            # Include both the original nucleotide and its uppercase version in the output
            out_handle.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(mutation[0], mutation[1], mutation[2], mutation[3], mutation[4].upper(), mutation[4]))

# Example usage
maf_file = "/path/to/input.maf"         # Input .maf or .maf.gz
output_file = "mutations_output.txt"    # Output mutation table
extract_mutations(maf_file, output_file)

