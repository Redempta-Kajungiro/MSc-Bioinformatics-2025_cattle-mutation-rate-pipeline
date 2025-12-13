from Bio import AlignIO

def extract_trinucleotide_context(maf_file, output_file):
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
            for column in range(1, alignment_length - 1):  # Start from 1 and end at length-1 to ensure room for context
                # Get the nucleotides for each species and the ancestral sequence, including one before and one after
                species1_trinucleotide_original = alignment[0, column-1:column+2].seq
                ancestral_trinucleotide_original = alignment[1, column-1:column+2].seq
                species2_trinucleotide_original = alignment[2, column-1:column+2].seq
                
                # Convert trinucleotides to uppercase for comparison
                species1_trinucleotide = species1_trinucleotide_original.upper()
                ancestral_trinucleotide = ancestral_trinucleotide_original.upper()
                species2_trinucleotide = species2_trinucleotide_original.upper()

                # Convert trinucleotides to strings for comparison
                species1_trinucleotide_str = str(species1_trinucleotide)
                ancestral_trinucleotide_str = str(ancestral_trinucleotide)
                species2_trinucleotide_str = str(species2_trinucleotide)

                # Extract chromosome information from sequence IDs
                species1_chromosome = alignment[0].id.split('.')[1]
                species2_chromosome = alignment[2].id.split('.')[1]

                # Check if all trinucleotides are valid (not containing gaps)
                if '-' not in species1_trinucleotide_str and '-' not in species2_trinucleotide_str and '-' not in ancestral_trinucleotide_str:
                    # Calculate the position within the alignment block including the start position
                    species1_position = species1_start + column
                    species2_position = species2_start + column

                    # Check for single nucleotide mutation in the context of trinucleotide between species 1 and the ancestral sequence
                    if species1_trinucleotide_str[1] != ancestral_trinucleotide_str[1]:
                        mutations.append((alignment[0].id, species1_chromosome, species1_position, ancestral_trinucleotide_str, species1_trinucleotide_str, str(ancestral_trinucleotide_original), str(species1_trinucleotide_original)))

                    # Check for single nucleotide mutation in the context of trinucleotide between species 2 and the ancestral sequence
                    if species2_trinucleotide_str[1] != ancestral_trinucleotide_str[1]:
                        mutations.append((alignment[2].id, species2_chromosome, species2_position, ancestral_trinucleotide_str, species2_trinucleotide_str, str(ancestral_trinucleotide_original), str(species2_trinucleotide_original)))

    # Write mutations to the output file
    with open(output_file, 'w') as out_handle:
        out_handle.write("Sequence ID\tChromosome\tPosition\tAncestral Trinucleotide\tMutation Trinucleotide\tOriginal Ancestral Trinucleotide\tOriginal Mutation Trinucleotide\n")
        for mutation in mutations:
            out_handle.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(mutation[0], mutation[1], mutation[2], mutation[3], mutation[4], mutation[5], mutation[6]))

# Example usage
maf_file = "output_CA.maf"
output_file = "trinucleotide_mutations_context.txt"
extract_trinucleotide_context(maf_file, output_file)

