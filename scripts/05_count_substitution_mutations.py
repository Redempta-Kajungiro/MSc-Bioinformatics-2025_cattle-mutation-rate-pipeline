# Define the order of nucleotides
nucleotides = ['A', 'C', 'G', 'T']

# Initialize a 4x4 matrix with zeros
substitution_matrix = [[0]*4 for _ in range(4)]

# Define counters for different types of substitutions
substitution_counts = {}

# Open the BED file for reading
with open("B_indi_filtered_mutations1.bed", "r") as infile:
    # Iterate over each line in the BED file
    for line in infile:
        # Split the line into columns
        columns = line.strip().split("\t")
        # Extract mutation type and ancestral nucleotide
        ancestral_nucleotide = columns[3]
        mutation_type = columns[4]

        # Update the corresponding count in the substitution_counts dictionary
        substitution = f"{ancestral_nucleotide}->{mutation_type}"
        substitution_counts[substitution] = substitution_counts.get(substitution, 0) + 1

# Populate the substitution_matrix with the counts
for i, ancestral_nucleotide in enumerate(nucleotides):
    for j, mutation_nucleotide in enumerate(nucleotides):
        substitution = f"{ancestral_nucleotide}->{mutation_nucleotide}"
        substitution_matrix[i][j] = substitution_counts.get(substitution, 0)

# Define a filename to save the results
output_filename = "B_indi_substitution_matrix.txt"

# Write the results to a file
with open(output_filename, "w") as outfile:
    # Write nucleotide labels
    outfile.write("\t" + "\t".join(nucleotides) + "\n")
    # Write the substitution matrix
    for i, row in enumerate(substitution_matrix):
        outfile.write(nucleotides[i] + "\t" + "\t".join(map(str, row)) + "\n")

# Print a message indicating successful file writing
print(f"Substitution matrix saved to {output_filename}")

