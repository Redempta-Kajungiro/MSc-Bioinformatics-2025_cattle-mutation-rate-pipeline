# Initialize a counter for mutations
mutation_count = 0

# Open the BED file for reading
with open("B_indi_filtered_mutations1.bed", "r") as infile:
    # Iterate over each line in the BED file
    for line in infile:
        # Increment the mutation count for each line
        mutation_count += 1

# Save the total number of mutations to a file
output_file = "B_indi_mutation_count.txt"
with open(output_file, "w") as outfile:
    outfile.write("Total number of mutations: " + str(mutation_count))

print("Mutation count saved to", output_file)

