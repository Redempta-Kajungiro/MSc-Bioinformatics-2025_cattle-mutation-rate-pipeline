from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import os
import gzip

# Folder containing MAF or MAF.GZ files
input_folder = "/path/to/maf_files"

# Output filtered file
output_file = "/path/to/output/filtered_output.maf"

# List of required species identifiers (partial name allowed)
required_species = [
    "species1_keyword",
    "species2_keyword",
    "species3_keyword"
]

filtered_alignments = []
files_checked = []

# Loop through all files in the folder
for filename in os.listdir(input_folder):
    full_path = os.path.join(input_folder, filename)

    # Skip non-MAF files
    if not (filename.endswith(".maf") or filename.endswith(".maf.gz")):
        continue

    # Open .maf or .maf.gz
    if filename.endswith(".gz"):
        handle = gzip.open(full_path, "rt")
    else:
        handle = open(full_path, "r")

    # Read all alignment blocks in the file
    for alignment in AlignIO.parse(handle, "maf"):

        # Temporary storage for required species sequences
        found_sequences = {}

        # Check each sequence in the alignment block
        for seq_rec in alignment:
            for sp in required_species:
                if sp.lower() in seq_rec.id.lower():  # case-insensitive matching
                    found_sequences[sp] = seq_rec

        # Keep block only if all required species are present
        if len(found_sequences) == len(required_species):
            # Preserve order of required_species list
            ordered_seqs = [found_sequences[sp] for sp in required_species]
            filtered_alignments.append(MultipleSeqAlignment(ordered_seqs))

    handle.close()
    files_checked.append(filename)

# Write all filtered alignments to output file
with open(output_file, "w") as out:
    AlignIO.write(filtered_alignments, out, "maf")

# Optional: print files processed
print("Processed files:")
for f in files_checked:
    print(" -", f)

print(f"\nTotal alignments kept: {len(filtered_alignments)}")
print(f"Filtered MAF saved to: {output_file}")
