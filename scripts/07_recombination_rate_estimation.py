import pandas as pd


# Input file 
input_file = 'input_map_file.txt'

# Column names in your file 
marker_col = 'marker'
chr_col = 'chr'
bp_col = 'position_bp'
cm_col = 'position_cM'

# Output file name
output_file = 'categorized_recombination_rates.txt'

# Load data

df = pd.read_csv(input_file, sep='\t')

# Rename columns internally so the rest of the script stays identical
df = df.rename(columns={
    marker_col: 'marker',
    chr_col: 'chr',
    bp_col: 'position_bp',
    cm_col: 'position_cM'
})

# Calculate recombination rate
# cM per bp
df['recombination_rate_cM_per_bp'] = df['position_cM'].diff() / df['position_bp'].diff()

# Remove the first NA row
df = df.dropna().reset_index(drop=True)

# Convert to cM per Mb
df['recombination_rate_cM_per_Mb'] = df['recombination_rate_cM_per_bp'] * 1e6

# Drop any remaining NA
df = df.dropna().reset_index(drop=True)

# Sort by recombination rate (not required but kept as original)
df_sorted = df.sort_values(by='recombination_rate_cM_per_Mb', ascending=False).reset_index(drop=True)

# Categorize (High / Medium / Low)
top_20_percentile = df['recombination_rate_cM_per_Mb'].quantile(0.8)
bottom_20_percentile = df['recombination_rate_cM_per_Mb'].quantile(0.2)

def categorize_recombination(rate):
    if rate >= top_20_percentile:
        return 'High'
    elif rate <= bottom_20_percentile:
        return 'Low'
    else:
        return 'Medium'

df['recombination_category'] = df['recombination_rate_cM_per_Mb'].apply(categorize_recombination)

# Save output
output_df = df[['marker', 'chr', 'position_bp', 'position_cM',
                'recombination_rate_cM_per_bp', 'recombination_rate_cM_per_Mb',
                'recombination_category']]

output_df.to_csv(output_file, sep='\t', index=False)

print(f"Output saved to {output_file}")

