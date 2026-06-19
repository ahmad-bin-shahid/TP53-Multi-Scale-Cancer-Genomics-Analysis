import pandas as pd

# File path
file_path = r"D:\TP53_project\data\mutations\variant_summary.txt.gz"

# Empty list to store TP53 rows
tp53_chunks = []

# Read file in chunks
chunk_size = 50000

for chunk in pd.read_csv(
    file_path,
    sep='\t',
    compression='gzip',
    chunksize=chunk_size,
    low_memory=False
):
    
    # Filter TP53 rows only
    tp53_chunk = chunk[chunk['GeneSymbol'] == 'TP53']
    
    # Store filtered rows
    tp53_chunks.append(tp53_chunk)

# Combine all TP53 rows
tp53_df = pd.concat(tp53_chunks, ignore_index=True)

# Display results
print("TP53 rows extracted successfully!\n")

print("Shape:")
print(tp53_df.shape)

print("\nColumns:")
print(tp53_df.columns.tolist())

# Inspect unique clinical significance values
print("\nUnique Clinical Significance Values:\n")
print(tp53_df['ClinicalSignificance'].unique())

# Inspect unique variant types
print("\nUnique Variant Types:\n")
print(tp53_df['Type'].unique())

# Define accepted clinical significance values
accepted_significance = [
    'Pathogenic',
    'Likely pathogenic',
    'Pathogenic/Likely pathogenic'
]

# Filter clinically relevant variants
filtered_df = tp53_df[
    tp53_df['ClinicalSignificance'].isin(accepted_significance)
]

# Keep only single nucleotide variants
filtered_df = filtered_df[
    filtered_df['Type'] == 'single nucleotide variant'
]

# Reset index
filtered_df = filtered_df.reset_index(drop=True)

# Display filtering results
print("\nFiltered TP53 dataset shape:")
print(filtered_df.shape)

print("\nRemaining Clinical Significance values:")
print(filtered_df['ClinicalSignificance'].unique())

print("\nRemaining Variant Types:")
print(filtered_df['Type'].unique())

# Save filtered dataset
output_path = r"D:\TP53_project\data\mutations\filtered_tp53_variants.csv"

filtered_df.to_csv(output_path, index=False)

print("\nFiltered dataset saved successfully!")
print(f"Saved to: {output_path}")