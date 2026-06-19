"""
PHASE 2 — TCGA BRCA TP53 LABELING (MEMORY SAFE VERSION)
=======================================================

PURPOSE
-------
1. Load MC3 PanCancer MAF safely using chunked reading
2. Extract TCGA-BRCA TP53 mutant samples
3. Generate TP53 mutant vs wildtype metadata
4. Save DESeq2-ready metadata

REQUIREMENTS
------------
Place this file here manually:

D:\TP53_project\data\expression\mc3.v0.2.8.PUBLIC.maf.gz

EXPECTED OUTPUT
---------------
~250–400 TP53 mutant samples
~700–850 TP53 wildtype samples
"""

import os
import gzip
import shutil
import pandas as pd

# ============================================================
# PATHS
# ============================================================

OUTPUT_DIR = r"D:\TP53_project\data\expression"

maf_gz_path = os.path.join(
    OUTPUT_DIR,
    "mc3.v0.2.8.PUBLIC.maf.gz"
)

maf_extract_path = os.path.join(
    OUTPUT_DIR,
    "mc3.v0.2.8.PUBLIC.maf"
)

counts_path = os.path.join(
    OUTPUT_DIR,
    "deseq2_ready_counts.csv.gz"
)

# ============================================================
# STEP 1 — CHECK FILES
# ============================================================

print("=" * 60)
print("STEP 1 — Checking MC3 MAF")
print("=" * 60)

if not os.path.exists(maf_gz_path):
    raise FileNotFoundError(
        f"\nMC3 file not found:\n{maf_gz_path}"
    )

print("\nMC3 gzip file found.")

# ============================================================
# STEP 2 — EXTRACT MAF
# ============================================================

print("\n" + "=" * 60)
print("STEP 2 — Extracting MAF")
print("=" * 60)

if not os.path.exists(maf_extract_path):

    print("\nExtracting...")

    with gzip.open(maf_gz_path, "rb") as f_in:
        with open(maf_extract_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    print("Extraction complete!")

else:
    print("\nExtracted MAF already exists.")

# ============================================================
# STEP 3 — LOAD BRCA MAF (MEMORY SAFE)
# ============================================================

print("\n" + "=" * 60)
print("STEP 3 — Loading BRCA rows from MAF")
print("=" * 60)

chunksize = 50000

usecols = [
    "Hugo_Symbol",
    "Tumor_Sample_Barcode",
    "Variant_Classification",
    "Variant_Type",
    "HGVSp_Short",
    "Chromosome",
    "Start_Position"
]

brca_chunks = []

print("\nReading MAF in chunks...\n")

for i, chunk in enumerate(pd.read_csv(
    maf_extract_path,
    sep="\t",
    comment="#",
    usecols=usecols,
    chunksize=chunksize,
    low_memory=True
)):

    # Keep TCGA samples only
    tcga_chunk = chunk[
        chunk["Tumor_Sample_Barcode"].str.startswith("TCGA", na=False)
    ].copy()

    brca_chunks.append(tcga_chunk)

    if i % 10 == 0:
        print(f"Processed chunk {i}")

maf_df = pd.concat(brca_chunks, ignore_index=True)

print("\nMAF loaded successfully!")
print("Shape:", maf_df.shape)

print("\nExample barcodes:")
print(maf_df["Tumor_Sample_Barcode"].head())

# ============================================================
# STEP 4 — FILTER TP53 NON-SILENT MUTATIONS
# ============================================================

print("\n" + "=" * 60)
print("STEP 4 — Extracting TP53 mutations")
print("=" * 60)

silent_classes = [
    "Silent",
    "3'UTR",
    "5'UTR",
    "3'Flank",
    "5'Flank",
    "Intron",
    "RNA",
    "IGR"
]

tp53_maf = maf_df[
    (maf_df["Hugo_Symbol"] == "TP53") &
    (~maf_df["Variant_Classification"].isin(silent_classes))
].copy()

print("\nTP53 variants found:", len(tp53_maf))

# ============================================================
# STEP 5 — EXTRACT CASE IDS
# ============================================================

print("\n" + "=" * 60)
print("STEP 5 — Building mutant sample list")
print("=" * 60)

tp53_maf["case_barcode"] = (
    tp53_maf["Tumor_Sample_Barcode"]
    .astype(str)
    .str[:12]
)

tp53_mutant_cases = set(
    tp53_maf["case_barcode"].unique()
)

print("\nUnique TP53 mutant cases:")
print(len(tp53_mutant_cases))

print("\nExample mutant cases:")
print(list(tp53_mutant_cases)[:10])

# ============================================================
# STEP 6 — SAVE TP53 VARIANTS
# ============================================================

tp53_output = os.path.join(
    OUTPUT_DIR,
    "tp53_variants_from_maf.csv"
)

tp53_maf.to_csv(tp53_output, index=False)

print("\nTP53 variant table saved:")
print(tp53_output)

# ============================================================
# STEP 7 — LOAD COUNT MATRIX SAMPLE IDS
# ============================================================

print("\n" + "=" * 60)
print("STEP 7 — Loading count matrix sample IDs")
print("=" * 60)

counts_header = pd.read_csv(
    counts_path,
    compression="gzip",
    nrows=0
)

sample_barcodes = list(counts_header.columns)

print("\nSamples in count matrix:")
print(len(sample_barcodes))

# ============================================================
# STEP 8 — BUILD METADATA
# ============================================================

print("\n" + "=" * 60)
print("STEP 8 — Building metadata")
print("=" * 60)

metadata_rows = []

for barcode in sample_barcodes:

    case_id = str(barcode)[:12]

    tp53_status = (
        "TP53_mutant"
        if case_id in tp53_mutant_cases
        else "TP53_wildtype"
    )

    metadata_rows.append({
        "sample_barcode": barcode,
        "case_id": case_id,
        "TP53_status": tp53_status
    })

meta_df = pd.DataFrame(metadata_rows)

# ============================================================
# STEP 9 — SUMMARY
# ============================================================

mutant_n = (
    meta_df["TP53_status"] == "TP53_mutant"
).sum()

wildtype_n = (
    meta_df["TP53_status"] == "TP53_wildtype"
).sum()

print("\nTP53-mutant samples :", mutant_n)
print("TP53-wildtype samples :", wildtype_n)

# ============================================================
# STEP 10 — SAVE METADATA
# ============================================================

print("\n" + "=" * 60)
print("STEP 10 — Saving outputs")
print("=" * 60)

metadata_path = os.path.join(
    OUTPUT_DIR,
    "sample_metadata.csv"
)

deseq_metadata_path = os.path.join(
    OUTPUT_DIR,
    "deseq2_ready_metadata.csv"
)

meta_df.to_csv(metadata_path, index=False)
meta_df.to_csv(deseq_metadata_path, index=False)

print("\nMetadata saved:")
print(metadata_path)

print("\nDESeq2 metadata saved:")
print(deseq_metadata_path)

# ============================================================
# STEP 11 — SAVE BARCODE LISTS
# ============================================================

mutant_txt = os.path.join(
    OUTPUT_DIR,
    "tp53_mutant_barcodes.txt"
)

wildtype_txt = os.path.join(
    OUTPUT_DIR,
    "tp53_wildtype_barcodes.txt"
)

with open(mutant_txt, "w") as f:
    for b in meta_df.loc[
        meta_df["TP53_status"] == "TP53_mutant",
        "sample_barcode"
    ]:
        f.write(str(b) + "\n")

with open(wildtype_txt, "w") as f:
    for b in meta_df.loc[
        meta_df["TP53_status"] == "TP53_wildtype",
        "sample_barcode"
    ]:
        f.write(str(b) + "\n")

print("\nBarcode lists saved.")

# ============================================================
# FINAL SUMMARY
# ============================================================

print("\n" + "=" * 60)
print("COMPLETE")
print("=" * 60)

print(f"\nTP53 mutant samples  : {mutant_n}")
print(f"TP53 wildtype samples: {wildtype_n}")

print("\nReady for DESeq2:")
print("  deseq2_ready_counts.csv.gz")
print("  deseq2_ready_metadata.csv")
print("=" * 60)