import pandas as pd
from Bio import SeqIO

# =========================
# INPUT FILES
# =========================
presence_file = "Ophiostomatales_rRNA_PCGS_concat_summary_file.tsv"
fasta_file = "Ophiostomatales_rRNA_PCGs_concat_renamed.msa.fasta"

# =========================
# LOAD TABLE
# =========================
df = pd.read_csv(presence_file, sep="\t")
seq_col = df.columns[0]  # first column = sequence name

# Rename sum column for clarity
df = df.rename(columns={"sum": "gene_count"})

# =========================
# LOAD GAP COUNTS
# =========================
gap_counts = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    seq = str(record.seq)
    gap_counts[record.id] = seq.count("-") + seq.count("?")

df["gap_count"] = df[seq_col].map(gap_counts)

# =========================
# DIAGNOSTIC: sequences missing from FASTA
# =========================
missing = df[df["gap_count"].isna()]
if len(missing) > 0:
    print(f"WARNING: {len(missing)} sequences missing from FASTA. Examples:")
    print(missing[seq_col].head())

# =========================
# EXTRACT SPECIES NAME
# =========================
def get_species(name):
    parts = name.split("_")
    return "_".join(parts[:2])

df["species"] = df[seq_col].apply(get_species)

# =========================
# CREATE GENE COMBINATION STRING
# =========================
gene_cols = [
    col for col in df.columns
    if col not in [seq_col, "sum", "gene_count", "gap_count", "species"]
]

def make_gene_string(row):
    present = [g for g in gene_cols if row[g]]
    return "_".join(sorted(present)) if present else "none"

df["gene_string"] = df.apply(make_gene_string, axis=1)

# =========================
# REMOVE SUBSET COMBINATIONS
# =========================
def remove_subset_combinations(sub_df):
    combo_sets = {}
    for gs in sub_df["gene_string"].unique():
        if gs == "none":
            combo_sets[gs] = set()
        else:
            combo_sets[gs] = set(gs.split("_"))

    keep = set(combo_sets.keys())

    for gs1, set1 in combo_sets.items():
        for gs2, set2 in combo_sets.items():
            if gs1 == gs2:
                continue
            if set1 < set2:  # strict subset
                keep.discard(gs1)
                break

    return sub_df[sub_df["gene_string"].isin(keep)]

# =========================
# SELECTION LOGIC
# =========================
keep_rows = []

for species, sub in df.groupby("species"):
    # Remove subset combinations FIRST
    sub = remove_subset_combinations(sub)

    # ===== SANITY CHECK =====
    combos = sorted(sub["gene_string"].unique())
    print(f"{species}: {combos}")
    # ========================

    # Now pick best per remaining combination
    for gs, gs_group in sub.groupby("gene_string"):
        gs_group_clean = gs_group.dropna(subset=["gap_count"])

        if len(gs_group_clean) == 0:
            keep_rows.extend(gs_group.to_dict("records"))
        else:
            best = gs_group_clean.loc[gs_group_clean["gap_count"].idxmin()]
            keep_rows.append(best.to_dict())

# Convert list of dicts to DataFrame
selected_df = pd.DataFrame(keep_rows)

# =========================
# WRITE TABLE OUTPUT
# =========================
output_table = "Ophiostomatales_selected_sequences.tsv"
selected_df.to_csv(output_table, sep="\t", index=False)
print(f"Kept {len(selected_df)} sequences. Table written to {output_table}")

# =========================
# FILTER FASTA
# =========================
selected_names = set(selected_df[seq_col])
output_fasta = "Ophiostomatales_selected_sequences.fasta"

with open(output_fasta, "w") as out_fasta:
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in selected_names:
            SeqIO.write(record, out_fasta, "fasta")

print(f"Filtered FASTA written to {output_fasta}")
