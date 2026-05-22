#!/usr/bin/env python3
"""
Find matching sequences in a reference FASTA based on query FASTA headers.

- First two underscore-separated tokens = organism (fuzzy match), unless the
  second token looks like a strain/accession code
- Tokens like sp1293/spR071 are treated as "sp." plus the code 1293/R071
- If an accession is included in the query name, accession matching is tried first
- Genus names can differ by one character if the species name still matches exactly
- Taxon spelling spr is treated as spor
- Other tokens = strain/isolate (exact preferred, punctuation-insensitive fallback)
- Outputs CSV: query_name, accession, header

Usage:
    python match_fasta_ref.py query.fasta reference.fasta output.csv
"""

import sys, csv, re

# ---------------- FUNCTIONS ----------------

def split_strain_token(token):
    """Split alphanumeric strain codes like CXY4044 into ['CXY', '4044']"""
    m = re.match(r"^([A-Za-z]+)(\d+)$", token)
    if m:
        return list(m.groups())
    else:
        return [token]

def normalize_for_match(text):
    """
    Normalize headers/strain names so formatting differences do not block matches.

    Example:
        DTO401F9, DTO 401-F9, and DTO_401_F9 all become dto401f9
    """
    return re.sub(r"[^a-z0-9]+", "", text.lower())

def edit_distance_at_most_one(a, b):
    """Return True if two strings are identical or one edit apart."""
    if a == b:
        return True
    if abs(len(a) - len(b)) > 1:
        return False

    if len(a) > len(b):
        a, b = b, a

    i = j = edits = 0
    while i < len(a) and j < len(b):
        if a[i] == b[j]:
            i += 1
            j += 1
            continue

        edits += 1
        if edits > 1:
            return False

        if len(a) == len(b):
            i += 1
        j += 1

    return True

def canonical_taxon_word(word):
    """Normalize known taxon spelling variants."""
    word = normalize_for_match(word)
    return word.replace("spr", "spor")

def organism_matches(org_pattern, organism_tokens, header_lower):
    """Match organism exactly first, then allow one genus typo with exact species."""
    if org_pattern.search(header_lower):
        return True

    if len(organism_tokens) == 1:
        genus = canonical_taxon_word(organism_tokens[0])
        header_words = [
            canonical_taxon_word(word)
            for word in re.findall(r"[A-Za-z0-9]+", header_lower)
        ]
        return genus in header_words

    if len(organism_tokens) < 2:
        return False

    genus = canonical_taxon_word(organism_tokens[0])
    species = canonical_taxon_word(organism_tokens[1])
    header_words = [
        canonical_taxon_word(word)
        for word in re.findall(r"[A-Za-z0-9]+", header_lower)
    ]

    for i in range(len(header_words) - 1):
        if header_words[i + 1] == species and edit_distance_at_most_one(genus, header_words[i]):
            return True

    return False

def token_looks_like_code(token):
    """Return True for strain/accession-like tokens such as CBS10319."""
    return any(c.isdigit() for c in token)

def split_sp_token(token):
    """Split tokens like sp1293 or spR071 into ('sp', '1293'/'R071')."""
    m = re.match(r"^(sp)([A-Za-z0-9]+)$", token, re.IGNORECASE)
    if m:
        return m.group(1), m.group(2)
    return None

def accession_from_query_name(query_name):
    """
    Return an accession-like value from a query name, if present.

    Example:
        HQ608102_1 -> hq6081021, matching HQ608102.1 after normalization
    """
    tokens = [t for t in re.split(r"[^A-Za-z0-9]+", query_name) if t]
    for i, token in enumerate(tokens):
        # Common GenBank-style nucleotide accessions: AY649780, JX244068, HQ608102.
        if re.match(r"^[A-Za-z]{2}\d{6}$", token):
            if i + 1 < len(tokens) and tokens[i + 1].isdigit():
                return normalize_for_match(token + tokens[i + 1])
            return normalize_for_match(token)

    return None

def accession_values(text):
    """Return normalized full/base accession forms for a reference ID/header token."""
    values = {normalize_for_match(text)}
    m = re.match(r"^([A-Za-z]{1,4}_?\d+)[._](\d+)", text)
    if m:
        values.add(normalize_for_match(m.group(1)))
        values.add(normalize_for_match(m.group(1) + m.group(2)))
    return values

def find_accession_match(accession_normalized, ref_records):
    """Find a reference record by normalized accession."""
    if not accession_normalized:
        return None

    for r in ref_records:
        ref_accession = r.description.split()[0] if r.description.split() else ""
        candidates = accession_values(r.id) | accession_values(ref_accession)
        if accession_normalized in candidates:
            return r

    return None

def is_unverified_record(record):
    """Return True for reference records marked as unverified."""
    return "unverified" in record.description.lower()

def query_parse_candidates(parts):
    """
    Build possible organism/strain parses from a query FASTA header.

    Most queries are Genus_species_strain, but some are Genus_strain_accession.
    Try the likely parse first, then fall back to the older two-token organism
    parse when useful.
    """
    if len(parts) < 2:
        return [([parts[0]], parts[1:])]

    sp_split = split_sp_token(parts[1])
    if sp_split:
        sp_token, code_token = sp_split
        return [
            ([parts[0], sp_token], [code_token] + parts[2:]),
            ([parts[0]], parts[1:]),
            (parts[:2], parts[2:]),
        ]

    if token_looks_like_code(parts[1]):
        return [
            ([parts[0]], parts[1:]),
            (parts[:2], parts[2:]),
        ]

    return [(parts[:2], parts[2:])]

def find_match_for_tokens(organism_tokens, strain_tokens, ref_records):
    """
    Find best matching sequence from reference FASTA for one query parse.
    """
    # Regex for organism: allow any non-word chars between the two words
    if len(organism_tokens) == 1:
        org_pattern = re.compile(re.escape(organism_tokens[0]), re.IGNORECASE)
    else:
        org_pattern = re.compile(
            rf"{re.escape(organism_tokens[0])}\W+{re.escape(organism_tokens[1])}",
            re.IGNORECASE
        )

    strain_exact = " ".join(strain_tokens).lower()
    strain_normalized = normalize_for_match(" ".join(strain_tokens))
    strain_tokens_split = []
    for t in strain_tokens:
        strain_tokens_split.extend(split_strain_token(t.lower()))
    strain_tokens_split_normalized = [
        normalize_for_match(t) for t in strain_tokens_split
        if normalize_for_match(t)
    ]

    best_partial = None

    for r in ref_records:
        header_lower = r.description.lower().replace("_", " ")
        header_normalized = normalize_for_match(r.description)
        exact_strain_matches = strain_exact and strain_exact in header_lower
        normalized_strain_matches = (
            strain_normalized and strain_normalized in header_normalized
        )
        partial_strain_matches = (
            strain_tokens_split
            and all(t in header_lower for t in strain_tokens_split)
        )
        normalized_partial_strain_matches = (
            strain_tokens_split_normalized
            and all(t in header_normalized for t in strain_tokens_split_normalized)
        )

        # Organism must match, except for UNVERIFIED records where the strain/code
        # is sometimes the only useful identifier in the header.
        if not organism_matches(org_pattern, organism_tokens, header_lower):
            if is_unverified_record(r) and (
                normalized_strain_matches or normalized_partial_strain_matches
            ):
                best_partial = r
            continue

        # Exact strain match preferred
        if exact_strain_matches:
            return r

        # Formatting-insensitive strain match, e.g. DTO401F9 vs DTO 401-F9
        if normalized_strain_matches:
            return r

        # Partial strain token fallback
        if partial_strain_matches:
            best_partial = r

        # Punctuation-insensitive partial fallback, e.g. CBS10319 vs CBS 103.19
        if normalized_partial_strain_matches:
            best_partial = r

    return best_partial

def match_sequence(query_name, ref_records):
    """
    Find best matching sequence from reference FASTA.
    """
    parts = query_name.split("_")
    accession_match = find_accession_match(
        accession_from_query_name(query_name),
        ref_records
    )
    if accession_match:
        return accession_match

    for organism_tokens, strain_tokens in query_parse_candidates(parts):
        best = find_match_for_tokens(organism_tokens, strain_tokens, ref_records)
        if best:
            return best

    return None

# ---------------- MAIN ----------------

def main():
    try:
        from Bio import SeqIO
    except ImportError:
        print("Error: Biopython is required. Install it with: pip install biopython")
        sys.exit(1)

    if len(sys.argv) != 4:
        print("Usage: python match_fasta_ref.py query.fasta reference.fasta output.csv")
        sys.exit(1)

    query_fasta = sys.argv[1]
    ref_fasta = sys.argv[2]
    output_csv = sys.argv[3]

    print(f"Loading reference FASTA: {ref_fasta} ...")
    ref_records = list(SeqIO.parse(ref_fasta, "fasta"))
    print(f"Loaded {len(ref_records)} reference sequences.\n")

    print(f"Processing query FASTA: {query_fasta}\n")
    total = sum(1 for _ in SeqIO.parse(query_fasta, "fasta"))

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["query_name", "accession", "sequence_title"])

        for i, record in enumerate(SeqIO.parse(query_fasta, "fasta"), start=1):
            query_name = record.id
            print(f"[{i}/{total}] {query_name} ...", end=" ")

            best = match_sequence(query_name, ref_records)
            if best:
                accession = best.id.split()[0]  # first word = accession
                writer.writerow([query_name, accession, best.description])
                print(f"Matched: {accession}")
            else:
                writer.writerow([query_name, "NO_MATCH", ""])
                print("NO MATCH")

    print(f"\nDone! Results saved to: {output_csv}")

if __name__ == "__main__":
    main()
