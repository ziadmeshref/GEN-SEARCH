import pandas as pd
from src.loaders import load_reference
from input import get_sequence, select_gene
from src.scoring import calculate_score
from src.alignment import needleman_wunsch, smith_waterman
from src.mutation_detector import detect_mutations
from src.variant_matcher import find_variant
GENE_DATABASE = {
    "HBB": "ncbi_dataset/data/HBB/gene.fna",
    "BRCA1": "ncbi_dataset/data/BRCA1/gene.fna",
    "TP53": "ncbi_dataset/data/TP53/gene.fna"
}

sample = get_sequence()



reference_path, clinvar_path = select_gene()
reference = load_reference(reference_path)
clinvar_data = pd.read_csv(clinvar_path, sep="\t")



print("MATCH = 2, MISMATCH = -1, GAP = -2")

print("result is:",calculate_score(reference, sample),"!!!")


print("================")

if len(sample) >= len(reference) * 0.8:

    print("Using Needleman-Wunsch (Global Alignment)")

    aligned_ref, aligned_sample, score = needleman_wunsch(
        reference, sample
    )

    start = 0
    end = len(reference)

else:

    print("Using Smith-Waterman (Local Alignment)")

    aligned_ref, aligned_sample, score, start, end = smith_waterman(
        reference, sample
    )


print("\nAlignment Results:")
print("Reference:", aligned_ref)
print("Sample:   ", aligned_sample)
print("Alignment Score:", score)
print("Reference Start:", start)
print("Reference End:", end)

print("====================================")

print("\nDetecting mutations...")

mutations_found = detect_mutations(
    aligned_ref,
    aligned_sample,
    start
)

if mutations_found:

    print("\nMutations Found:")

    for mutation in mutations_found:

        print("-------------------")
        print("Position:", mutation["position"])
        print("Type:", mutation["type"])
        print("Reference:", mutation["reference"])
        print("Sample:", mutation["sample"])

        print("\nChecking ClinVar database...")

        matches = find_variant(clinvar_data, mutation)

        if matches is not None and not matches.empty:
            print(matches[[
                "Name",
                "Gene(s)",
                "Variant type",
                "Germline classification"
            ]])

        else:
            print(
                f"No matching ClinVar record found at position "
                f"{mutation['position']}"
            )

else:

    print("No mutations detected in the aligned region.")

