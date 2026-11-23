def get_sequence():

    print("DNA Mutation Detection System")
    print("Enter your DNA sequence.")
    print("Press ENTER on an empty line when finished:\n")

    lines = []

    while True:
        line = input().strip()

        if not line:
            break

        if line.startswith(">"):
            continue

        lines.append(line)

    sequence = "".join(lines).upper()

    allowed = set("ATCGN")

    if not sequence:
        print("Error: Empty sequence")
        return None

    if not set(sequence).issubset(allowed):
        print("Error: Invalid DNA sequence")
        return None

    print("\nSequence accepted!")
    print("Sequence length:", len(sequence))

    return sequence


GENE_DATABASE = {
    "HBB": "ncbi_dataset/data/HBB/gene.fna",
    "BRCA1": "ncbi_dataset/BRCA1/BRCA1.fna",
    "TP53": "ncbi_dataset/TP53/TP53.fna"
}

CLINVAR_DATABASE = {
    "HBB": "ncbi_dataset/data/HBB/clinvar_download_20261114_160109.tsv",
    "BRCA1": "ncbi_dataset/BRCA1/clinvar_download_20261116_183351.tsv",
    "TP53": "ncbi_dataset/TP53/clinvar_download_20261116_183510.tsv"
}
def select_gene():

    print("\nAvailable Genes:")

    for gene in GENE_DATABASE:
        print("-", gene)

    gene = input("\nEnter gene name: ").strip().upper()

    if gene not in GENE_DATABASE:
        print("Gene not found!")
        return None
    


    return GENE_DATABASE[gene], CLINVAR_DATABASE[gene]
