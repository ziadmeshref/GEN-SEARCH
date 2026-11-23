
def find_variant(mutations, mutation):
    reference = mutation["reference"]
    sample = mutation["sample"]

    # Determine the mutation type
    if mutation["type"] == "Substitution":
        matches = mutations[
            mutations["Name"].astype(str).str.contains(
                f"{reference}>{sample}",
                regex=False,
                na=False
            )
        ]

    else:
        return None

    return matches
