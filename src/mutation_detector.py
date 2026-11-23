def detect_mutations(aligned_ref, aligned_sample, start_position=0):

    mutations = []

    reference_position = start_position

    for i in range(len(aligned_ref)):

        ref_base = aligned_ref[i]
        sample_base = aligned_sample[i]

        # Insertion
        if ref_base == "-" and sample_base != "-":

            mutations.append({
                "position": reference_position,
                "type": "Insertion",
                "reference": "-",
                "sample": sample_base
            })

        # Deletion
        elif ref_base != "-" and sample_base == "-":

            mutations.append({
                "position": reference_position,
                "type": "Deletion",
                "reference": ref_base,
                "sample": "-"
            })

            reference_position += 1

        # Substitution
        elif ref_base != sample_base:

            mutations.append({
                "position": reference_position,
                "type": "Substitution",
                "reference": ref_base,
                "sample": sample_base
            })

            reference_position += 1

        else:
            reference_position += 1

    return mutations