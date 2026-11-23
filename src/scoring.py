MATCH = 2
MISMATCH = -1
GAP = -2


def calculate_score(seq1, seq2):

    score = 0

    for a, b in zip(seq1, seq2):

        if a == "-" or b == "-":
            score += GAP

        elif a == b:
            score += MATCH

        else:
            score += MISMATCH

    return score


