MATCH = 2
MISMATCH = -1
GAP = -2


def needleman_wunsch(seq1, seq2):

    rows = len(seq1) + 1
    cols = len(seq2) + 1

    # Create an empty scoring matrix
    matrix = [[0 for j in range(cols)] for i in range(rows)]

    # Initialize the first column
    for i in range(rows):
        matrix[i][0] = i * GAP

    # Initialize the first row
    for j in range(cols):
        matrix[0][j] = j * GAP

    # Fill the scoring matrix
    for i in range(1, rows):
        for j in range(1, cols):

            if seq1[i - 1] == seq2[j - 1]:
                diagonal_score = matrix[i - 1][j - 1] + MATCH
            else:
                diagonal_score = matrix[i - 1][j - 1] + MISMATCH

            up_score = matrix[i - 1][j] + GAP
            left_score = matrix[i][j - 1] + GAP

            matrix[i][j] = max(
                diagonal_score,
                up_score,
                left_score
            )

    # Traceback: reconstruct the alignment
    aligned_seq1 = ""
    aligned_seq2 = ""

    i = len(seq1)
    j = len(seq2)

    while i > 0 or j > 0:

        if i > 0 and j > 0:

            if seq1[i - 1] == seq2[j - 1]:
                score = MATCH
            else:
                score = MISMATCH

            if matrix[i][j] == matrix[i - 1][j - 1] + score:
                aligned_seq1 = seq1[i - 1] + aligned_seq1
                aligned_seq2 = seq2[j - 1] + aligned_seq2

                i -= 1
                j -= 1
                continue

        if i > 0 and matrix[i][j] == matrix[i - 1][j] + GAP:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2

            i -= 1

        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2

            j -= 1

    final_score = matrix[len(seq1)][len(seq2)]

    return aligned_seq1, aligned_seq2, final_score




def smith_waterman(seq1, seq2):

    rows = len(seq1) + 1
    cols = len(seq2) + 1

    matrix = [[0 for j in range(cols)] for i in range(rows)]

    max_score = 0
    max_position = (0, 0)

    # Fill the matrix
    for i in range(1, rows):
        for j in range(1, cols):

            if seq1[i - 1] == seq2[j - 1]:
                diagonal_score = matrix[i - 1][j - 1] + MATCH
            else:
                diagonal_score = matrix[i - 1][j - 1] + MISMATCH

            up_score = matrix[i - 1][j] + GAP
            left_score = matrix[i][j - 1] + GAP

            matrix[i][j] = max(
                0,
                diagonal_score,
                up_score,
                left_score
            )

            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_position = (i, j)

    # Traceback
    aligned_seq1 = ""
    aligned_seq2 = ""

    i, j = max_position

    # Save the ending position before traceback
    end_position = i

    while i > 0 and j > 0 and matrix[i][j] > 0:

        if seq1[i - 1] == seq2[j - 1]:
            score = MATCH
        else:
            score = MISMATCH

        if matrix[i][j] == matrix[i - 1][j - 1] + score:

            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2

            i -= 1
            j -= 1

        elif matrix[i][j] == matrix[i - 1][j] + GAP:

            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2

            i -= 1

        else:

            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2

            j -= 1

    # Save the starting position
    start_position = i

    return (
        aligned_seq1,
        aligned_seq2,
        max_score,
        start_position,
        end_position
    )