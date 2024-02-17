from Bio import pairwise2
import numpy as np
from matplotlib import pyplot as plt

def score(str1, str2):
    score=0
    for i in range(len(str1)):
        if str1[i] == str2[i]:
            score += match
        elif str1[i] == "-" or str2[i] == "-":
            score += gap_t
        elif str1[i] != str2[i]:
            score += mismatch
    print("str1: ", str1)
    print("str2: ", str2)
    print("score = ", score)


def pairwise_aligment(str_1, str_2):
    for a in pairwise2.align.globalms(str_1, str_2, match, mismatch, gap_t, extend):
        print(a)

def threeAligment():
    def printMatrix(m, s1, s2, s3):
        rows = len(m)
        columns = len(m[0])
        depth = len(m[0][0])
        x = []
        y = []
        z = []
        c = []
        for i in range(1, rows):
            for j in range(1, columns):
                for k in range(1, depth):
                    x.append(i)
                    y.append(j)
                    z.append(k)
                    c.append(m[i][j][k])

        plt.rcParams["figure.figsize"] = [4, 4]
        plt.rcParams["figure.autolayout"] = True
        fig = plt.figure(figsize=(8, 8), dpi=80)
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xticks(np.arange(len(s1) + 1))
        ax.set_yticks(np.arange(len(s2) + 1))
        ax.set_zticks(np.arange(len(s3) + 1))
        tx = list(s1)
        ty = list(s2)
        tz = list(s3)
        tx.insert(0, "*")
        ty.insert(0, "*")
        tz.insert(0, "*")
        ax.set_xticklabels(tx)
        ax.set_yticklabels(ty)
        ax.set_zticklabels(tz)
        # Use this to ass scores lables on the cube
        # for i in range(len(x)):
        #  ax.text(x[i],y[i],z[i],  '%s' % (str(c[i])), size=10, zorder=1, color='k')
        sc = ax.scatter(x, y, z, color='b', s=3)  # [i])
        return plt, ax, fig, sc

    alphabet = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u",
                "v", "w", "x", "y", "z", "è", "é", "ò", "à", "ù", "ì", " "]
    # Edit the substitution matrix as you like
    sub_matrix = [[0 for col in alphabet] for row in alphabet]
    for i in range(len(alphabet)):
        for j in range(len(alphabet)):
            if alphabet[i] == alphabet[j]:
                sub_matrix[i][j] = match
            else:
                sub_matrix[i][j] = mismatch

    # 3D Dynamic Programming

    matrix = [[[0 for col in range(len(s3) + 1)] for row in range(len(s2) + 1)] for depth in range(len(s1) + 1)]
    print()
    # Score calculation
    max_position = (0, 0, 0)
    max_value = -1
    gap = -gap_t
    rows = len(matrix)  # s2
    columns = len(matrix[0])  # s3
    depth = len(matrix[0][0])  # s1

    # EDIT Here to make it Smith and Waterman
    i = 0
    for j in range(0, columns):  # s3
        for k in range(0, depth):  # s1
            matrix[i][j][k] = -1 - j - k
    j = 0
    for i in range(0, rows):  # s3
        for k in range(0, depth):  # s1
            matrix[i][j][k] = -1 - i - k
    k = 0
    for j in range(0, columns):  # s3
        for i in range(0, rows):  # s1
            matrix[i][j][k] = -1 - j - i

    matrix[0][0][0] = 0

    for i in range(1, rows):  # s2
        for j in range(1, columns):  # s3
            for k in range(1, depth):  # s1
                # On the same layer
                up = matrix[i - 1][j][k] - gap - gap
                left = matrix[i][j - 1][k] - gap - gap
                diag = matrix[i - 1][j - 1][k] + sub_matrix[alphabet.index(s1[i - 1])][alphabet.index(s2[j - 1])] - gap
                # Previous layer
                back_up = matrix[i - 1][j][k - 1] + sub_matrix[alphabet.index(s1[i - 1])][
                    alphabet.index(s3[k - 1])] - gap
                back_left = matrix[i][j - 1][k - 1] + sub_matrix[alphabet.index(s2[j - 1])][
                    alphabet.index(s3[k - 1])] - gap
                back_diag = matrix[i - 1][j - 1][k - 1] + (
                        sub_matrix[alphabet.index(s1[i - 1])][alphabet.index(s2[j - 1])] +
                        sub_matrix[alphabet.index(s1[i - 1])][alphabet.index(s3[k - 1])] +
                        sub_matrix[alphabet.index(s2[j - 1])][alphabet.index(s3[k - 1])]
                )
                back = matrix[i][j][k - 1] - gap - gap
                matrix[i][j][k] = max(up, left, diag, back_diag, back, back_up, back_left)
                if matrix[i][j][k] >= max_value:
                    max_position = (i, j, k)
                    max_value = matrix[i][j][k]

    final_plot, final_ax, final_fig, final_sc = printMatrix(matrix, s1, s2, s3)

    max_position = (0, 0, 0)
    max_value = -1
    for r in range(rows):
        for c in range(columns):
            if matrix[r][c][depth - 1] >= max_value:
                max_position = (r, c, depth - 1)
                max_value = matrix[r][c][depth - 1]

    for d in range(depth):
        for c in range(columns):
            if matrix[rows - 1][c][d] >= max_value:
                max_position = (rows - 1, c, d)
                max_value = matrix[rows - 1][c][d]

    for d in range(depth):
        for r in range(rows):
            if matrix[r][columns - 1][d] >= max_value:
                max_position = (r, columns - 1, d)
                max_value = matrix[r][columns - 1][d]

    r, c, d = (rows - 1, columns - 1, depth - 1)  # max_position#
    l1 = ""
    l2 = ""
    l3 = ""
    while True:
        up = matrix[r - 1][c][d]
        left = matrix[r][c - 1][d]
        diag = matrix[r - 1][c - 1][d]

        back = matrix[r][c][d - 1]
        back_up = matrix[r - 1][c][d - 1]
        back_left = matrix[r][c - 1][d - 1]
        back_diag = matrix[r - 1][c - 1][d - 1]
        if (back_diag >= up and
                back_diag >= diag and
                back_diag >= back and
                back_diag >= back_up and
                back_diag >= back_left and
                back_diag >= left):
            final_plot.quiver(
                r, c, d,  # <-- starting point of vector
                -1, -1, -1,  # <-- directions of vector
                color='red', alpha=.8, lw=2,
            )
            l1 = s1[r - 1] + l1
            l2 = s2[c - 1] + l2
            l3 = s3[d - 1] + l3
            r = r - 1
            c = c - 1
            d = d - 1

        elif (back_left >= up and
              back_left >= diag and
              back_left >= back and
              back_left >= back_diag and
              back_left >= back_up and
              back_left >= left):
            final_plot.quiver(
                r, c, d,  # <-- starting point of vector
                0, -1, -1,  # <-- directions of vector
                color='red', alpha=.8, lw=2,
            )
            l1 = "-" + l1
            l2 = s2[c - 1] + l2
            l3 = s3[d - 1] + l3
            c = c - 1
            d = d - 1
        elif (back_up >= up and
              back_up >= diag and
              back_up >= back and
              back_up >= back_diag and
              back_up >= back_left and
              back_up >= left):
            final_plot.quiver(
                r, c, d,  # <-- starting point of vector
                -1, 0, -1,  # <-- directions of vector
                color='red', alpha=.8, lw=2,
            )
            l1 = s1[r - 1] + l1
            l2 = "-" + l2
            l3 = s3[d - 1] + l3
            r = r - 1
            d = d - 1
        elif (diag >= left and
              diag >= up and
              diag >= back and
              diag >= back_up and
              diag >= back_left and
              diag > back_diag):
            final_plot.quiver(
                r, c, d,  # <-- starting point of vector
                -1, -1, 0,  # <-- directions of vector
                color='red', alpha=.8, lw=2,
            )
            l1 = s1[r - 1] + l1
            l2 = s2[c - 1] + l2
            l3 = "-" + l3
            r = r - 1
            c = c - 1

        elif (back >= up and
              back >= diag and
              back >= back_up and
              back >= back_diag and
              back >= back_left and
              back >= left):
            final_plot.quiver(
                r, c, d,  # <-- starting point of vector
                0, 0, -1,  # <-- directions of vector
                color='red', alpha=.8, lw=2,
            )
            l1 = "-" + l1
            l2 = "-" + l2
            l3 = s3[d - 1] + l3
            d = d - 1
        elif (up >= left and
              up >= diag and
              up >= back and
              up >= back_up and
              up >= back_left and
              up > back_diag):
            final_plot.quiver(
                r, c, d,  # <-- starting point of vector
                -1, 0, 0,  # <-- directions of vector
                color='red', alpha=.8, lw=2,
            )
            l1 = s1[r - 1] + l1
            l2 = "-" + l2
            l3 = "-" + l3
            r = r - 1
        elif (left >= up and
              left >= diag and
              left >= back and
              left >= back_up and
              left >= back_left and
              left > back_diag):
            final_plot.quiver(
                r, c, d,  # <-- starting point of vector
                0, -1, 0,  # <-- directions of vector
                color='red', alpha=.8, lw=2,
            )
            l1 = "-" + l1
            l2 = s2[c - 1] + l2
            l3 = "-" + l3
            c = c - 1

        if ((r <= 0 or c <= 0) or
                (r <= 0 or d <= 0) or
                (c <= 0 or d <= 0)):
            break

    # Adding missing parts to results
    ss1 = s1[0:r]
    ss2 = s2[0:c]
    ss3 = s3[0:d]
    l1 = ss1 + l1
    l2 = ss2 + l2
    l3 = ss3 + l3

    needed = max(len(l1), len(l2), len(l3)) - len(l1)
    l1 = '-' * needed + l1
    needed = max(len(l1), len(l2), len(l3)) - len(l2)
    l2 = '-' * needed + l2
    needed = max(len(l1), len(l2), len(l3)) - len(l3)
    l3 = '-' * needed + l3
    return l1, l2, l3




# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    match = 0
    mismatch = -15
    gap_t = -10
    extend = gap_t
    s1 = "acab"
    s2 = "dbad"
    s3 = "dabad"
    l1, l2, l3 = threeAligment()
    pairwise_aligment(s1, s2)
    score(l1, l2)
    pairwise_aligment(s1, s3)
    score(l1, l3)
    pairwise_aligment(s2, s3)
    score(l2, l3)


