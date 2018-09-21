from Bio.SubsMat import MatrixInfo
import numpy

# Ahmad Fauzi Wibowo (16/394073/PA/17164)
# Program implementasi algoritma Needleman-Wunsch
# Menggunakan blosum62 sebagai substitution matrix
blosum = MatrixInfo.blosum62


# Function utama algoritma Needleman-Wunsch
def needleman(s1, s2, blosum, gap_penalty):
    rows, col = len(s2) + 1, len(s1) + 1
    M = numpy.zeros((rows, col), int)  # Matriks skor
    P = numpy.zeros((rows, col), int)  # Matriks pointer traceback


    for i in range(rows):
        M[i][0] = i * gap_penalty
        P[i][0] = 3
    for j in range(col):
        M[0][j] = j * gap_penalty
        P[0][j] = 2



     
    for i in range(1, rows):
        for j in range(1, col):
            k1, k2 = s1[j-1], s2[i-1]
            if (k1, k2) in blosum:
                d = M[i-1][j-1] + blosum[(k1, k2)]
            else:
                d = M[i-1][j-1] + blosum[(k2, k1)]
            u = M[i-1][j] + gap_penalty
            l = M[i][j-1] + gap_penalty
            M[i][j] = max(d, u, l)

            # angka 1 = diagonal, angka 2 = kiri, angka 3 = atas
            if M[i][j] == d: P[i][j] = 1
            elif M[i][j] == l: P[i][j] = 2
            elif M[i][j] == u: P[i][j] = 3
    
    print(M)
    print(P)


    # Traceback
    i, j = rows - 1, col - 1
    al_s1 = al_s2 = ""
    al_score = M[i][j]


    while i > 0 and j > 0:
        if P[i][j] == 1: # arah diagonal
            al_s1 = s1[j-1] + al_s1
            al_s2 = s2[i-1] + al_s2
            i -= 1
            j -= 1
        elif P[i][j] == 2: # arah kiri
            al_s1 = s1[j-1] + al_s1
            al_s2 = "-" + al_s2
            j -= 1
        elif P[i][j] == 3: # arah atas
            al_s1 = "-" + al_s1
            al_s2 = s2[i-1] + al_s2
            i -= 1

    while i > 0:
        al_s1 = "-" + al_s1
        al_s2 = s2[i-1] + al_s2
        i -= 1

    while j > 0:
        al_s1 = s1[j-1] + al_s1
        al_s2 = "-" + al_s2
        j -= 1


    al_symbol = ''
    for i in range(len(al_s1)):
        al_symbol += "|" if al_s1[i] == al_s2[i] else " "


    return al_score, al_s1, al_s2, al_symbol

if __name__ == '__main__':
    s1 = "MNALQM"
    s2 = "NALMSQA"

    gap_penalty = -6

    score, al_s1, al_s2, al_sym = needleman(s1, s2, blosum, gap_penalty)

    print al_s1
    print al_sym
    print al_s2
    print "Score =", score