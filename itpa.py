itpa = 'MCSRNIKISVVLFLVLIPIFAALPHNHNLSKRSNFFDLECKGIFNKTMFFRLDRICEDCYQLFRETSIHRLCKQECFGSPFFNACIEALQLHEEMDKYNEWRDTLGRK'

cys_coord = []
for coord, aminoAcid in enumerate(itpa):
    if aminoAcid == 'C':
        cys_coord.append(coord)
cys_coord.append(len(itpa))

print(cys_coord) # [1, 39, 55, 58, 71, 75, 84]

difference = []
n = 0
while n < len(cys_coord) - 1:
    x1 = n + 1
    x0 = n
    difference.append(cys_coord[x1] - cys_coord[x0])
    n += 1

print(difference) # [38, 16, 3, 13, 4, 9]

motif = []
m = 0
while m < len(cys_coord) - 1:
    pre_motif = ''
    y1 = m + 1
    y0 = m
    for i in range(cys_coord[y0] + 1, cys_coord[y1]):
        pre_motif += itpa[i]
    motif.append(pre_motif)
    m += 1

print(motif)

length_motif = []
for i in motif:
    length_motif.append(len(i))

print(length_motif)

"""
                37                             15        2      12       3     8
 |                                     |               |  |            |   |        |
MCSRNIKISVVLFLVLIPIFAALPHNHNLSKRSNFFDLECKGIFNKTMFFRLDRICEDCYQLFRETSIHRLCKQECFGSPFFNACIEALQLHEEMDKYNEWRDTLGRK
"""