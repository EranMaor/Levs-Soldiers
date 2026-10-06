import numpy as np
import itertools

def create_alphabet(m, w, L):
    Theta = np.zeros((4, m))

    bases = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    base_combos = list((itertools.product(bases, repeat=L)))
    
    if sum(w) == 1.0:
        for symbol in range(m):
            e = base_combos[symbol]
            for l in range(L):
                    Theta[:, symbol] += w[l] * e[l]

    return Theta

w_2_1 = np.array([0.333, 0.667])
alphabet = create_alphabet(16, w_2_1, 2)
for row in alphabet:
    for element in row:
        print(element, end='\t')
    print()