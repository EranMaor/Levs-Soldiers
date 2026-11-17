import numpy as np
import itertools

# weight vectors of Table 3 (exact fractions)
WEIGHTS = {
    "w41": [1/15, 2/15, 4/15, 8/15],   # L=4, |Theta|=256, rate 8.00
    "w42": [1/8, 1/8, 3/8, 3/8],       # L=4, |Theta|=100, rate 6.64
    "w21": [1/3, 2/3],                 # L=2, |Theta|=16,  rate 4.00
}


def get_weights(name):
    return np.array(WEIGHTS[name], dtype=float)


def create_alphabet(w):
    # Returns Theta as an m x 4 matrix
    # Different base combinations can produce the same vector (e.g. w42), those are merged,
    # so m is the number of distinct vectors.
    w = np.asarray(w, dtype=float)
    # Tolerance for small floating point errors
    if not np.isclose(w.sum(), 1.0, atol=1e-9):
        raise ValueError(f"weights must sum to 1 (got {w.sum():.12f})")
    L = len(w)

    bases = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    base_combos = itertools.product(bases, repeat=L)

    Theta = []
    seen = set()
    for e in base_combos:
        theta = np.zeros(4)
        for l in range(L):
            theta += w[l] * e[l]
        key = tuple(np.round(theta, 9))
        if key not in seen:
            seen.add(key)
            Theta.append(theta)

    return np.array(Theta)

# Debug print
if __name__ == "__main__":
    for name in WEIGHTS:
        alphabet = create_alphabet(get_weights(name))
        print(f"{name}: |Theta| = {len(alphabet)}, rate = {np.log2(len(alphabet)):.2f}")

    alphabet = create_alphabet(get_weights("w21"))
    for row in alphabet:
        for element in row:
            print(f"{element:.3f}", end='\t')
        print()
