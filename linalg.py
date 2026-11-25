"""
Step 2: linear algebra over a finite field F (GF16 or GF31).

Matrices are plain Python lists of rows: M[i][j] = row i, column j.
Vectors are plain lists. Every operation goes through F.add / F.sub / F.mul / F.inv,
so the same code works for both fields.
"""


# ---------------------------------------------------------------- basic operations

def transpose(A):
    """Swap rows and columns."""
    return [list(col) for col in zip(*A)]


def dot(F, u, v):
    """Dot product u . v = u1*v1 + u2*v2 + ... (computed in the field)."""
    s = 0
    for a, b in zip(u, v):
        s = F.add(s, F.mul(a, b))
    return s


def mat_vec(F, M, x):
    """Matrix times vector: each output entry is (row of M) . x"""
    return [dot(F, row, x) for row in M]


def vec_mat(F, x, M):
    """Row vector times matrix: x^T M."""
    return mat_vec(F, transpose(M), x)


def mat_mul(F, A, B):
    """Matrix times matrix: entry (i, j) = (row i of A) . (column j of B)."""
    Bt = transpose(B)
    return [[dot(F, row, col) for col in Bt] for row in A]


def mat_add(F, A, B):
    return [[F.add(a, b) for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def vec_add(F, u, v):
    return [F.add(a, b) for a, b in zip(u, v)]


def vec_sub(F, u, v):
    return [F.sub(a, b) for a, b in zip(u, v)]


# ---------------------------------------------------------------- random matrices

def random_matrix(F, rows, cols, rng):
    return [[rng.randrange(F.q) for _ in range(cols)] for _ in range(rows)]


def random_upper(F, n, rng):
    """Random upper-triangular n x n matrix (zeros below the diagonal)."""
    return [[rng.randrange(F.q) if j >= i else 0 for j in range(n)] for i in range(n)]


def random_vector(F, n, rng):
    return [rng.randrange(F.q) for _ in range(n)]


# ---------------------------------------------------------------- Gaussian elimination

def solve(F, M, t):
    """Solve M y = t for a square matrix M over F.

    Returns the solution y, or None if M is singular (no unique solution).
    In UOV, None means "bad luck with the vinegar, pick new vinegar and retry".
    """
    n = len(M)
    # Augmented matrix [M | t]: we work on a copy, never on the caller's data.
    A = [row[:] + [t[i]] for i, row in enumerate(M)]

    for c in range(n):                       # for each column c
        # 1. Find a pivot: a row (from c downwards) with a nonzero entry in column c.
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            return None                      # whole column is zero -> M singular
        A[c], A[p] = A[p], A[c]              # move the pivot row into position c

        # 2. Normalise the pivot row so the pivot becomes 1 (divide by the pivot).
        iv = F.inv(A[c][c])
        A[c] = [F.mul(iv, x) for x in A[c]]

        # 3. Eliminate column c from every other row: row_r = row_r - f * row_c.
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [F.sub(x, F.mul(f, y)) for x, y in zip(A[r], A[c])]

    # Now the left part is the identity, so the last column is the solution.
    return [A[i][n] for i in range(n)]


# ---------------------------------------------------------------- quadratic forms

def quad_form(F, A, x):
    """Evaluate x^T A x = sum over i <= j of A[i][j] * x_i * x_j.

    This is how a public equation P_k is evaluated during verification.
    """
    n = len(x)
    s = 0
    for i in range(n):
        if x[i] == 0:
            continue
        for j in range(i, n):
            if A[i][j] and x[j]:
                s = F.add(s, F.mul(A[i][j], F.mul(x[i], x[j])))
    return s


def upper(F, M):
    """Fold a square matrix into upper-triangular form with the same quadratic form.

    x^T M x only depends on M[i][j] + M[j][i] for i != j, so we put that sum
    above the diagonal and 0 below. Diagonal entries stay as they are.
    """
    n = len(M)
    U = [[0] * n for _ in range(n)]
    for i in range(n):
        U[i][i] = M[i][i]
        for j in range(i + 1, n):
            U[i][j] = F.add(M[i][j], M[j][i])
    return U