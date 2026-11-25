import random
import pytest
from field import GF16, GF31
from linalg import (solve, mat_vec, mat_mul, transpose, random_matrix,
                    random_vector, quad_form, upper, dot)

FIELDS = [GF16(), GF31()]


@pytest.mark.parametrize("F", FIELDS, ids=repr)
def test_solve_recovers_solution(F):
    """Pick a random y, compute t = M y, then solve(M, t) must give y back."""
    rng = random.Random(1)
    solved = 0
    for _ in range(300):
        n = rng.randint(1, 10)
        M = random_matrix(F, n, n, rng)
        y = random_vector(F, n, rng)
        t = mat_vec(F, M, y)
        sol = solve(F, M, t)
        if sol is not None:          # M invertible -> unique solution
            assert sol == y
            solved += 1
    assert solved > 250              # most random matrices are invertible


@pytest.mark.parametrize("F", FIELDS, ids=repr)
def test_singular_matrix_returns_none(F):
    rng = random.Random(2)
    M = random_matrix(F, 4, 4, rng)
    M[3] = M[0][:]                   # two identical rows -> singular
    assert solve(F, M, [1, 2, 3, 4]) is None
    assert solve(F, [[0, 0], [0, 0]], [0, 0]) is None


@pytest.mark.parametrize("F", FIELDS, ids=repr)
def test_mat_mul_and_transpose(F):
    rng = random.Random(3)
    A = random_matrix(F, 3, 4, rng)
    B = random_matrix(F, 4, 2, rng)
    x = random_vector(F, 2, rng)
    # (A B) x == A (B x)
    assert mat_vec(F, mat_mul(F, A, B), x) == mat_vec(F, A, mat_vec(F, B, x))
    assert transpose(transpose(A)) == A


@pytest.mark.parametrize("F", FIELDS, ids=repr)
def test_quad_form_and_upper(F):
    """x^T M x computed directly == quad_form(upper(M), x)."""
    rng = random.Random(4)
    for _ in range(50):
        n = rng.randint(1, 8)
        M = random_matrix(F, n, n, rng)
        x = random_vector(F, n, rng)
        direct = dot(F, x, mat_vec(F, M, x))
        assert quad_form(F, upper(F, M), x) == direct