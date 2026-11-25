import itertools
import pytest
from field import GF16, GF31

FIELDS = [GF16(), GF31()]


def test_gf16_reference_value():
    # From the guide: 2 * 8 = x * x^3 = x^4 = x + 1 = 3
    assert GF16().mul(2, 8) == 3


@pytest.mark.parametrize("F", FIELDS, ids=repr)
def test_field_axioms(F):
    q = range(F.q)
    for a, b in itertools.product(q, q):
        assert F.add(a, b) == F.add(b, a)
        assert F.mul(a, b) == F.mul(b, a)
        assert F.add(F.sub(a, b), b) == a
    for a, b, c in itertools.product(q, q, q):
        assert F.add(F.add(a, b), c) == F.add(a, F.add(b, c))
        assert F.mul(F.mul(a, b), c) == F.mul(a, F.mul(b, c))
        assert F.mul(a, F.add(b, c)) == F.add(F.mul(a, b), F.mul(a, c))
    for a in range(1, F.q):
        assert F.mul(a, F.inv(a)) == 1


@pytest.mark.parametrize("F", FIELDS, ids=repr)
def test_zero_has_no_inverse(F):
    with pytest.raises(ZeroDivisionError):
        F.inv(0)
