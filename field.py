"""
Step 1: finite fields GF(31) and GF(16).

Every operation (add, sub, mul, inv) is a lookup in a precomputed table,
as the project guide asks ("tables d'arithmétique").
"""


class Field:
    """A finite field with q elements, numbered 0..q-1.

    Subclasses only fill the tables in _build(); the rest is shared,
    so GF(16) and GF(31) can be swapped anywhere in the code.
    """

    def __init__(self, q):
        self.q = q
        self.ADD = [[0] * q for _ in range(q)]
        self.SUB = [[0] * q for _ in range(q)]
        self.MUL = [[0] * q for _ in range(q)]
        self.INV = [0] * q          # INV[0] stays 0 (no inverse)
        self._build()

    def add(self, a, b): return self.ADD[a][b]
    def sub(self, a, b): return self.SUB[a][b]
    def mul(self, a, b): return self.MUL[a][b]
    def neg(self, a):    return self.SUB[0][a]

    def inv(self, a):
        if a == 0:
            raise ZeroDivisionError("0 has no inverse")
        return self.INV[a]

    def div(self, a, b): return self.mul(a, self.inv(b))

    def __repr__(self):
        return f"GF({self.q})"


class GF31(Field):
    """Prime field: plain arithmetic modulo 31 (like a clock with 31 hours)."""

    def __init__(self):
        super().__init__(31)

    def _build(self):
        p = 31
        for a in range(p):
            for b in range(p):
                self.ADD[a][b] = (a + b) % p
                self.SUB[a][b] = (a - b) % p
                self.MUL[a][b] = (a * b) % p
        for a in range(1, p):
            self.INV[a] = pow(a, p - 2, p)   # Fermat: a^(p-2) = a^-1 mod p


class GF16(Field):
    """GF(2^4): each element is a 4-bit number = polynomial of degree < 4 over {0,1}.

    Example: 11 = 0b1011 = x^3 + x + 1.
    Addition = XOR (coefficients are added mod 2).
    Multiplication = polynomial product, reduced modulo x^4 + x + 1 (0b10011).
    """

    MODULUS = 0b10011  # x^4 + x + 1

    def __init__(self):
        super().__init__(16)

    @staticmethod
    def _slow_mul(a, b):
        """Schoolbook polynomial multiplication with reduction (used to build tables)."""
        result = 0
        while b:
            if b & 1:              # if the current bit of b is 1, add (XOR) a
                result ^= a
            b >>= 1
            a <<= 1                # a = a * x
            if a & 0b10000:        # degree reached 4: replace x^4 by x + 1
                a ^= GF16.MODULUS
        return result

    def _build(self):
        for a in range(16):
            for b in range(16):
                self.ADD[a][b] = a ^ b
                self.SUB[a][b] = a ^ b       # in characteristic 2, -b = b
                self.MUL[a][b] = self._slow_mul(a, b)
        for a in range(1, 16):
            self.INV[a] = next(b for b in range(1, 16) if self.MUL[a][b] == 1)
