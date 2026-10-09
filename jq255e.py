"""
jq255e.py -- Reference implementation of the jq255e Jacobi-quartic elliptic curve.

Paper
-----
A Genus-1 Reformulation of an alpha-beta Modified Edwards-Type Curve
for Secure Telemedicine Data Transmission

Authors: P. Adeola, A. T. Akinwale, B. B. Olofin, A. E. Babalola
         Department of Computing, Anchor University, Lagos, Nigeria

Curve parameters
----------------
Equation : E^2 = 8*u^4 + 1   (Jacobi quartic)
Field    : F_q, q = 2^255 - 18651
Order    : #E(F_q) = 2*r
Subgroup : r = 2^254 - 131528281291764213006042413802501683931  (prime)
Cofactor : h = 2
Generator: G = (E=3, u=1)  -- full curve, order 2r
           H = 2G           -- prime-order subgroup generator, order r

Source   : https://doubleodd.group/params-do255e.html
Reference: T. Pornin, "Double-Odd Jacobi Quartic," IACR ePrint 2022/1052.

Addition formulas (manuscript Equations 23-25)
----------------------------------------------
Given P1 = (E1, u1) and P2 = (E2, u2):

    D  = 1 - 8*u1^2*u2^2
    u3 = (u1*E2 + u2*E1) / D
    E3 = ((1 + 8*u1^2*u2^2)*E1*E2 + 16*u1*u2*(u1^2 + u2^2)) / D^2

Identity  : (E=1, u=0)
Inverse   : -(E, u) = (E, -u)

IMPORTANT
---------
This is a REFERENCE IMPLEMENTATION intended for research validation only.

- It is written in pure Python using affine coordinates.
- The affine addition formulas are NOT complete: the exceptional case D = 0
  is not handled and will raise ValueError. A complete, exception-free group
  law requires projective or extended Jacobi-quartic coordinates, which is
  left for future work.
- The implementation is NOT constant-time and is NOT resistant to
  side-channel attacks.
- It must NOT be used in production or clinical deployment.
"""

# ---------------------------------------------------------------------------
# Field and subgroup parameters
# ---------------------------------------------------------------------------

q = 2**255 - 18651
r = 2**254 - 131528281291764213006042413802501683931


def _modinv(a, m=q):
    """Modular multiplicative inverse via Fermat's little theorem (m prime)."""
    return pow(a, m - 2, m)


# ---------------------------------------------------------------------------
# Point class
# ---------------------------------------------------------------------------

class JQ255ePoint:
    """
    Affine point (E, u) on E^2 = 8*u^4 + 1 over F_q.

    All arithmetic is performed modulo q.  Private scalars are reduced
    modulo r (the prime subgroup order) before scalar multiplication.
    """

    def __init__(self, E, u):
        self.E = int(E) % q
        self.u = int(u) % q

    # ------------------------------------------------------------------
    # Predicates
    # ------------------------------------------------------------------

    def is_identity(self):
        """Return True if this point is the identity element (E=1, u=0)."""
        return self.u == 0 and self.E == 1

    def on_curve(self):
        """Return True if (E, u) satisfies E^2 = 8*u^4 + 1 mod q."""
        lhs = (self.E * self.E) % q
        rhs = (8 * pow(self.u, 4, q) + 1) % q
        return lhs == rhs

    def in_prime_subgroup(self):
        """Return True if r*self == identity (prime-subgroup membership)."""
        return (self * r).is_identity()

    # ------------------------------------------------------------------
    # Group operations
    # ------------------------------------------------------------------

    def __neg__(self):
        return JQ255ePoint(self.E, (-self.u) % q)

    def __eq__(self, other):
        if not isinstance(other, JQ255ePoint):
            return NotImplemented
        return self.E == other.E and self.u == other.u

    def __hash__(self):
        return hash((self.E, self.u))

    def __add__(self, other):
        if not isinstance(other, JQ255ePoint):
            return NotImplemented

        # Identity shortcuts
        if self.is_identity():
            return JQ255ePoint(other.E, other.u)
        if other.is_identity():
            return JQ255ePoint(self.E, self.u)

        # Inverse shortcut: P + (-P) = identity
        if self == -other:
            return JQ255ePoint(1, 0)

        E1, u1 = self.E, self.u
        E2, u2 = other.E, other.u

        u1sq = (u1 * u1) % q
        u2sq = (u2 * u2) % q
        A    = (8 * u1sq % q * u2sq) % q   # 8*u1^2*u2^2

        D = (1 - A) % q
        if D == 0:
            raise ValueError(
                "Exceptional case in affine Jacobi-quartic addition (D = 0). "
                "Projective coordinates are required to handle this case."
            )

        D_inv  = _modinv(D)
        D_inv2 = (D_inv * D_inv) % q

        # Equation (24): u3 = (u1*E2 + u2*E1) / D
        u3 = ((u1 * E2 + u2 * E1) % q * D_inv) % q

        # Equation (25): E3 = ((1+A)*E1*E2 + 16*u1*u2*(u1^2+u2^2)) / D^2
        term1 = ((1 + A) % q * E1 % q * E2) % q
        term2 = (16 * u1 % q * u2 % q * ((u1sq + u2sq) % q)) % q
        E3    = ((term1 + term2) % q * D_inv2) % q

        return JQ255ePoint(E3, u3)

    def __mul__(self, k):
        """
        Scalar multiplication k*P using binary double-and-add.

        The scalar k is reduced modulo r before use.  This method requires
        approximately 254 point doublings and on average 127 additional point
        additions for a 254-bit scalar.
        """
        k = int(k) % r
        if k == 0:
            return JQ255ePoint(1, 0)

        result = JQ255ePoint(1, 0)   # identity
        addend = JQ255ePoint(self.E, self.u)

        while k:
            if k & 1:
                result = result + addend
            addend = addend + addend
            k >>= 1

        return result

    def __rmul__(self, k):
        return self.__mul__(k)

    def encode(self):
        """
        Canonical 64-byte encoding of the point: E || u, each 32 bytes
        little-endian.  Used as input to the BLAKE2s key-derivation step.
        """
        return (self.E.to_bytes(32, 'little') +
                self.u.to_bytes(32, 'little'))

    def __repr__(self):
        return f"JQ255ePoint(E={self.E}, u={self.u})"


# ---------------------------------------------------------------------------
# Published constants
# ---------------------------------------------------------------------------

IDENTITY = JQ255ePoint(1, 0)

# Conventional generator G (full curve, order 2r).
# Coordinates from the published jq255e specification.
G = JQ255ePoint(E=3, u=1)

# Prime-order subgroup generator H = 2G (order r).
H = G + G


# ---------------------------------------------------------------------------
# Public-key validation
# ---------------------------------------------------------------------------

def validate_public_key(point):
    """
    Validate a received public key or ephemeral point.

    Checks (in order):
      1. The point satisfies the curve equation E^2 = 8*u^4 + 1 mod q.
      2. The point is not the identity element.
      3. The point belongs to the prime-order subgroup (r*P == identity).

    Returns True on success; raises ValueError with a descriptive message
    on failure.  Note: step 3 requires one scalar multiplication and is
    the most expensive check.
    """
    if not point.on_curve():
        raise ValueError("Received point does not satisfy the curve equation.")
    if point.is_identity():
        raise ValueError("Received point is the identity element.")
    if not point.in_prime_subgroup():
        raise ValueError(
            "Received point does not belong to the prime-order subgroup."
        )
    return True
