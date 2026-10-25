#!/usr/bin/env python3
"""Exact fixed-witness certificate accompanying the finite-group paper.

Checks a Hamilton cycle in the nonisotropic Hermitian orthogonality graph
over F_64 using E=F_2[T]/(T^18+T^7+1) and h(x,y)=Tr_{E/F_64}(x*y^512).
There is no search, random choice, floating-point arithmetic, or dependency.
"""

from array import array
from math import isqrt

Q = 8
DEGREE = 18
POLYNOMIAL = (1 << 18) | (1 << 7) | 1
ORDER = (1 << DEGREE) - 1
PROJECTIVE_ORDER = Q**4 + Q**2 + 1
L = Q**2 - Q + 1
R = Q**2 + Q + 1
PRIMES = (3, 7, 19, 73)
CYCLE = (
    0, 23, 67, 25, 38, 54, 12, 36, 34, 45, 61, 9, 15, 70, 53, 66,
    10, 56, 44, 41, 55, 64, 2, 4, 57, 6, 40, 19, 30, 8, 68, 71,
    35, 24, 42, 21, 14, 32, 3, 65, 27, 11, 22, 46, 39, 28, 5, 33,
    17, 49, 18, 20, 50, 72, 47, 1, 48, 59, 7, 69, 37, 16, 60, 63,
)
SHIFTS = (
    23, 30, 22, 40, 1, 11, 13, 8, 26, 21, 4, 14, 18, 2, 16, 24,
    3, 6, 43, 20, 2, 24, 23, 40, 44, 4, 38, 4, 33, 31, 26, 10,
    16, 19, 35, 25, 5, 5, 6, 36, 5, 32, 10, 29, 17, 27, 37, 25,
    11, 4, 3, 6, 14, 54, 24, 16, 30, 2, 38, 37, 23, 32, 6, 2,
)


def multiply(a, b):
    """Multiplication in the binary polynomial quotient."""
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & (1 << DEGREE):
            a ^= POLYNOMIAL
    return result


def power(a, exponent):
    result = 1
    while exponent:
        if exponent & 1:
            result = multiply(result, a)
        a = multiply(a, a)
        exponent >>= 1
    return result


def polynomial_remainder(a, b):
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def polynomial_gcd(a, b):
    while b:
        a, b = b, polynomial_remainder(a, b)
    return a


def audit_field():
    """Rabin irreducibility criterion, followed by a full order certificate."""
    assert POLYNOMIAL.bit_length() - 1 == DEGREE
    assert power(2, 1 << DEGREE) == 2
    for prime in (2, 3):  # all prime divisors of 18
        assert polynomial_gcd(
            POLYNOMIAL, power(2, 1 << (DEGREE // prime)) ^ 2
        ) == 1
    remaining = ORDER
    for prime in PRIMES:
        assert prime >= 2 and all(prime % d for d in range(2, isqrt(prime) + 1))
        while remaining % prime == 0:
            remaining //= prime
        assert power(2, ORDER // prime) != 1
    assert remaining == 1 and power(2, ORDER) == 1


def exponential_table():
    table = array("I", [0]) * ORDER
    value = 1
    for i in range(ORDER):
        table[i] = value
        value <<= 1
        if value & (1 << DEGREE):
            value ^= POLYNOMIAL
    assert value == 1
    return table


def trace_power(table, exponent):
    """Trace of gamma^exponent from E to F_64."""
    return (
        table[exponent % ORDER]
        ^ table[(Q**2 * exponent) % ORDER]
        ^ table[(Q**4 * exponent) % ORDER]
    )


def main():
    if not __debug__:
        raise RuntimeError("Run this certificate with assertions enabled")
    audit_field()
    table = exponential_table()
    assert PROJECTIVE_ORDER == L * R
    assert ORDER == (Q**2 - 1) * PROJECTIVE_ORDER
    assert (R * (Q**3 + 1) * (Q - 1)) % ORDER == 0
    print("PASS: irreducible field polynomial and primitive element of order 262143")

    nonisotropic = [
        trace_power(table, (Q**3 + 1) * i) != 0
        for i in range(PROJECTIVE_ORDER)
    ]
    representatives = [i for i in range(R) if nonisotropic[i]]
    assert len(representatives) == Q**2
    assert len(CYCLE) == len(SHIFTS) == Q**2
    assert sorted(CYCLE) == representatives
    assert sum(nonisotropic) == Q**2 * L
    print("PASS: the 64 listed representatives are exactly the nonisotropic K-orbits")

    for i, a in enumerate(CYCLE):
        b = CYCLE[(i + 1) % len(CYCLE)]
        assert 0 <= SHIFTS[i] < L
        assert trace_power(table, a + Q**3 * (b + R * SHIFTS[i])) == 0
    total = sum(SHIFTS) % L
    assert total == 1
    print("PASS: all 64 base edges are orthogonal; total voltage is 1 modulo 57")

    visited = bytearray(PROJECTIVE_ORDER)
    lifted = []
    shift = 0
    for j in range(L):
        for i, a in enumerate(CYCLE):
            code = (a + R * (shift % L)) % PROJECTIVE_ORDER
            assert not visited[code] and nonisotropic[code]
            visited[code] = 1
            lifted.append(code)
            shift += SHIFTS[i]
        assert shift % L == (j + 1) % L
    assert len(lifted) == Q**2 * L
    assert all(bool(visited[i]) == nonisotropic[i] for i in range(PROJECTIVE_ORDER))
    for i, a in enumerate(lifted):
        b = lifted[(i + 1) % len(lifted)]
        assert trace_power(table, a + Q**3 * b) == 0
    print("PASS: all 3648 vertices and all 3648 edges of the closed lift audited")


if __name__ == "__main__":
    main()
