from fractions import Fraction as QQ
from math import gcd
from numpy import prod
from sympy import (
    totient as euler_phi,
    divisors,
    primefactors as prime_divisors,
    legendre_symbol,
)


def kronecker(n, p):
    """
    The Kronecker symbol.

    Parameters
    ----------
    n : int
        An integer.
    p : int
        A prime number.

    Returns
    -------
    int
        The Kronecker symbol (`n`/`p`) which is the Legendre symbol when `p` is
        odd, and if
        `p` == 2 is:
            0  if `n` == 0 mod 2
            1  if `n` == ± 1 mod 8
            -1 if `n` == ± 3 mod 8
    """
    if p % 2 == 1:
        return legendre_symbol(n, p)
    elif p == 2:
        if n % 2 == 0:
            return 0
        elif n % 8 in [1, 7]:
            return 1
        else:
            return -1
    else:
        raise ValueError("p is not a prime number")


def possible_r(N):
    """
    The smallest positive representatives for (ℤ/Nℤ)* modulo squares.

    Parameters
    ----------
    N : int
        A positive integer

    Returns
    -------
    list
        A list consisting of the smallest reps for (ℤ/Nℤ)* / ((ℤ/Nℤ)*)^2.

    Notes
    -----
    Do not use for anything non-trivial, inefficient algorithm is probably O(N)
    (exponential in the input length). Works by enumerating all square classes
    by writing out every element of the group (ℤ/Nℤ)*.
    """
    left = [x for x in range(1, N) if gcd(N, x) == 1]
    sqrs = [x**2 % N for x in left]
    classes = [sqrs]
    left = [x for x in left if not x in sqrs]
    done = False
    while not done:
        x = left[0]
        sq_cl = [(x * a) % N for a in sqrs]
        classes.append(sq_cl)
        left = [x for x in left if not x in sq_cl]
        if len(left) == 0:
            done = True
    reps = [min(cl) for cl in classes]
    return reps


def HJ_step(n, q):
    """
    Subprocess to computing the Hirzebruch--Jung contined fraction of a
    rational number.

    Parameters
    ----------
    n : int
        A positive integer.
    q : int
        A positive integer 1 < `q` < `n` coprime to `n`.

    Returns
    -------
    a : int
        The first entry in the HJ continued fraction of `n`/`q`
    new_n : int
        An updated n, so that the continued fraction of `n`/`q` is given by
        that of `new_n`/`new_q`.
    new_q : int
        An updated q, so that the continued fraction of `n`/`q` is given by
        that of `new_n`/`new_q`.
    """
    x = n % q
    a = ((n - x) // q) + 1
    new_n = q
    new_q = q - x
    return a, new_n, new_q


def HJ_continued_fraction(n, q):
    """Computes the Hirzebruch--Jung continued fraction of a rational number.

    Parameters
    ----------
    n : int
        A positive integer
    q : int
        A positive integer 1 <= `q` < `n`.

    Returns
    -------
    list
        A list consisting of the continued fraction of the rational number
        `n`/`q`.
    """
    if n <= 0 or q <= 0 or n <= q:
        raise ValueError("n,q should be positive and q < n")

    d = gcd(n, q)
    n = n // d
    q = q // d

    ret = []
    done = False
    while not done:
        if q == 1:
            done = True
            ret.append(n)
        else:
            step = HJ_step(n, q)
            if step[2] == 1:
                done = True
                ret += [step[0], step[1]]
            else:
                ret.append(step[0])
                n = step[1]
                q = step[2]
    return ret


def rat_part(x):
    """
    The rational part of a rational number. That is, its difference with its
    floor.

    Parameters
    ----------
    x : fractions.Fraction
        A rational number `x`

    Returns
    -------
    fractions.Fraction
        x - floor(x)
    """
    n = x.numerator
    d = x.denominator
    r = n % d
    f = (n - r) // d
    return x - f


def index_m(N):
    """
    The index of the modular curve X(N)

    Parameters
    ----------
    N : int
        A positive integer

    Returns
    -------
    int
        The index of the modular curve X(N), that is the order of the group
        SL_2(ℤ/Nℤ)/{±1}
    """
    pr_div = prime_divisors(N)
    prod_n = prod([p**2 - 1 for p in pr_div])
    prod_d = prod([p**2 for p in pr_div])
    return (N**3 * prod_n) // (2 * prod_d)


def r_0(N):
    """
    The class number of the quadratic order of discriminant -4*N**2.

    Parameters
    ----------
    N : int
        A positive integer

    Returns
    -------
    int
        The class number h(-4 * N**2) which Kani--Schanz call r_0(N)
    """
    ret = [QQ(int(p - kronecker(-4, p)), int(p)) for p in prime_divisors(N)]
    return int(N * prod(ret)) // 2


def r_1(N):
    """
    The class number of the quadratic order of discriminant -3*N**2

    Parameters
    ----------
    N : int
        A positive integer

    Returns
    -------
    int
        The class number h(-3 * N**2) which Kani--Schanz call r_1(N)
    """
    ret = [QQ(int(p - kronecker(-3, p)), int(p)) for p in prime_divisors(N)]
    return int(N * prod(ret)) // 3


def r_inf(N):
    """
    The number of quotient singularities on Z_(N,r) at the cusp.

    Parameters
    ----------
    N : int
        A positive integer

    Returns
    -------
    int
        The difference between the number of cusps on X_1(N) and 1/2*φ(N), this
        is Kani--Schanz's function r_∞(N)
    """
    ret = sum([euler_phi(d) * euler_phi(N // d) for d in divisors(N)]) 
    ret -= euler_phi(N)
    return ret // 2


def s_21(N, r):
    """
    The number of non-cuspidal quotient singularities of type (2,1) on Z_(N,r).
    Equivalently, the number of components of E_(2,1).

    Parameters
    ----------
    N : int
        The level of the HMS Z_(`N`,`r`).
    r : int
        The power of the congruence.

    Returns
    -------
    int
        The number of components of E_(2,1)

    Notes
    -----
    This is the Kani--Schanz function s_(0,1,r)
    """
    if N % 4 != 0:
        ret = r_0(N) // 2
    else:
        if r % 4 == 3:
            ret = r_0(N)
        else:
            ret = 0
    return ret


def s_31(N, r):
    """
    The number of non-cuspidal quotient singularities of type (3,1) on Z_(N,r).
    Equivalently, the number of components of E_(3,1).

    Parameters
    ----------
    N : int
        The level of the HMS Z_(`N`,`r`).
    r : int
        The power of the congruence.

    Returns
    -------
    int
        The number of components of E_(3,1)

    Notes
    -----
    This is the Kani--Schanz function s_(1,1,r)
    """
    if N % 3 != 0:
        return r_1(N) // 2
    else:
        if r % 3 == 1:
            return r_1(N)
        else:
            return 0


def s_32(N, r):
    """
    The number of non-cuspidal quotient singularities of type (3,2) on Z_(N,r).
    Equivalently, the number of disconnected (analytic topology) components of
    E_(3,2).

    Parameters
    ----------
    N : int
        The level of the HMS Z_(`N`,`r`).
    r : int
        The power of the congruence.

    Returns
    -------
    int
        The number of components of E_(3,2) (analytic topology).

    Notes
    -----
    This is the Kani--Schanz function s_(1,2,r)
    """
    if N % 3 != 0:
        ret = r_1(N) // 2
    else:
        if r % 3 == 2:
            ret = r_1(N)
        else:
            ret = 0
    return ret


def bbL(n, q):
    """
    The length of the continued fraction of n/[q] where 0 < [q]  < n.

    Parameters
    ----------
    n : int
        Positive integer
    q : int
        Positive integer

    Returns
    -------
    int
        The length of the continued fraction expansion of `n`/[q] where
        0 < [q] < `n` is the same as `q` modulo `n`.

    Notes
    -----
    This is the Kani--Schanz function 𝕃_(n,q).
    """
    d = gcd(n, q)
    n = n // d
    q = q // d
    q = q % n
    return len(HJ_continued_fraction(n, q))


def bbL_infr(N, r):
    """
    A certain sum of lengths of continued fraction expansions of cyclic
    quotient singularities occuring on Z_(N,r).

    Parameters
    ----------
    N : int
        The level of the HMS.
    r : int
        A positive integer coprime to `N`, the power of the congruence.

    Returns
    -------
    int

    Notes
    -----
    This is the function Kani--Schanz call 𝕃_(∞,r).
    """
    ret = 0
    for d in divisors(N):
        if d != N:
            term = sum(
                [
                    bbL(N, r * n**2 * d)
                    for n in range(1, (N // d))
                    if gcd(n, N // d) == 1
                ]
            )
            ret += euler_phi(d) * term
    return ret // 2


def bbL_r(N, r):
    """
    A specific sum involving r_0, r_1, s_31, and 𝕃_(∞,r).

    Parameters
    ----------
    N : int
        The level of the HMS.
    r : int
        A positive integer coprime to `N`, the power of the congruence.

    Returns
    -------
    int

    Notes
    -----
    This is the function Kani--Schanz call 𝕃_(r)
    """
    return r_0(N) + 2 * r_1(N) - s_31(N, r) + bbL_infr(N, r)


def sawtooth(x):
    """
    Computes the sawtooth function of a rational number which appears in a
    Dedekind sum.

    Parameters
    ----------
    x : fractions.Fraction
        A rational number

    Returns
    -------
    fractions.Fraction
        The sawtooth function, typically denoted ((x)) = `x` - floor(`x`) - 1/2
    """
    if x.denominator == 1:
        return 0
    else:
        return rat_part(x) - QQ(1, 2)


def bbS(q, n):
    """
    The typical Dedekind sum over products of sawtooth functions with n in the
    denominator.

    Parameters
    ----------
    q : int
        A positive integer
    n : int
        An integer coprime to `q`

    Returns
    -------
    fractions.Fraction
        The Dedekind sum Σ_(i=1)^(n-1) ((k/n)) * ((k*q/n))

    Notes
    -----
    Kani--Schanz call this function 𝕊(q,n). It may also be called D(q, 1 ;  
    elsewhere (e.g., wikipedia as of 20/12/2023)
    """
    ret = sum(
        [   
            sawtooth(QQ(k, n)) * sawtooth(QQ(k * q, n)) 
            for k in range(1, n)
        ]   
    )
    return ret


def bbS_infr(N, r):
    """
    A certain sum over Dedekind sums which is indexed by cusp quotient
    singularities of Z_(N,r).

    Parameters
    ----------
    N : int
        The level of the HMS.
    r : int
        A positive integer coprime to `N`, the power of the congruence.

    Returns
    -------
    int

    Notes
    -----
    Kani--Schanz call this function 𝕊_(∞,r).
    """
    ret = 0
    for d in [int(d) for d in divisors(N)]:
        if d != N:
            term = sum(
                [
                    bbS(r * n**2, N // d)
                    for n in range(1, (N // d) + 1)
                    if gcd(n, N // d) == 1
                ]
            )
            ret += euler_phi(d) * term
    return QQ(ret, 2)


def R_infr(N, r):
    """
    A certain sum of 𝕊_(∞,r), 𝕃_(∞,r), and r_∞(N).

    Parameters
    ----------
    N : int
        The level of the HMS.
    r : int
        A positive integer coprime to `N`, the power of the congruence.

    Returns
    -------
    int

    Notes
    -----
    This is the function Kani--Schanz call R_(∞,r).
    """
    return int(12 * bbS_infr(N, r)) + bbL_infr(N, r) + r_inf(N)


def bbS_r(N, r):
    """A certain sum of s_(3,1), r_1(N), 𝕊_(∞,r).

    Parameters
    ----------
    N : int
        The level of the HMS.
    r : int
        A positive integer coprime to `N`, the power of the congruence.

    Returns
    -------
    int

    Notes
    -----
    This is the function Kani--Schanz call 𝕊_(r)
    """
    return (2 * s_31(N, r) - r_1(N) + 18 * bbS_infr(N, r)) // 18


def p_adic_val(n, p):
    """
    The p-adic valuation of an integer.

    Parameters
    ----------
    n : int
        A positive integer
    p : int
        A prime number

    Returns
    -------
    k : int
        The greatest integer such that `p`^`k` | `n`

    Notes
    -----
    This is a very inefficient implementation, do not use for anything
    non-trivial
    """
    done = False
    k = 0
    while not done:
        if n % p != 0:
            done = True
        else:
            n = n // p
            k += 1
    return k


def rho(N, r):
    """
    The invariant ρ(N, r) = #{x ∈ ℤ/Nℤ : x^2 = r} when gcd(N,r) == 1 and 0 o/w

    Parameters
    ----------
    N : int
        A positive integer
    r : int
        A positive integer

    Returns
    -------
    int
        The order of {x ∈ ℤ/Nℤ : x^2 = r} when gcd(N,r) == 1 and 0 otherwise.

    Notes
    -----
    Inefficient implementation due to 2-adic valuation being bad, shouldn't
    matter for anything reasonably small (we'll only need N < 6000 at worst).
    """
    if gcd(N, r) != 1:
        return 0
    k = p_adic_val(N, 2)
    M = N // (2**k)
    # compute the odd part of the formula
    if M == 1:
        ret = 1
    else:
        ret = prod([1 + int(kronecker(r, p)) for p in prime_divisors(M)])
    # compute the 2-part of the formula
    if k == 2:
        if r % 4 == 1:
            ret = 2 * ret
        else:
            ret = 0
    elif k >= 3:
        if r % 8 == 1:
            ret = 4 * ret
        else:
            ret = 0
    return ret


def rho_2(N, r):
    """
    The modified version of rho(N,r) that is, ρ_2(N, r) = ρ(N/2, r) if
    N = 2 (4) and 0 otherwise.

    Parameters
    ----------
    N : int
        A positive integer
    r : int
        A positive integer

    Returns
    -------
    int
        An integer ρ(N/2, r) if N = 2 (4) and 0 otherwise.
    """
    if N % 4 == 2:
        ret = rho(N // 2, r)
    else:
        ret = 0
    return ret


def rho_3(N, r):
    """
    The modified version of rho(N,r) that is, ρ_3(N, r) = ρ(N/3, r) if
    N = 3,6 (9) and 0 otherwise.

    Parameters
    ----------
    N : int
        A positive integer
    r : int
        A positive integer

    Returns
    -------
    int
        ρ_3(N, r) = ρ(N/3, r) if N = 3,6 (9) and 0 otherwise.
    """
    if N % 9 == 3 or N % 9 == 6:
        ret = rho(N // 3, r)
    else:
        ret = 0
    return ret
