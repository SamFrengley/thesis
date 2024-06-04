from fractions import Fraction as QQ
from math import gcd
from sympy import (
    totient as euler_phi, 
    divisors, 
    primefactors as prime_divisors
)
from numpy import prod
from .arithmetic import *


def RH_genus(mu, e_2, e_3, e_inf):
    """
    Compute the genus of a modular curve using the Riemann--Hurwitz formula.

    Parameters
    ----------
    mu : int
        The index of the subgroup H < GL_2(ℤ/Nℤ).
    e_2 : int
        The number of ell. points of order 2.
    e_3 : int
        The number of ell. points of order 3.
    e_inf : int
        The number of cusps.

    Returns
    -------
    int
       The genus of the modular curve with these parameters.
    """
    g_times_12 = 12 + (mu - 3 * e_2 - 4 * e_3 - 6 * e_inf)
    g, r = divmod(g_times_12 % 12)
    if r != 0:
        raise ValueError("This cannot possibly be a modular curve")
    return g


class X_1:
    """
    Class of the modular curve X_1(N)

    Parameters
    ----------
    N : int
        The level of the modular curve

    Attributes
    ----------
    N : int
        The level of the modular curve
    """

    def __init__(self, N):
        self.N = N

    def genus(self):
        """
        The genus of the modular curve X_1(N).

        Returns
        -------
        int
           The genus of the modular curve X_1(N).
        """
        N = self.N
        if N <= 5:
            ret = 0
        else:
            ret = 1
            pr_div = prime_divisors(N)
            prod_n = prod([p**2 - 1 for p in pr_div])
            prod_d = prod([p**2 for p in pr_div])
            ret += QQ(N**2, 24) * QQ(prod_n, prod_d)
            sum_term = sum(
                [
                    euler_phi(v) * euler_phi(N // v) 
                    for v in divisors(N)
                ]
            )
            ret += QQ(-int(sum_term), 4)
        return int(ret)


def sq_free_divisors(N):
    """
    Computes the sq. free divisors of an integer.

    Parameters
    ----------
    N : int
        An integer.

    Returns
    -------
    list
        A list of those squarefree integers d which divide `N`.
    """
    dd = divisors(N)
    ret = []
    for d in dd:
        pp = prime_divisors(d)
        sq_f = True
        i = 0
        while sq_f and i < len(pp):
            if d % (pp[i] ** 2) == 0:
                sq_f = False
            i += 1
        if sq_f:
            ret.append(d)
    return [int(x) for x in ret]


def get_rat_X0_plus():
    """
    The integers such that X_0^+(N) has genus 0.

    Returns
    -------
    list
        A list of integers such that the Fricke quotient of X_0(N) is a
        rational curve.

    Note
    ----
    Implemented as a hardcoded list.
    """
    ret_list = [
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        23,
        24,
        25,
        26,
        27,
        29,
        31,
        32,
        35,
        36,
        39,
        41,
        47,
        49,
        50,
        59,
        71
    ]
    return ret_list


def fricke_sort(N, cusps, widths=[]):
    """
    Make the Fricke involution on the cusps of X_0(N) reverses the list.
    Helper for class X_0.

    Parameters
    ----------
    N : int
        The level of the modular curve X_0(N)
    cusps : list
        A list of tuples (d,x) where d|N and x runs over reps for
        ℤ/gcd(N,N/d)ℤ.
    widths=[] : list
        If provided should be a list of length len(`cusps`) of the
        corresponding cusps widths.
    """
    new1 = []
    new2 = []
    fixed = []
    for c in cusps:
        if not c in new1 + new2:
            d, x = c
            delta = gcd(d, N // d)
            new_d = N // d
            new_x = -x % delta
            tau_c = (new_d, new_x)
            if c != tau_c:
                new1.append(c)
                new2.append(tau_c)
            else:
                fixed.append(c)
    new2.reverse()
    assert len(fixed) <= 1
    new_cusps = new1 + fixed + new2
    new_wd = []
    for c in new_cusps:
        i = cusps.index(c)
        new_wd.append(widths[i])
    cusps = new_cusps
    widths = new_wd
    return cusps, widths


class X_0:
    """
    Class of the modular curve X_0(N).

    Parameters
    ----------
    N : int
        The level of the modular curve.

    Attributes
    ----------
    N : int
        The level of the modular curve.
    mu : int
        The index of the modular curve.
    e_2 : int
        The number of elliptic points of order 2.
    e_3 : int
        The number of elliptic points of order 3.
    e_inf : int
        The number of cusps.
    cusps : list
        A list of tuples (d,x) where d|N and x runs over reps for
        ℤ/gcd(N,N/d)ℤ. The order is chosen so that the Fricke involution
        acts as cusps.reverse().
    widths : list
        A list of cusp widths corresponding directly to the cusps.
    """

    def __init__(self, N):
        self.N = N
        self.mu, self.e_2, self.e_3, self.e_inf = self.init_ramification()
        self.cusps, self.widths = self.make_cusps()

    def init_ramification(self):
        """
        Gets data about the map X_0(N) -> X(1)

        Returns
        -------
        mu : int
            The index of the modular curve.
        e_2 : int
            The number of elliptic points of order 2.
        e_3 : int
            The number of elliptic points of order 3.
        e_inf : int
            The number of cusps.
        """
        N = self.N
        mu = int(N * sum([QQ(1, a) for a in sq_free_divisors(N)]))
        e_inf = sum([euler_phi(gcd(d, N // d)) for d in divisors(N)])
        e_3 = len([x for x in range(0, N) if (x**2 + x + 1) % N == 0])
        e_2 = len([x for x in range(0, N) if (x**2 + 1) % N == 0])
        return mu, e_2, e_3, e_inf

    def genus(self):
        """
        Computes the genus of X_0(N).

        Returns
        -------
        int
            The genus of X_0(N)
        """
        return RH_genus(self.mu, self.e_2, self.e_3, self.e_inf)

    def make_cusps(self):
        """
        Representatives for classes of cusps.

        Returns
        -------
        cusps : list
            A list of tuples (d,x) where d|N and x runs over reps for
            ℤ/gcd(N,N/d)ℤ. The order is chosen so that the Fricke involution
            acts as cusps.reverse().
        widths : list
            A list of cusp widths corresponding directly to the cusps.

        Note
        ----
        The list `cusps` of pairs (d,x) bijects with Γ_0(N) \\ ℙ^1(ℚ) under
        (d, x) -> [x/d]. The width of such a cusp is easy to check as the least
        element 0 < w < N such that w*d**2 == 0 (N)
        """
        N = self.N
        cusps = []
        widths = []
        for d in [int(d) for d in divisors(N)]:
            delta = gcd(d, N // d)
            w = N // gcd(N, d**2)
            for x in [x for x in range(0, delta) if gcd(x, delta) == 1]:
                cusps.append((d, x))
                widths.append(w)
        cusps, widths = fricke_sort(N, cusps, widths=widths)
        return cusps, widths


class X_w:
    """
    Class of the modular curve X_w for odd N

    Parameters
    ----------
    N : int
        An odd integer which is the level of the modular curve X_w
    r : int
        An integer coprime to `N`

    Attributes
    ----------
    N : int
        The level of X_w
    r : int
        An integer coprime to `N`
    mu : int
        The index of the modular curve.
    e_2 : int
        The number of elliptic points of order 2.
    e_3 : int
        The number of elliptic points of order 3.
    e_inf : int
        The number of cusps.
    """

    def __init__(self, N, r):
        if N % 2 == 0:
            raise ValueError("N should be odd")
        if gcd(N, r) != 1:
            raise ValueError("You need N,r to be coprime")
        self.N = N
        self.r = r
        self.mu, self.e_2, self.e_3, self.e_inf = self.init_ramification()

    def init_ramification(self):
        """
        Gets data about the map X_w -> X(1)

        Returns
        -------
        mu : int
            The index of the modular curve.
        e_2 : int
            The number of elliptic points of order 2.
        e_3 : int
            The number of elliptic points of order 3.
        e_inf : int
            The number of cusps.
        """
        N = self.N
        r = self.r
        if N == 1:
            mu = 1
            e_2 = 1
            e_3 = 1
            e_inf = 1
        else:
            mu = N**2 * prod(
                [QQ(int(p + kronecker(-r, p)), int(p)) for p in prime_divisors(N)]
            )
            e_2 = rho(N, r)
            e_3 = rho(N, 3 * r)
            e_inf = 1
            for p in prime_divisors(N):
                s = sum(
                    [
                        (1 + kronecker(-r, p)) * euler_phi(p**l)
                        for l in range(0, p_adic_val(N, p))
                    ]
                )
                e_inf = e_inf * (euler_phi(p ** p_adic_val(N, p)) + s)
            e_inf = int(e_inf)
        return mu, e_2, e_3, e_inf

    def genus(self):
        """
        Computes the genus of X_w.

        Returns
        -------
        int
            The genus of X_w
        """
        return RH_genus(self.mu, self.e_2, self.e_3, self.e_inf)


class Xg_plus:
    """
    Class of the modular curves X_g^+.

    Parameters
    ----------
    N : int
        An odd integer which is the level of the modular curves X_g^+.
    r : int
        An integer coprime to `N`

    Attributes
    ----------
    N : int
        The maximum level of X_g^+, should not be 1.
    r : int
        An integer coprime to `N`.
    M : int
        The largest odd divior of `N`.
    k : int
        An integer so that `N` == `M` * 2**`k`
    mu : int
        The index of the modular curve on X_w of level M.
    e_2 : int
        The number of elliptic points of order 2 on X_w of level M.
    e_3 : int
        The number of elliptic points of order 3 on X_w of level M.
    e_inf : int
        The number of cusps on X_w of level M.
    mup : fractions.Fraction
        The index of the modular curve X_w^+ of level M.
    e_2p : fractions.Fraction
        The number of elliptic points of order 2 on X_w^+ of level M.
    e_3p : fractions.Fraction
        The number of elliptic points of order 3 on X_w^+ of level M.
    e_infp : fractions.Fraction
        The number of cusps on X_w^+ of level M.

    Notes
    -----
    The rational numbers `mup`, `e_2p`, `e_3p`, `e_infp` are all elements of
    1/2*ℤ and are integers if `M` > 1.
    """

    def __init__(self, N, r):
        if N <= 1:
            raise ValueError("N should be at least 2")
        if gcd(N, r) != 1:
            raise ValueError("You need N,r to be coprime")
        self.N = N
        self.r = r
        self.k = p_adic_val(N, 2)
        self.M = N // (2**self.k)
        Xw_M = X_w(self.M, self.r)
        self.mu = Xw_M.mu
        self.e_2 = Xw_M.e_2
        self.e_3 = Xw_M.e_3
        self.e_inf = Xw_M.e_inf
        self.mup, self.e_2p, self.e_3p, self.e_infp = self.init_ramification()

    def init_ramification(self):
        """
        Computes the rational numbers μ+, e_2+, e_3+, e_∞+ (these integers when
        N != 1).

        Returns
        -------
        mup : fractions.Fraction
            The index of the modular curve.
        e_2p : fractions.Fraction
            The number of elliptic points of order 2.
        e_3p : fractions.Fraction
            The number of elliptic points of order 3.
        e_infp : fractions.Fraction
            The number of cusps.
        """
        M = self.M
        mup = QQ(self.mu, 2)
        e_2p = QQ(self.e_2, 2) + r_0(M)
        e_3p = QQ(self.e_3, 2)
        e_infp = QQ(self.e_inf, 2)
        return mup, e_2p, e_3p, e_infp

    def genus(self):
        """
        Computes the genus of X_g^+ for each g such that g^2 = -I.

        Returns
        -------
        list
            The genera of the modular curves X_g^+. The order of the components
            is given as follows:
            - `self.k` == 0 : [w]
            - `self.k` == 1 : [(I,w), (c,w)]
            - `self.k` == 2 and r == 1 : [(c,w)]
                            and r == 3 : [(ns, w), (s, w), (c,w)]
            - `self.k` >= 3 and r == 1 : [(c,w)]
            - `self.k` >= 3 and r == 3 : [(ns, w), (ωns, w), (c,w)]
            - `self.k` >= 3 and r == 5 : [(c,w)]
            - `self.k` >= 3 and r == 7 : [(s, w), (ωs, w), (c,w)]
        """
        N = self.N
        r = self.r
        k = self.k
        M = self.M
        mu = self.mu
        e_2 = self.e_2
        e_3 = self.e_3
        e_inf = self.e_inf
        mup = self.mup
        e_2p = self.e_2p
        e_3p = self.e_3p
        e_infp = self.e_infp

        if k == 0:
            ret = [RH_genus(mup, e_2p, e_3p, e_infp)]

        elif k == 1:
            X_Iw = Xg_plus(M, r)
            g_Iw = X_Iw.genus()[0]
            g_cw = RH_genus(
                3 * mup, 
                e_2p, 
                0, 
                2 * e_infp
            )
            ret = [g_Iw, g_cw]

        elif k == 2:
            if r % 4 == 1:
                gw = RH_genus(
                    3 * 2 ** (2 * k - 2) * mup, 
                    2 * e_2, 
                    0, 
                    4 * e_infp
                )
                ret = [gw]
            elif r % 4 == 3:
                gns = RH_genus(
                    2 * mup, 
                    r_0(2 * M), 
                    e_3, 
                    e_infp
                )
                gs = RH_genus(
                    6 * mup, 
                    r_0(2 * M), 
                    0, 
                    3 * e_infp
                )
                gw = RH_genus(
                    12 * mup, 
                    r_0(N), 
                    0, 
                    4 * e_infp
                )
                ret = [gns, gs, gw]

        elif k >= 3:
            if r % 8 == 1:
                gw = RH_genus(
                    3 * 2 ** (2 * k - 2) * mup, 
                    4 * e_2, 
                    0, 
                    (2**k) * e_infp
                )
                ret = [gw]
            elif r % 8 == 3:
                gns = RH_genus(
                    2 ** (2 * k - 3) * mup,
                    r_0(2 ** (k - 1) * M),
                    e_3,
                    2 ** (k - 2) * e_infp,
                )
                gw = RH_genus(
                    3 * 2 ** (2 * k - 2) * mup, 
                    r_0(N), 
                    0, 
                    (2**k) * e_infp
                )
                ret = [gns, gns, gw]
            elif r % 8 == 5:
                gw = RH_genus(
                    3 * 2 ** (2 * k - 2) * mup, 
                    0, 
                    0, 
                    (2**k) * e_infp
                )
                ret = [gw]
            elif r % 8 == 7:
                gs = RH_genus(
                    3 * 2 ** (2 * k - 3) * mup,
                    r_0(2 ** (k - 1) * M),
                    0,
                    3 * 2 ** (k - 2) * e_infp,
                )
                gw = RH_genus(
                    3 * 2 ** (2 * k - 2) * mup, 
                    r_0(N), 
                    0, 
                    (2**k) * e_infp
                )
                ret = [gs, gs, gw]
        return ret