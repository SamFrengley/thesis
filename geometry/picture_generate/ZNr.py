FIGURE_DIMENSIONS = [6.6, 2]

from math import gcd
from fractions import Fraction as QQ
from sympy import totient as euler_phi, divisors
import numpy as np
from .arithmetic import *
from .modularcurves import *
from .texassist import *

prod = np.prod

import matplotlib
from matplotlib.backends.backend_pgf import FigureCanvasPgf

matplotlib.backend_bases.register_backend("pdf", FigureCanvasPgf)
matplotlib.rcParams.update(
    {
        "pgf.texsystem": "pdflatex",
        "font.family": "serif",
        "text.usetex": True,
        "pgf.rcfonts": False,
        "figure.figsize": FIGURE_DIMENSIONS,
    }
)
import matplotlib.pyplot as plt

plt.close("all")


def find_small_r(N, q):
    """
    Finds the smallest integer in (ℤ/Nℤ)* which in the same class as
    q modulo ((ℤ/Nℤ)*)^2

    Parameters
    ----------
    N : int
        Modulus of the ring ℤ/`N`ℤ
    q : int
        An integer coprime to `N`

    Returns
    -------
    r : int
        The smallest integer in the same square class as `q` modulo `N`
    """
    found = False
    r = 0
    while not found:
        r += 1
        if rho(N, r * q) != 0:
            found = True
    return r


def format_table_same_col_widths(*ls, sep="|"):
    """
    Aligns lists in a table with columns of the same widths.

    Parameters
    ----------
    *ls :
        Lists of the same length consisting of objects (with __str__ methods).
    sep : str, default="|"
        Column separator for the output table

    Returns
    -------
    str
        A string which is a table formatted with the entries as the elements of
        `ls`[i], and the rows as the `ls`[i].
    """
    if len(set([len(l) for l in ls])) != 1:
        raise ValueError("All lists should be same length")
    ls = [[str(a) for a in l] for l in ls]
    ret = ["" for l in ls]
    max_l = max([max([len(a) for a in l]) for l in ls])
    max_l = 2 * (max_l // 2 + 1)

    for i in range(0, len(ls[0])):
        frmt = f"|{{:^{str(max_l + 1)}}}"
        for j in range(0, len(ret)):
            ret[j] += frmt.format(ls[j][i])
    for j in range(0, len(ret)):
        ret[j] += sep
    return ret


##################################################
######### HELPER FUNCTIONS FOR F_m ###############
##################################################


def comp_mults(N, d, q):
    """
    The multiplicities of the components of a cusp resolution on ~Z_(N,r)

    Parameters
    ----------
    N : int
        Level of the HMS
    d : int
        Integer dividing N, singularity of type (d,q)
    q : int
        Integer coprime to `d`, singularity of type (d,q)

    Returns
    -------
    x, y : list, list
        Lists of integers, the multiplicities of the components of the
        resolution in the divisors j*(∞) and (j')*(∞) respecitively
    """
    qq = q % d
    c = HJ_continued_fraction(d, qq)
    M = [[1] + [0 for i in range(0, len(c) + 1)]]
    for i in range(0, len(c)):
        M.append(
            [0 for j in range(1, i + 1)]
            + [-1, c[i], -1]
            + [0 for j in range(1, len(c) - i)]
        )
    M.append([0 for i in range(0, len(c) + 1)] + [1])
    M = np.array(M)
    S = [N] + [0 for i in range(0, len(c) + 1)]
    S = np.array(S)
    x = np.linalg.solve(M, S)
    x = [round(a) for a in x]
    S = [0 for i in range(0, len(c) + 1)] + [N]
    S = np.array(S)
    y = np.linalg.solve(M, S)
    y = [round(a) for a in y]
    return x, y


def err(N, m):
    """
    The correction term 𝖊 in the formula for the intersection of the canonical
    and a HZ divisor on ~Z_(N,r)

    Parameters
    ----------
    N : int
        Level of the HMS
    m : int
        Integer coprime to `m`, the level of the HZ divisor F_(m,λ)

    Returns
    -------
    int
        The intersection number K.~F_(m,λ)
    """
    ST = []
    kks = []
    for n in divisors(m):
        q = [q for q in range(1, N) if m % N == (q * n**2) % N][0]
        a, ap = comp_mults(N, N, q)
        ap2 = [0] + ap[1 : len(ap)]
        k = [
            k
            for k in range(1, len(a))
            if (a[k - 1] * n**2 > m * ap2[k - 1]) and (a[k] * n**2 <= m * ap2[k])
        ]
        k = k[0]
        kks.append((n, q, k))
        M = np.array([[a[k], a[k - 1]], [ap[k], ap[k - 1]]])
        S = np.array([m // n, n])
        x = np.linalg.solve(M, S)
        ST.append([round(a) for a in x])
    ret = sum(
        [
            euler_phi(gcd(st[0], st[1])) * (QQ(st[0] + st[1], gcd(st[0], st[1])) - 1)
            for st in ST
        ]
    )
    return (int(ret), ST, kks)


class F_m:
    """
    Class of the modular curve F_(m,λ) on ~Z_(N,r)

    Parameters
    ----------
    N : int
        Level of HMS
    m : int
        Level of HZ divisor (uniquely determines the power of the congruence)

    Arguments
    ---------
    N : int
        Level of the HMS
    m : int
        Level of the HZ divisor
    cusps : list
        List of tuples (`N`, `q`) where F_m meets resolution of a cusp
        singularity of type (`N`, `q`)
    comps : list
        List of tuples (k) or (k,k+1) where the corresponding cusp of F_m meets
        the resolution at the components E^(k) (resp. at the intersection of
        E^(k) and E^(k+1)
    int_mults : list
        List of tuples (same length as corresponding one in `comps`) which are
        the multiplicities of the relevant intersections
    X0_cusps : list
        Cusps on X_0(`m`), see modularcurves.X_0.make_cusps()
    X0_widths : list
        Cusp widths on X_0(`m`), see modularcurves.X_0.make_cusps()
    """

    def __init__(self, N, m):
        self.N = N
        self.m = m
        X = X_0(self.m)
        self.X0_cusps = X.cusps
        self.X0_widths = X.widths
        self.make_cusps()

    def make_cusps(self):
        """
        Computes data about the cusps on this HZ divisor.

        Returns
        -------
        None
        """
        self.cusps = []
        for n, _ in self.X0_cusps:
            q = self.m * pow(n, -2, self.N) % self.N
            self.cusps.append((self.N, q))

        self.comps = []
        self.int_mults = []
        for i in range(0, len(self.cusps)):
            a, ap = comp_mults(self.N, self.N, self.cusps[i][1])
            found = False
            k = -1
            while not found:
                k += 1
                M = np.array([[a[k], a[k + 1]], [ap[k], ap[k + 1]]])
                S = np.array([self.X0_widths[i], self.X0_widths[-1 - i]])
                st = np.linalg.solve(M, S)
                st_int = [round(x) for x in st]
                if st[0] > -0.01 and st[1] > -0.01:  # check are positive
                    st_diff = [st[jj] - st_int[jj] for jj in [0, 1]]
                    if (
                        -0.01 < st_diff[0] < 0.01 and -0.01 < st_diff[1] < 0.01
                    ):  # check are ints
                        found = True
                        if st_int[1] == 0:
                            kk = [k]
                            inters = [st_int[0]]
                        elif st_int[0] == 0:
                            kk = [k + 1]
                            inters = [st_int[1]]
                        else:
                            kk = [k, k + 1]
                            inters = st_int

            self.comps.append(kk)
            self.int_mults.append(inters)


##################################################
############ THE CLASS OF THE CUSPS ##############
##################################################


class ZtilCusp:
    """
    Cusp resolution E_(∞,d,q) on ~Z(N,r).

    Parameters
    ----------
    N : int
        Level of HMS.
    r : int
        Power of the congruence.
    d : int
        Integer dividing `N`, singularity of type (`d`,`q`).
    q : int
        Integer coprime to `d` and so that `q`*`r` is a square mod d.
    makesketch : bool, default=True
        If True then a sketch is auto-generated in the background.

    Arguments
    ---------
    N : int
        Level of HMS.
    r : int
        Power of the congruence.
    d : int
        Integer dividing `N`, singularity of type (`d`,`q`).
    q : int
        Integer coprime to `d` and so that `q`*`r` is a square mod d.
    ctd_frac : list
        The Hirzebruch--Jung continued fraction of `d`/`q`.
    mults : list
        List of tuples consisting of the multiplicities of the components of
        the resolution in the divisors j*(∞) and (j')*(∞) respecitively.
    fig : matplotlib.Figures.figure
        A figure representing the intersections near this resolution.
    ax : matplotlib.pyplot.axis
        An (the only) axis on `fig` representing the intersections near this
        resolution.
    title : str
        A LaTex title for this resolution.
    """

    def __init__(self, N, r, d, q, makesketch=True):
        if d <= 1 or N <= 1 or N % d != 0:
            raise ValueError("need d|N and positive")
        if gcd(r, N) != 1:
            raise ValueError("N and r need to be coprime")
        if gcd(q, d) != 1:
            raise ValueError("d and q need to be coprime")
        if rho(d, q * r) == 0:
            raise ValueError("q and r need to be same up to sq mod d")
        self.N = N
        self.r = r
        self.d = d
        q = q % d
        self.q = q
        self.ctd_frac = HJ_continued_fraction(d, q)
        self.mults = comp_mults(N, d, q)

        if makesketch:
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot()
            self.ax.axis("off")
            srf = f"$\\widetilde{{Z}}_{{{self.N},{self.r}}}$"
            self.title = f"Intersection behaviour in a neighbourhood of (a chunk of) $E_{{\\infty,{self.d},{self.q}}}$ on {srf}"
            self.fig.suptitle(self.title)
            if len(self.ctd_frac) % 2 == 1:
                self.oddstyle = True
            else:
                self.oddstyle = False

    def __str__(self, sketch=False):
        N = self.N
        r = self.r
        d = self.d
        q = self.q
        Z = Ztil(N, r)
        ret = f"\nThe divisor E_(∞,{d},{q}) on ~Z_({N},{r})\n"
        chunks = (rho(d, q * r) * int(euler_phi(N // d))) // 2
        ret += f"    <> It has {chunks} chunk(s)\n"
        ret += f"    <> The curves C_(∞,i) have genus {X_1(N).genus()}\n"
        ret += "\n" + self.intersection_table()
        if sketch:
            self.sketch()
        return ret

    def is_fixed(self):
        """
        Checks this cusp is fixed by the involution 𝜏

        Returns
        -------
        bool
            true if and only if E_(∞,d,q) is fixed by 𝜏
        """
        N = self.N
        r = self.r
        d = self.d
        q = self.q
        ret = False
        if q == d - 1:
            ret = True
        elif d == N and q == 1:
            ret = True
        elif N % 2 == 0:
            if d == N // 2 and q == 1:
                ret = True
        return ret

    def meets_Fm(self, m):
        """
        Asks if this cusp meets the resolution of the HZ divisor F_(m)

        Parameters
        ----------
        m : int
            Level of the HZ divisor F_(m)

        Returns
        -------
        bool
            True when E_(∞,d,q) meets F_(m)
        """
        if self.d == self.N:
            for n in divisors(m):
                if self.q * int(n) ** 2 % self.d == m:
                    return True
        return False

    def where_meets_Fm(self, m):
        """
        The components of E_(∞,d,q) meeting F_m.

        Parameters
        ----------
        m : int
            Level of the HZ divisor

        Returns
        -------
        list
            A list of list of length 1 or 2. If the length is 2 then it  meets
            the resolution at an intersection of E^(k) and E^(k+1).

        list
            Multiplicities of the intersections.
        """
        F = F_m(self.N, m)
        cusps = F.cusps
        comps = F.comps
        mults = F.int_mults
        inds = [i for i, x in enumerate([y for _, y in cusps]) if x == self.q]
        return [comps[i] for i in inds], [mults[i] for i in inds]

    def intersection_table(self):
        """
        Self intersections and mults of divisors

        Returns
        -------
        str
            A formatted string which is a table of components and their
            respective self intersections.
        """
        N = self.N
        r = self.r
        d = self.d
        q = self.q
        Z = Ztil(N, r)

        names = (
            ["C_(∞,1)"]
            + [f"E^({i})" for i in range(1, len(self.ctd_frac) + 1)]
            + ["C_(∞,2)"]
        )
        si = (
            [str(Z.Cinf_Cinf())]
            + [str(-a) for a in self.ctd_frac]
            + [str(Z.Cinf_Cinf())]
        )
        muls = [
            str((self.mults[0][i], self.mults[1][i]))
            for i in range(0, len(self.mults[0]))
        ]

        strs = format_table_same_col_widths(names, si, muls)
        strs[0] = "|{:^20}|{}".format("Curve", strs[0])
        strs[1] = "|{:^20}|{}".format("Self intersections", strs[1])
        strs[2] = "|{:^20}|{}".format("Multiplicities", strs[2])

        ret = len(strs[0]) * ("-") + "\n"
        ret += strs[0] + "\n"
        ret += len(strs[0]) * ("-") + "\n"
        for s in strs[1:]:
            ret += s + "\n"
        ret += len(strs[0]) * ("-") + "\n"
        return ret

    def sketch(
        self,
        E_labels=False,
        SI=True,
        F_g=True,
        F_m=True,
        m_list=get_rat_X0_plus(),
        use_oddstyle=True,
        eps=0.1,
    ):
        """
        A sketch of the intersections of the components

        Parameters
        ----------
        E_labels : bool, default=False
            If True then the figure has labels on all the exceptional divisors
        SI : bool, default=True
            If True then the self intersections of the exceptional divisors are
            printed.
        F_g : bool, default=True
            If True then the fixed curves F_g appear on the figure
        m_list : list, default=get_rat_X0_plus()
            A list of integers m > 0 for which we display the HZ divisors F_m
            on the figure. By default we only show those with X_0^+(m) (Fricke
            quotient) is rational.

        Returns
        -------
        fig : matplotlib.figure.Figure
            A matplotlib figure type, being the picture of the intersection
            around this cusp resolution.

        Other parameters
        ----------------
        use_oddstyle : bool, default=True
            If True, then all output figures will be symmetric across the
            y-axis
        eps : float, default=0.1
            A spacing parameter for how far the labels should be from the lines
            in the figures
        """
        num = len(self.ctd_frac)
        up = (num // 2) + (num % 2)
        down = num // 2
        if num % 2 == 1 and use_oddstyle:
            self.oddstyle = True
        else:
            self.oddstyle = False
        self.printed_on_components = [0 for i in range(0, num)]
        self.sketch_cyc_quot_res(E_labels=E_labels, eps=eps, SI=SI)
        self.sketch_Cinf(eps=eps)
        if F_g:
            self.sketch_Fg(eps=eps)
            if self.is_fixed() and num % 2 == 1:
                self.printed_on_components[num // 2] = 2
        if F_m:
            for m in [m for m in m_list if gcd(m, self.N) == 1]:
                if self.meets_Fm(m):
                    try:
                        self.sketch_Fm(m, eps=eps)
                    except NotImplementedError:
                        print(
                            f"We ran out of space on a component when we got to m={m}"
                        )
        extra = 1.2 * 0.25 * int(self.oddstyle)
        self.ax.set_xlim(-0.1, 0.6 * num + 0.4 + extra + 0.1)
        self.ax.set_ylim(-1.75, 2.75)
        return self.fig

    def sketch_cyc_quot_res(self, E_labels=False, SI=True, eps=0.1):
        """
        Sketches the resolution of the component of E_(∞,d,q)

        Parameters
        ----------
        E_labels : bool, default=False
            If True then the figure has labels on all the exceptional divisors.
        SI : bool, default=True
            If True then the self intersections of the exceptional divisors are
            printed.

        Returns
        -------
        None

        Other parameters
        ----------------
        eps : float, default=0.1
            A spacing parameter for how far the labels should be from the lines
            in the figures.
        """
        if len(self.ctd_frac) % 2 == 1 and self.oddstyle:
            return self.sketch_cyc_quot_res_odd(E_labels=E_labels, SI=SI, eps=eps)
        else:
            return self.sketch_cyc_quot_res_even(E_labels=E_labels, SI=SI, eps=eps)

    def sketch_cyc_quot_res_even(self, E_labels=False, SI=True, eps=0.1):
        """
        Subroutine for `sketch_cyc_quot_res` in the case when the number of
        components of E_(∞,d,q) has even length

        Parameters
        ----------
        E_labels : bool, default=False
            If True then the figure has labels on all the exceptional divisors
        SI : bool, default=True
            If True then the self intersections of the exceptional divisors are
            printed

        Returns
        -------
        None

        Other parameters
        ----------------
        eps : float, default=0.1
            A spacing parameter for how far the labels should be from the lines
            in the figures
        """
        num = len(self.ctd_frac)
        up = (num // 2) + (num % 2)
        down = num // 2
        for i in range(0, up):
            x = (1.2 * i, 1.2 * i + 1)
            y = (0, 1)
            self.ax.plot(x, y, color="black")
            if E_labels:
                s = f"$E^{{({2*i + 1})}}$"
                self.ax.text(x[1] + eps, y[1], s)
            if SI:
                s = f"${-self.ctd_frac[2*i]}$"
                self.ax.text(1.2 * i + 0.5, 0.5 - 0.05 - 2 * eps, s)
        for i in range(0, down):
            x = (1.2 * i + 0.6, 1.2 * i + 1.6)
            y = (1, 0)
            self.ax.plot(x, y, color="black")
            if E_labels:
                s = "$E^{{({2*i + 2})}}$"
                self.ax.text(x[1] + eps, y[1], s)
            if SI:
                s = f"${-self.ctd_frac[2*i + 1]}$"
                self.ax.text(1.2 * i + 1.1, 0.5 + 2 * eps, s)

    def sketch_cyc_quot_res_odd(self, E_labels=False, SI=True, eps=0.1):
        """
        Subroutine for `sketch_cyc_quot_res` in the case when the number of
        components of E_(∞,d,q) has odd length

        Parameters
        ----------
        E_labels : bool, default=False
            If True then the figure has labels on all the exceptional divisors
        SI : bool, default=True
            If True then the self intersections of the exceptional divisors are
            printed

        Returns
        -------
        None

        Other parameters
        ----------------
        eps : float, default=0.1
            A spacing parameter for how far the labels should be from the lines
            in the figures
        """
        self.sketch_cqr_odd_part(0, E_labels=E_labels, SI=SI, eps=eps)
        self.sketch_cqr_odd_part(1, E_labels=E_labels, SI=SI, eps=eps)
        halfway = len(self.ctd_frac) // 2
        at_top = (len(self.ctd_frac) - 1) % 4 != 0
        x = (0.6 * halfway, 0.6 * halfway + 1 + 1.2 * 0.25)
        if at_top:
            y = (0.8, 0.8)
        else:
            y = (0.2, 0.2)
        self.ax.plot(x, y, color="black")
        s = f"${-self.ctd_frac[halfway]}$"
        self.ax.text(0.5 * (x[1] + x[0]), y[0] + 2 * eps, s, ha="center")

    def sketch_cqr_odd_part(self, part, E_labels=False, SI=True, eps=0.1):
        """
        Subroutine for `sketch_cyc_quot_res_odd`. Sketches the first or second
        part of E_(∞,d,q) when the continued fraction has odd length.

        Parameters
        ----------
        part : int
            Either 0 or 1 depending on whether we are sketching the first or
            second half of the picture
        E_labels : bool, default=False
            If True then the figure has labels on all the exceptional divisors
        SI : bool, default=True
            If True then the self intersections of the exceptional divisors are
            printed

        Returns
        -------
        None

        Other parameters
        ----------------
        eps : float, default=0.1
            A spacing parameter for how far the labels should be from the lines
            in the figures
        """
        num = len(self.ctd_frac)
        num = num // 2
        part_c = self.ctd_frac[part * num + part : (part + 1) * num + part]
        up = (num // 2) + (num % 2)
        down = num // 2
        start = 1.2 * part * (down + 0.5 + 0.25 * ((len(self.ctd_frac) - 1) % 4) + 0.25)
        start_up = part * (num % 2) == 0  # chooses if second one goes up or down
        shift_dir = (-1) ** ((len(self.ctd_frac) - 1) % 4)

        for i in range(0, up):
            x = (start + 1.2 * i, start + 1.2 * i + 1)
            if start_up:
                y = (0, 1)
            else:
                y = (1, 0)
            self.ax.plot(x, y, color="black")
            if E_labels:
                s = f"$E^{{({2*i + 1})}}$"
                self.ax.text(x[1] + eps, y[1], s)
            if SI:
                s = f"${-part_c[2*i]}$"
                self.ax.text(
                    start + 1.2 * i + 0.5 - shift_dir * 0.05,
                    0.5 + shift_dir * (0.05 + 2 * eps),
                    s,
                )
        for i in range(0, down):
            x = (start + 1.2 * i + 0.6, start + 1.2 * i + 1.6)
            if not start_up:
                y = (0, 1)
            else:
                y = (1, 0)
            self.ax.plot(x, y, color="black")
            if E_labels:
                s = f"$E^{{({2*i + 2})}}$"
                self.ax.text(x[1] + eps, y[1], s)
            if SI:
                s = f"${-part_c[2*i + 1]}$"
                self.ax.text(
                    start + 1.2 * i + 1.1 - shift_dir * 0.05,
                    0.5 + shift_dir * (2 * eps),
                    s,
                )

    def sketch_Cinf(self, eps=0.1):
        """
        Sketches C_(∞).

        Returns
        -------
        None

        Other parameters
        ----------------
        eps : float, default=0.1
            A spacing parameter for how far the labels should be from the lines
            in the figures
        """
        num = len(self.ctd_frac)
        # C_(inf,1)
        x = (0.2, 0.2)
        y = (2, -1)
        self.ax.plot(x, y, color="black")
        self.ax.text(
            x[0], y[0] + eps, "$\\widetilde{C}_{\\infty, 1}$", ha="right", va="bottom"
        )
        # C_(inf,2)
        if self.oddstyle:
            x = (0.6 * num + 0.2 + 1.2 * 0.25, 0.6 * num + 0.2 + 1.2 * 0.25)
        else:
            x = (0.6 * num + 0.2, 0.6 * num + 0.2)
        y = (2, -1)
        self.ax.plot(x, y, color="black")
        self.ax.text(
            x[0], y[0] + eps, "$\\widetilde{C}_{\\infty, 2}$", ha="left", va="bottom"
        )

    def sketch_Fg(self, FPs=True, non_fixed=True, eps=0.1):
        """
        Sketch components of F_g meeting E_(∞,d,q).

        Parameters
        ----------
        FPs : bool, default=True
            True if and only if we display dots for isolated fixed points
        non_fixed : bool, default=True
            If True then print the non-fixed curves which get contracted to
            remove isolated fixed points

        Returns
        -------
        None

        Other parameters
        ----------------
        eps : float, default=0.1
            A spacing parameter for how far the labels should be from the lines
            in the figures

        Notes
        -----
        This is just a huge lookup table, and is implimented as such. You would
        do better to refer to Proposition ?? from the text.
        """
        N = self.N
        num = len(self.ctd_frac)
        k = p_adic_val(N, 2)
        M = N // 2**k
        if self.is_fixed():
            if self.oddstyle:
                x = 0.5 * (0.6 * num + 0.2 + 1.2 * 0.25 + 0.2) - 0.15
            else:
                x = 0.5 * (0.6 * num + 0.2 + 0.2)
            ##########
            if self.q != 1 or self.d == 2:
                if self.oddstyle:
                    self.ax.plot((x, x), (1.5, -0.5), color="black", linewidth=4)
                    self.ax.plot(
                        (x + 0.3, x + 0.3), (1.5, -0.5), color="black", linewidth=4
                    )
                else:
                    self.ax.plot((x, x), (1.5, -0.5), color="black", linewidth=4)
            ##########
            if k == 0:
                if self.q != 1:
                    self.ax.text(x, 1.5 + eps, format_Fg("w"), ha="center", va="bottom")
                if self.q == 1:
                    self.ax.plot((x, x), (1.5, -0.5), color="black", linewidth=4)
                    self.ax.text(
                        x, 1.5 + eps, "$\\widetilde{F}_{1}$", ha="center", va="bottom"
                    )
                    if FPs:
                        self.ax.plot(x + 0.3, 0.2, "o", color="black")
            ##########
            elif k == 1:
                if self.q != 1 or self.d == 2:
                    if self.oddstyle:
                        self.ax.text(
                            x, 1.5 + eps, format_Fg("I", "w"), ha="center", va="bottom"
                        )
                        self.ax.text(
                            x + 0.3,
                            1.5 + eps,
                            format_Fg("c", "w"),
                            ha="center",
                            va="bottom",
                        )
                    else:
                        self.ax.plot((x, x), (1.5, -0.5), color="black", linewidth=4)
                        self.ax.text(
                            x, 1.5 + eps, format_Fg("c", "w"), ha="center", va="bottom"
                        )

                if self.q == 1:
                    if self.d == self.N:
                        self.ax.plot((x, x), (1.5, -0.5), color="black", linewidth=4)
                        self.ax.plot(
                            (x + 0.3, x + 0.3), (1.5, -0.5), color="black", linewidth=4
                        )
                        self.ax.text(
                            x,
                            1.5 + eps,
                            "$\\widetilde{F}_{1, \\lambda}$",
                            ha="center",
                            va="bottom",
                        )
                        self.ax.text(
                            x + 0.3,
                            1.5 + eps,
                            format_Fg("b", "I", lam="\\lambda"),
                            ha="center",
                            va="bottom",
                        )
                    elif self.d == self.N // 2:
                        self.ax.plot((x, x), (1.5, -0.5), color="black", linewidth=4)
                        self.ax.text(
                            x,
                            1.5 + eps,
                            format_Fg("b", "I", lam="\\lambda"),
                            ha="center",
                            va="bottom",
                        )
                        if FPs:
                            self.ax.plot(x + 0.3, 0.2, "o", color="black")
                        if non_fixed:
                            self.ax.plot((x + 0.3, x + 0.3), (1.5, -0.5), color="black")
                            self.ax.text(
                                x + 0.3,
                                1.5 + eps,
                                format_Fg("I^{\\sharp}", "I", lam="\\lambda"),
                                ha="center",
                                va="bottom",
                            )
            ##########
            elif k == 2:
                if self.q != 1 or self.d == 2:
                    if self.d % 2 == 1:
                        self.ax.text(
                            x, 1.5 + eps, format_Fg("c", "w"), ha="center", va="bottom"
                        )
                    elif self.d % 4 == 2:
                        self.ax.text(
                            x, 1.5 + eps, format_Fg("c", "w"), ha="center", va="bottom"
                        )
                        self.ax.text(
                            x + 0.3,
                            1.5 + eps,
                            format_Fg("c", "w"),
                            ha="center",
                            va="bottom",
                        )
                    elif self.d % 4 == 0:
                        self.ax.text(
                            x, 1.5 + eps, format_Fg("s", "w"), ha="center", va="bottom"
                        )
                        both_labs = "50\\% {} \n50\\% {}".format(
                            format_Fg("s", "w"), format_Fg("ns", "w")
                        )
                        self.ax.text(
                            x + 0.3, 1.5 + eps, both_labs, ha="center", va="bottom"
                        )
                else:
                    if self.d == self.N:
                        self.ax.plot((x, x), (1.5, -0.5), color="black", linewidth=4)
                        self.ax.plot(
                            (x + 0.3, x + 0.3), (1.5, -0.5), color="black", linewidth=4
                        )
                        self.ax.text(
                            x,
                            1.5 + eps,
                            "$\\widetilde{F}_{1, \\lambda}$",
                            ha="center",
                            va="bottom",
                        )
                        self.ax.text(
                            x + 0.3,
                            1.5 + eps,
                            format_Fg("b", "I", lam="\\lambda"),
                            ha="center",
                            va="bottom",
                        )
                    elif self.d == self.N // 2:
                        if self.r % 4 == 1:
                            self.ax.plot(
                                (x, x), (1.5, -0.5), color="black", linewidth=4
                            )
                            self.ax.plot(
                                (x + 0.3, x + 0.3),
                                (1.5, -0.5),
                                color="black",
                                linewidth=4,
                            )
                            self.ax.text(
                                x,
                                1.5 + eps,
                                format_Fg("b", "I", lam="\\lambda"),
                                ha="center",
                                va="bottom",
                            )
                            self.ax.text(
                                x + 0.3,
                                1.5 + eps,
                                format_Fg("-\\mathrm{b}", "I", lam="\\lambda"),
                                ha="center",
                                va="bottom",
                            )
                        elif self.r % 4 == 3:
                            if FPs:
                                self.ax.plot(x, 0.2, "o", color="black")
                                self.ax.plot(x + 0.3, 0.2, "o", color="black")
                            if non_fixed:
                                self.ax.plot(
                                    (x + 0.3, x + 0.3), (1.5, -0.5), color="black"
                                )
                                self.ax.text(
                                    x + 0.3,
                                    1.5 + eps,
                                    format_Fg(
                                        "\\mathrm{b}^{\\sharp}", "I", lam="\\lambda"
                                    ),
                                    ha="center",
                                    va="bottom",
                                )
                                self.ax.plot((x, x), (1.5, -0.5), color="black")
                                self.ax.text(
                                    x,
                                    1.5 + eps,
                                    format_Fg("I^{\\sharp}", "I", lam="\\lambda"),
                                    ha="center",
                                    va="bottom",
                                )
            ##########
            else:
                if self.q != 1 or self.d == 2:
                    if self.d % 2 == 1:
                        self.ax.text(
                            x, 1.5 + eps, format_Fg("c", "w"), ha="center", va="bottom"
                        )
                    elif self.d % 4 == 2:
                        self.ax.text(
                            x, 1.5 + eps, format_Fg("c", "w"), ha="center", va="bottom"
                        )
                        self.ax.text(
                            x + 0.3,
                            1.5 + eps,
                            format_Fg("c", "w"),
                            ha="center",
                            va="bottom",
                        )
                    elif self.d % 4 == 0:
                        if self.r % 8 == 3:
                            self.ax.text(
                                x,
                                1.5 + eps,
                                format_Fg("ns", "w"),
                                ha="center",
                                va="bottom",
                            )
                            self.ax.text(
                                x + 0.3,
                                1.5 + eps,
                                format_Fg("\\omega\\mathrm{ns}", "w"),
                                ha="center",
                                va="bottom",
                            )
                        elif self.r % 8 == 7:
                            self.ax.text(
                                x,
                                1.5 + eps,
                                format_Fg("s", "w"),
                                ha="center",
                                va="bottom",
                            )
                            self.ax.text(
                                x + 0.3,
                                1.5 + eps,
                                format_Fg("\\omega\\mathrm{s}", "w"),
                                ha="center",
                                va="bottom",
                            )
                else:
                    if self.d == self.N:
                        self.ax.plot((x, x), (1.5, -0.5), color="black", linewidth=4)
                        self.ax.plot(
                            (x + 0.3, x + 0.3), (1.5, -0.5), color="black", linewidth=4
                        )
                        self.ax.text(
                            x,
                            1.5 + eps,
                            "$\\widetilde{F}_{1, \\lambda}$",
                            ha="center",
                            va="bottom",
                        )
                        self.ax.text(
                            x + 0.3,
                            1.5 + eps,
                            format_Fg("b", "I", lam="\\lambda"),
                            ha="center",
                            va="bottom",
                        )
                    elif self.d == self.N // 2:
                        if self.r % 8 == 1:
                            self.ax.plot(
                                (x, x), (1.5, -0.5), color="black", linewidth=4
                            )
                            self.ax.plot(
                                (x + 0.3, x + 0.3),
                                (1.5, -0.5),
                                color="black",
                                linewidth=4,
                            )
                            if k == 3:
                                self.ax.text(
                                    x,
                                    1.5 + eps,
                                    format_Fg("b", "I", lam="\\lambda"),
                                    ha="center",
                                    va="bottom",
                                )
                                self.ax.text(
                                    x + 0.3,
                                    1.5 + eps,
                                    format_Fg("-\\mathrm{b}", "I", lam="\\lambda"),
                                    ha="center",
                                    va="bottom",
                                )
                            elif k >= 4:
                                self.ax.text(
                                    x,
                                    1.5 + eps,
                                    "50\\%\n{}".format(
                                        format_Fg("b", "I", lam="\\lambda")
                                    ),
                                    ha="center",
                                    va="bottom",
                                )
                                self.ax.text(
                                    x + 0.3,
                                    1.5 + eps,
                                    "50\\%\n{}".format(
                                        format_Fg("-\\mathrm{b}", "I", lam="\\lambda")
                                    ),
                                    ha="center",
                                    va="bottom",
                                )
                                self.ax.text(
                                    x + 0.3,
                                    -0.8,
                                    "50\\%\n{}\n(unfixed)".format(
                                        format_Fg(
                                            "\\mathrm{b}^{\\sharp}", "I", lam="\\lambda"
                                        )
                                    ),
                                    ha="center",
                                    va="top",
                                )
                                self.ax.text(
                                    x,
                                    -0.8,
                                    "50\\%\n{}\n(unfixed)".format(
                                        format_Fg("I^{\\sharp}", "I", lam="\\lambda")
                                    ),
                                    ha="center",
                                    va="top",
                                )
                        if self.r % 8 == 5 and k == 3:
                            if FPs:
                                self.ax.plot(x, 0.2, "o", color="black")
                                self.ax.plot(x + 0.3, 0.2, "o", color="black")
                            if non_fixed:
                                self.ax.plot(
                                    (x + 0.3, x + 0.3), (1.5, -0.5), color="black"
                                )
                                self.ax.text(
                                    x + 0.3,
                                    1.5 + eps,
                                    format_Fg(
                                        "\\mathrm{b}^{\\sharp}", "I", lam="\\lambda"
                                    ),
                                    ha="center",
                                    va="bottom",
                                )
                                self.ax.plot((x, x), (1.5, -0.5), color="black")
                                self.ax.text(
                                    x,
                                    1.5 + eps,
                                    format_Fg("I^{\\sharp}", "I", lam="\\lambda"),
                                    ha="center",
                                    va="bottom",
                                )

    def sketch_Fm(self, m, eps=0.1):
        """
        Sketch components of F_m meeting E_(∞,d,q).

        Parameters
        ----------
        m : int
            Level of the HZ divisor F_`m`

        Returns
        -------
        None

        Other parameters
        ----------------
        eps : float, default=0.1
            A spacing parameter for how far the labels should be from the lines
            in the figures
        """
        if m <= 1:
            raise ValueError("m should be > 1")
        if self.meets_Fm(m):
            num = len(self.ctd_frac)
            # specific case where lands on fixed point
            if m == 4 and self.q == 1:
                x = 0.5 * (0.6 * num + 0.2 + 1.2 * 0.25 + 0.2) + 0.15
                self.ax.plot((x, x), (2, -1), color="black", linewidth=2)
                self.ax.text(
                    x, 2 + eps, "$\\widetilde{F}_{4}$", ha="center", va="bottom"
                )
            else:
                comps, mults = self.where_meets_Fm(m)
                for i in range(0, len(comps)):
                    self.sketch_Fm_comp(m, comps[i], mults[i])

    def sketch_Fm_comp(self, m, comp, mult, eps=0.1, Fm_col="black"):
        """
        Subroutine for sketch_Fm

        Parameters
        ----------
        m : int
            Level of the HZ divisor F_`m`.
        comp : tuple
            An int or a pairs of ints indicating the component(s) of the
            resolution which is met by F_`m`.
        mult : tuple
            An int or pair of ints indicating the multiplicity with which the
            components are met.
        Fm_col="black" : {"black"}, optional
            A matplotlib viable string for the colour of this particular curve.

        Returns
        -------
        None

        Other parameters
        ----------------
        eps : float, default=0.1
            A spacing parameter for how far the labels should be from the lines
            in the figures.
        """
        num = len(self.ctd_frac)
        if self.oddstyle and comp[0] > num // 2:
            adj = 1.2 * 0.25
        else:
            adj = 0
        label = f"$\\widetilde{{F}}_{{{m}}}$"
        if len(comp) == 2:
            if comp[0] == 0:
                x = (-0.1, 0.2 + 0.7 * 0.3 * (1.3) ** (-1))
                y = (1.5, -0.5)
                self.ax.plot(x, y, color=Fm_col, linewidth=2)
                self.ax.text(x[0] - eps, y[0] + eps, label, ha="center", va="bottom")
            elif comp[0] == num:
                tranl = 0.6 * num + adj
                x = (tranl, tranl + 1.3 * 0.3 * (0.7) ** (-1))
                y = (-0.5, 1.5)
                self.ax.plot(x, y, color=Fm_col, linewidth=2)
                self.ax.text(x[1], y[1] + eps, label, ha="center", va="bottom")
            else:
                x = (0.6 * comp[0] + 0.2 + adj, 0.6 * comp[0] + 0.2 + adj)
                y = (-0.5, 1.5)
                self.ax.plot(x, y, color=Fm_col, linewidth=2)
                self.ax.text(x[1], y[1] + eps, label, ha="center", va="bottom")
        else:
            sgn = (comp[0] - 1) >= (num // 2)
            sgn = (-1) ** (1 - int(sgn))
            if self.printed_on_components[comp[0] - 1] == 0:
                if self.oddstyle and comp[0] - 1 == num // 2:
                    mais = 0.4 - 0.5 * adj
                elif comp[0] == 1:
                    mais = 0.5 + sgn * 0.15
                elif comp[0] == num:
                    mais = 0.5 + sgn * 0.15
                else:
                    mais = 0.5 + sgn * 0.05
            elif self.printed_on_components[comp[0] - 1] == 1:
                if self.oddstyle and comp[0] - 1 == num // 2:
                    mais = 0.6 - 0.5 * adj
                else:
                    mais = 0.5 - sgn * 0.05
            elif self.printed_on_components[comp[0] - 1] == 2:
                if self.oddstyle and comp[0] - 1 == num // 2:
                    mais = 0.3 - 0.5 * adj
                else:
                    mais = 0.5
            elif self.printed_on_components[comp[0] - 1] == 3 and comp[0] == 1:
                mais = 0.7 - 0.5 * adj
            elif num == 1 and self.printed_on_components[comp[0] - 1] == 4:
                mais = 0.25 - 0.5 * adj
            elif num == 1 and self.printed_on_components[comp[0] - 1] == 5:
                mais = 0.75 - 0.5 * adj
            else:
                raise NotImplementedError("Sorry, too much stuff on this component")
            x = 0.6 * (comp[0] - 1) + adj + mais
            x = (x, x)
            y = (2, -1)
            self.ax.plot((x, x), (2, -1), color=Fm_col, linewidth=2)
            self.printed_on_components[comp[0] - 1] += 1
            self.ax.text(x[0], y[0] + eps, label, ha="center", va="bottom")


class Ztil:
    """
    This is the class of the minimal desingularisation of the HMS ~Z_{N,r}.

    Parameters
    ----------
    N : int
        The level `N` of the HMS
    r : int
        The power `r` of the congruence

    Attributes
    ----------
    N : int
        HMS of level N**2
    r : int
        Power of congruence
    M : int
        Odd part of N
    k : int
        So that N = M * 2**k
    """

    def __init__(self, N, r):
        self.N = N
        self.r = r
        self.k = p_adic_val(N, 2)
        self.M = N // (2**self.k)

    def __str__(self):
        ret = f"The surface ~Z_({self.N},{self.r})\n"
        ret += f"κ     = {self.kodaira_dim()}\n"
        ret += f"p_g   = {self.p_g()}\n"
        ret += f"K^2   = {self.K_K()}\n"
        ret += f"χ_top = {self.chi_top()}"
        return ret

    def p_g(self):
        """
        Compute the geometric genus of ~Z_{N,r}.

        Returns
        -------
        int
           The geometric genus of ~Z_{N,r}.
        """
        N = self.N
        r = self.r
        ret = index_m(N) * (N - 12)
        ret += (18 * N) * euler_phi(N)
        ret += (18 * N) * r_0(N)
        ret += (16 * N) * 2 * r_1(N)
        ret += -(16 * N) * s_31(N, r)
        ret += (48 * N) * r_inf(N)
        ret += (12 * N) * bbL_infr(N, r)
        ret += -(12 * N) * R_infr(N, r)
        ret += -(144 * N)
        return ret // (144 * N)

    def K_K(self):
        """
        Compute the self intersection of the canonical on ~Z_{N,r}.

        Returns
        -------
        int
            Self intersection K.K of the (a) canonical on ~Z_{N,r}.
        """
        N = self.N
        r = self.r
        ret = 3 * (index_m(N) * (N - 12))
        ret += 3 * (18 * N) * euler_phi(N)
        ret += -(18 * N) * s_31(N, r)
        ret += 3 * (18 * N) * 3 * r_inf(N)
        ret += -3 * (18 * N) * R_infr(N, r)
        return ret // (3 * 18 * N)

    def chi_top(self):
        """
        Compute the topological Euler characteristic of ~Z_{N,r}

        Returns
        -------
        int
            The topological Euler characteristic of ~Z_{N,r}
        """
        N = self.N
        r = self.r
        ret = index_m(N) * (N - 12)
        ret += (18 * N) * euler_phi(N)
        ret += (18 * N) * 3 * r_0(N)
        ret += (12 * N) * 8 * r_1(N)
        ret += -(36 * N) * s_31(N, r)
        ret += (36 * N) * r_inf(N)
        ret += (36 * N) * bbL_infr(N, r)
        return ret // (36 * N)

    def kodaira_dim(self):
        """
        Computes the Kodaira dimension of ~Z_{N,r}.

        Returns
        -------
        int
            The Kodaira dimension of ~Z_{N,r}.
        """
        N = self.N
        if N >= 13:  # not needed, but could get slow for large N otherwise
            kappa = 2
        else:
            kappa = min([2, self.p_g() - 1])
        return kappa

    def Cinf_Cinf(self):
        """
        Computes the self intersection of the curve ~C_{∞,i}.

        Returns
        -------
        int
            The self intersection of the curve ~C_{∞,i} (for either i=1,2).
        """
        N = self.N
        r = self.r
        ret = [
            euler_phi(gcd(v, N)) * rat_part(QQ(v**2 * r, N * gcd(N, v)))
            for v in range(1, N)
        ]
        return -sum(ret) // 2

    def K_Cinf(self):
        """
        Computes the intersection number K.~C_{∞,i} on ~Z_{N,r}.

        Returns
        -------
        int
            The intersection number K.~C_{∞,i} on ~Z_{N,r}.
        """
        N = self.N
        return 2 * X_1(N).genus() - 2 - self.Cinf_Cinf()

    def K_Fg(self):
        """
        Computes K.F_g for each F_g in the strict transform of F°.

        Returns
        -------
        list
            The intersection K.F_g for each F_g in the strict transform of F°.
            The order of the components
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

        X = Xg_plus(N, r)
        mup = X.mup
        e_2p = X.e_2p
        e_3p = X.e_3p
        e_infp = X.e_infp

        if k == 0:
            ret = [(mup - 6 * e_infp - e_3p) // 3]
        elif k == 1:
            ret = [(mup - 3 * e_infp - e_3p) // 3, mup - 3 * e_infp]
        elif k == 2:
            K_Fw = 2 ** (2 * k - 2) * mup - 3 * 2 ** (k - 1) * e_infp
            K_Fns = (2 ** (2 * k - 3) * mup - 3 * 2 ** (k - 2) * e_infp - 2 * e_3p) // 3
            K_Fs = 2 ** (2 * k - 3) * mup - 3 * 2 ** (k - 2) * e_infp
            if r % 4 == 1:
                ret = [K_Fw]
            elif r % 4 == 3:
                ret = [K_Fns, K_Fs, K_Fw]
        elif k >= 3:
            K_Fw = 2 ** (2 * k - 2) * mup - 3 * 2 ** (k - 1) * e_infp
            K_Fns = (2 ** (2 * k - 3) * mup - 3 * 2 ** (k - 2) * e_infp - 2 * e_3p) // 3
            K_Fs = 2 ** (2 * k - 3) * mup - 3 * 2 ** (k - 2) * e_infp
            if r % 4 == 1:
                ret = [K_Fw]
            elif r % 8 == 3:
                ret = [K_Fns, K_Fns, K_Fw]
            elif r % 8 == 7:
                ret = [K_Fs, K_Fs, K_Fw]
        return [int(x) for x in ret]

    def Fg_Fg(self):
        """
        Computes F_g^2 for each F in the strict transform of F°.

        Returns
        -------
        list
            The self intersection F_g^2 for each F_g in the strict transform
            of F°. The order of the components is given as follows:
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
        X = Xg_plus(N, r)
        Kdot = self.K_Fg()
        genera = X.genus()
        return [2 * genera[i] - 2 - Kdot[i] for i in range(0, len(Kdot))]

    def K_Fm(self, m):
        """
        Computes K.~F_m

        Parameters
        ----------
        m : int
            An integer `m` coprime to `self.N` for which we are intersecting
            with the HZ divisor ~F_`m`

        Returns
        -------
        int
            The intersection number K.~F_`m` on the surface.
        """
        N = self.N
        r = self.r
        if gcd(m, N) != 1:
            raise ValueError("m should be coprime to N")
        if not (m * r) % N in [
            x**2 % N for x in range(1, N) if gcd(x, N) == 1
        ]:  # should change this is O(N)
            raise ValueError("m*r is not a square mod N")
        X = X_0(m)
        mu = X.mu
        nu_inf = X.e_inf
        nu_2 = X.e_2
        nu_3 = X.e_3
        return ((mu - 3 * nu_inf - nu_3) // 3) - err(N, m)[0]

    def sketch_cusps(
        self, tex=True, compile_tex=False, open_pdf=False, disp=False, delay=1
    ):
        """
        Sketch intersection behaviour at each cusp of ~Z_(N,r).

        Parameters
        ----------
        tex : bool, default=True
            Whether tex code is generated in ./figs/N-r/
        compile_tex : bool, default=False
            Whether the tex code is then compiled (may cause error on Windows
            machines?)
        open_pdf : bool, default=False
            If tex is generated display pdf
        disp : bool, default=False
            Display matplotlib figures
        delay=1 : bool, default=False
            Delay in secs between figure displays

        Returns
        -------
        list
            A list of matplotlib.pyplot figures indexed by the cusps.
        """
        N = self.N
        r = self.r
        qq = [q for q in range(1, N) if rho(N, q * r) != 0]
        ret = []
        for d in [int(d) for d in divisors(N) if d != 1]:
            sing = list(set([q % d for q in qq]))
            sing.sort()
            for q in sing:
                c = ZtilCusp(N, r, d, q)
                here = c.sketch()
                if disp:
                    plt.show(block=False)
                    plt.pause(delay)
                here.suptitle("")
                ret.append((here, d, q, c.title))
        if tex:
            dir_name = f"./figs/Z{N}-{r}/"
            f_names = []
            caps = []
            try:
                os.mkdir(dir_name)
            except FileExistsError:
                pass
            for c_fig, d, q, caption in ret:
                filename = f"{d}-{q}.pgf"
                c_fig.savefig(dir_name + filename)
                f_names.append(filename)
                caps.append(caption)
            make_tex_file(dir_name, f_names, captions=caps)
            if compile_tex:
                compile_tex_file(dir_name, open_pdf=open_pdf)
        return [c_fig for c_fig, _, _, _ in ret]
