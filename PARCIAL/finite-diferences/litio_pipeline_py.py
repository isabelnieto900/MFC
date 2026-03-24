import math
import sys
import time
from pathlib import Path

import numpy as np

NPT = 2200
NBMAX = 1000
HARTREE_TO_EV = 27.211386245988
ITMAX_ZBRENT = 100
EPS_ZBRENT = 3.0e-8

spec = {
    "n": NPT,
    "h": 120.0 / NPT,
    "ell": 0.0,
    "z": 3.0,
    "alpha": 2.535930,
}


def func_vec(x: np.ndarray, z: float) -> float:
    a = x[0]
    b = x[1]
    y = 2.0 * a * b / (2.0 * a + b) ** 5
    return (
        a * a
        - 2.0 * z * a
        + (5.0 / 8.0) * a
        + (1.0 / 8.0) * b * b
        - (z * b / 4.0)
        + y * (8.0 * a**4 + 20.0 * a**3 * b + 12.0 * a * a * b * b + 10.0 * a * b**3 + b**4)
    )


def ionization_energy_li(alpha: float, beta: float, z: float) -> float:
    y = 2.0 * alpha * beta / (2.0 * alpha + beta) ** 5
    return (
        beta * beta / 8.0
        - z * beta / 4.0
        + y * (8.0 * alpha**4 + 20.0 * alpha**3 * beta + 12.0 * alpha * alpha * beta * beta + 10.0 * alpha * beta**3 + beta**4)
    )


def set_spec_params(z: float, alpha: float, ell: float, n: int) -> None:
    spec["n"] = n
    spec["h"] = 120.0 / n
    spec["ell"] = ell
    spec["z"] = z
    spec["alpha"] = alpha


def pot_l(x: float) -> float:
    p = -spec["z"] / x + 2.0 / x * (1.0 - (1.0 + spec["alpha"] * x) * math.exp(-2.0 * spec["alpha"] * x))
    p += spec["ell"] * (spec["ell"] + 1.0) / (2.0 * x * x)
    return p


def determ_l(energy: float) -> float:
    h2 = spec["h"] * spec["h"]
    p0 = 0.0
    p1 = 2.0 + 2.0 * h2 * pot_l(spec["h"]) - 2.0 * energy * h2
    p2 = p1

    for i in range(2, spec["n"] + 1):
        x = spec["h"] * i
        p2 = (2.0 + 2.0 * h2 * pot_l(x) - 2.0 * h2 * energy) * p1 - p0
        p0 = p1
        p1 = p2

    return p2


def zbrak(fx, x1: float, x2: float, n: int, nbmax: int):
    xb1 = []
    xb2 = []
    dx = (x2 - x1) / n
    x = x1
    fp = fx(x)

    for _ in range(n):
        x += dx
        fc = fx(x)
        if fc * fp <= 0.0:
            xb1.append(x - dx)
            xb2.append(x)
            if len(xb1) == nbmax:
                break
        fp = fc

    return xb1, xb2


def zbrent(func, x1: float, x2: float, tol: float) -> float:
    a = x1
    b = x2
    c = x2
    d = 0.0
    e = 0.0

    fa = func(a)
    fb = func(b)
    if (fa > 0.0 and fb > 0.0) or (fa < 0.0 and fb < 0.0):
        raise RuntimeError("Root must be bracketed in zbrent")

    fc = fb

    for _ in range(ITMAX_ZBRENT):
        if (fb > 0.0 and fc > 0.0) or (fb < 0.0 and fc < 0.0):
            c = a
            fc = fa
            d = b - a
            e = d

        if abs(fc) < abs(fb):
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa

        tol1 = 2.0 * EPS_ZBRENT * abs(b) + 0.5 * tol
        xm = 0.5 * (c - b)

        if abs(xm) <= tol1 or fb == 0.0:
            return b

        if abs(e) >= tol1 and abs(fa) > abs(fb):
            s = fb / fa
            if a == c:
                p = 2.0 * xm * s
                q = 1.0 - s
            else:
                q1 = fa / fc
                r = fb / fc
                p = s * (2.0 * xm * q1 * (q1 - r) - (b - a) * (r - 1.0))
                q = (q1 - 1.0) * (r - 1.0) * (s - 1.0)

            if p > 0.0:
                q = -q
            p = abs(p)
            min1 = 3.0 * xm * q - abs(tol1 * q)
            min2 = abs(e * q)

            if 2.0 * p < min(min1, min2):
                e = d
                d = p / q
            else:
                d = xm
                e = d
        else:
            d = xm
            e = d

        a = b
        fa = fb

        if abs(d) > tol1:
            b += d
        else:
            b += math.copysign(tol1, xm)

        fb = func(b)

    raise RuntimeError("Maximum iterations in zbrent")


def tridag(a: np.ndarray, b: np.ndarray, c: np.ndarray, r: np.ndarray) -> np.ndarray:
    n = len(a)
    gam = np.zeros(n)
    u = np.zeros(n)

    if b[0] == 0.0:
        raise RuntimeError("Error 1 in tridag")

    bet = b[0]
    u[0] = r[0] / bet

    for j in range(1, n):
        gam[j] = c[j - 1] / bet
        bet = b[j] - a[j] * gam[j]
        if bet == 0.0:
            raise RuntimeError("Error 2 in tridag")
        u[j] = (r[j] - a[j] * u[j - 1]) / bet

    for j in range(n - 2, -1, -1):
        u[j] = u[j] - gam[j + 1] * u[j + 1]

    return u


def brent_min(f1d, ax: float, bx: float, cx: float, tol: float):
    itmax = 200
    cgold = 0.3819660
    zeps = 1.0e-12

    a = min(ax, cx)
    b = max(ax, cx)
    x = bx
    w = bx
    v = bx
    fw = f1d(w)
    fv = fw
    fx = fw
    d = 0.0
    e = 0.0

    for _ in range(itmax):
        xm = 0.5 * (a + b)
        tol1 = tol * abs(x) + zeps
        tol2 = 2.0 * tol1

        if abs(x - xm) <= (tol2 - 0.5 * (b - a)):
            return x, fx

        if abs(e) > tol1:
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q2 = 2.0 * (q - r)
            if q2 > 0.0:
                p = -p
            q2 = abs(q2)
            etemp = e
            e = d

            if abs(p) >= abs(0.5 * q2 * etemp) or p <= q2 * (a - x) or p >= q2 * (b - x):
                e = a - x if x >= xm else b - x
                d = cgold * e
            else:
                d = p / q2
                u = x + d
                if u - a < tol2 or b - u < tol2:
                    d = math.copysign(tol1, xm - x)
        else:
            e = a - x if x >= xm else b - x
            d = cgold * e

        u = x + d if abs(d) >= tol1 else x + math.copysign(tol1, d)
        fu = f1d(u)

        if fu <= fx:
            if u >= x:
                a = x
            else:
                b = x
            v, fv = w, fw
            w, fw = x, fx
            x, fx = u, fu
        else:
            if u < x:
                a = u
            else:
                b = u
            if fu <= fw or w == x:
                v, fv = w, fw
                w, fw = u, fu
            elif fu <= fv or v == x or v == w:
                v, fv = u, fu

    return x, fx


def linmin(p: np.ndarray, xi: np.ndarray, func):
    p0 = p.copy()
    direc = xi.copy()

    def f1d(alpha: float) -> float:
        return func(p0 + alpha * direc)

    gold = 1.618034
    ax = 0.0
    bx = 1.0
    fa = f1d(ax)
    fb = f1d(bx)

    if fb > fa:
        ax, bx = bx, ax
        fa, fb = fb, fa

    cx = bx + gold * (bx - ax)
    fc = f1d(cx)

    for _ in range(100):
        if fb <= fc:
            break
        ax, fa = bx, fb
        bx, fb = cx, fc
        cx = bx + gold * (bx - ax)
        fc = f1d(cx)

    xmin, fmin = brent_min(f1d, ax, bx, cx, 1.0e-10)
    p_new = p0 + xmin * direc
    xi_new = xmin * direc

    return p_new, xi_new, fmin


def powell(p: np.ndarray, xi: np.ndarray, ftol: float, func):
    itmax = 200
    tiny = 1.0e-25
    n = len(p)

    pt = p.copy()
    fret = func(p)

    for it in range(1, itmax + 1):
        fp = fret
        ibig = 0
        delta = 0.0

        for i in range(n):
            xit = xi[:, i].copy()
            fptt = fret
            p, xit, fret = linmin(p, xit, func)
            if fptt - fret > delta:
                delta = fptt - fret
                ibig = i

        if 2.0 * (fp - fret) <= ftol * (abs(fp) + abs(fret)) + tiny:
            return p, it, fret

        ptt = 2.0 * p - pt
        xit = p - pt
        pt = p.copy()
        fptt = func(ptt)

        if fptt < fp:
            t = 2.0 * (fp - 2.0 * fret + fptt) * (fp - fret - delta) ** 2 - delta * (fp - fptt) ** 2
            if t < 0.0:
                p, xit, fret = linmin(p, xit, func)
                xi[:, ibig] = xi[:, n - 1]
                xi[:, n - 1] = xit

    raise RuntimeError("powell exceeding maximum iterations")


def build_potential_table(z: float, alpha: float, n: int):
    rmin = 1.0e-4
    rmax = 120.0
    h = (rmax - rmin) / (n - 1)
    r = np.array([rmin + i * h for i in range(n)], dtype=float)
    vc = -(z - 2.0) / r
    veff = -z / r + 2.0 / r * (1.0 - (1.0 + alpha * r) * np.exp(-2.0 * alpha * r))
    return r, vc, veff


def solve_l_spectrum(z: float, alpha: float, ell: float, n: int, x1: float, x2: float, nscan: int):
    set_spec_params(z, alpha, ell, n)
    xb1, xb2 = zbrak(determ_l, x1, x2, nscan, NBMAX)
    energies = np.zeros(NBMAX)

    for i in range(len(xb1)):
        tol = 1.0e-12 * abs(0.5 * (xb1[i] + xb2[i]))
        if tol == 0.0:
            tol = 1.0e-12
        energies[i] = zbrent(determ_l, xb1[i], xb2[i], tol)

    return energies, len(xb1)


def wave_for_energy(z: float, alpha: float, ell: float, n: int, energy: float):
    set_spec_params(z, alpha, ell, n)
    h = spec["h"]

    c = np.zeros(n)
    diag = np.zeros(n)
    superd = np.zeros(n)
    subd = np.zeros(n)

    for k in range(n):
        c[k] = pot_l(h * (k + 1))
        diag[k] = 2.0 + 2.0 * h * h * c[k] - 2.0 * h * h * energy

    superd[:-1] = -1.0
    superd[-1] = 0.0
    subd[0] = 0.0
    subd[1:] = -1.0

    vk = np.zeros(n)
    vk[0] = 1.0

    for _ in range(10):
        vkm1 = tridag(subd, diag, superd, vk)
        nrm = math.sqrt(np.sum(vkm1 * vkm1) * h)
        if nrm > 0.0:
            vk = vkm1 / nrm
        else:
            vk = vkm1

    return vk


def rel_error(value: float, reference: float) -> float:
    return abs((value - reference) / reference)


def main() -> None:
    z = float(sys.argv[1]) if len(sys.argv) > 1 else 3.0
    t0 = time.perf_counter()

    p0 = np.array([1.5, 1.5], dtype=float)
    xi = np.eye(2)
    p_opt, it_var, e0 = powell(p0, xi, 1.0e-8, lambda x: func_vec(x, z))
    alpha_opt, beta_opt = p_opt[0], p_opt[1]
    ei = ionization_energy_li(alpha_opt, beta_opt, z)

    r, vc, veff = build_potential_table(z, 2.535930, NPT)
    e_s, nfound_s = solve_l_spectrum(z, 2.535930, 0.0, NPT, -10.0, 0.0, 2500)
    e_p, nfound_p = solve_l_spectrum(z, 2.535930, 1.0, NPT, -10.0, 0.0, 2500)

    if nfound_s < 4 or nfound_p < 2:
        raise RuntimeError("No se encontraron estados suficientes")

    rs2 = wave_for_energy(z, 2.535930, 0.0, NPT, e_s[1])
    rs3 = wave_for_energy(z, 2.535930, 0.0, NPT, e_s[2])
    rs4 = wave_for_energy(z, 2.535930, 0.0, NPT, e_s[3])
    rp2 = wave_for_energy(z, 2.535930, 1.0, NPT, e_p[0])
    rp3 = wave_for_energy(z, 2.535930, 1.0, NPT, e_p[1])

    elapsed = time.perf_counter() - t0

    np.savetxt(
        "potential_litio_py.dat",
        np.column_stack([r, vc, veff]),
        header="r(a.u.)  Vcoulomb(Ha)  Veff(Ha)",
        comments="# ",
        fmt="%.14e",
    )

    np.savetxt(
        "funciones_radiales_litio_py.dat",
        np.column_stack([r, rs2, rs3, rs4, rp2, rp3]),
        header="r(a.u.)  R_2s  R_3s  R_4s  R_2p  R_3p",
        comments="# ",
        fmt="%.14e",
    )

    e2s_ev = e_s[1] * HARTREE_TO_EV
    e3s_ev = e_s[2] * HARTREE_TO_EV
    e4s_ev = e_s[3] * HARTREE_TO_EV
    e2p_ev = e_p[0] * HARTREE_TO_EV
    e3p_ev = e_p[1] * HARTREE_TO_EV

    with Path("energias_litio_py.dat").open("w", encoding="utf-8") as f:
        f.write("# estado  l  E(Ha)  E(eV)\n")
        f.write(f"2s 0 {e_s[1]:.14e} {e2s_ev:.14e}\n")
        f.write(f"3s 0 {e_s[2]:.14e} {e3s_ev:.14e}\n")
        f.write(f"4s 0 {e_s[3]:.14e} {e4s_ev:.14e}\n")
        f.write(f"2p 1 {e_p[0]:.14e} {e2p_ev:.14e}\n")
        f.write(f"3p 1 {e_p[1]:.14e} {e3p_ev:.14e}\n")

    e0_ev = e0 * HARTREE_TO_EV
    ei_ev = ei * HARTREE_TO_EV

    err_rel_2s = rel_error(e2s_ev, -5.39)
    err_rel_ei = rel_error(ei_ev, -4.58064)

    with Path("metrics_litio_py.dat").open("w", encoding="utf-8") as f:
        f.write("# language time_s var_iter e0_ev ion_ev e2s_ev err_rel_ion err_rel_2s nfound_s nfound_p\n")
        f.write(
            f"py {elapsed:.10f} {it_var:d} {e0_ev:.10f} {ei_ev:.10f} {e2s_ev:.10f} "
            f"{err_rel_ei:.10f} {err_rel_2s:.10f} {nfound_s:d} {nfound_p:d}\n"
        )

    print(f"Iteraciones (variacional): {it_var}")
    print(f"E0 [eV]: {e0_ev:.10f}")
    print(f"Eion [eV]: {ei_ev:.10f}")
    print(f"E2s [eV]: {e2s_ev:.10f}")
    print(f"Tiempo [s]: {elapsed:.6f}")


if __name__ == "__main__":
    main()
