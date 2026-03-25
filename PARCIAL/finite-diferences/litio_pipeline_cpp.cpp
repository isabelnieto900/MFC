#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

    constexpr int NPT = 16000;
    constexpr int NBMAX = 1000;
    constexpr double HARTREE_TO_EV = 27.211386245988;
    constexpr int ITMAX_ZBRENT = 100;
    constexpr double EPS_ZBRENT = 3.0e-8;

    struct SpecParams
    {
        int n = NPT;
        double h = 120.0 / static_cast<double>(NPT);
        double ell = 0.0;
        double z = 3.0;
        double alpha = 2.535930;
    } g_spec;

    double func_vec(const std::vector<double> &x, double z)
    {
        const double a = x[0];
        const double b = x[1];
        const double y = 2.0 * a * b / std::pow(2.0 * a + b, 5);

        return a * a - 2.0 * z * a + (5.0 / 8.0) * a + (1.0 / 8.0) * b * b - (z * b / 4.0) +
               y * (8.0 * std::pow(a, 4) + 20.0 * std::pow(a, 3) * b + 12.0 * a * a * b * b +
                    10.0 * a * std::pow(b, 3) + std::pow(b, 4));
    }

    double ionization_energy_li(double alpha, double beta, double z)
    {
        const double y = 2.0 * alpha * beta / std::pow(2.0 * alpha + beta, 5);
        return beta * beta / 8.0 - z * beta / 4.0 +
               y * (8.0 * std::pow(alpha, 4) + 20.0 * std::pow(alpha, 3) * beta +
                    12.0 * alpha * alpha * beta * beta + 10.0 * alpha * std::pow(beta, 3) +
                    std::pow(beta, 4));
    }

    void set_spec_params(double z, double alpha, double ell, int n)
    {
        g_spec.n = n;
        g_spec.h = 120.0 / static_cast<double>(n);
        g_spec.ell = ell;
        g_spec.z = z;
        g_spec.alpha = alpha;
    }

    double pot_l(double x)
    {
        double p = -g_spec.z / x + 2.0 / x * (1.0 - (1.0 + g_spec.alpha * x) * std::exp(-2.0 * g_spec.alpha * x));
        p += g_spec.ell * (g_spec.ell + 1.0) / (2.0 * x * x);
        return p;
    }

    double determ_l(double energy)
    {
        const double h2 = g_spec.h * g_spec.h;
        double p0 = 0.0;
        double p1 = 2.0 + 2.0 * h2 * pot_l(g_spec.h) - 2.0 * energy * h2;
        double p2 = p1;

        for (int i = 2; i <= g_spec.n; ++i)
        {
            const double x = g_spec.h * static_cast<double>(i);
            p2 = (2.0 + 2.0 * h2 * pot_l(x) - 2.0 * h2 * energy) * p1 - p0;
            p0 = p1;
            p1 = p2;
        }

        return p2;
    }

    void zbrak(const std::function<double(double)> &fx, double x1, double x2, int n,
               std::vector<double> &xb1, std::vector<double> &xb2, int &nb)
    {
        int nbb = 0;
        const int nbmax = nb;
        const double dx = (x2 - x1) / static_cast<double>(n);
        double x = x1;
        double fp = fx(x);

        for (int i = 1; i <= n; ++i)
        {
            x += dx;
            const double fc = fx(x);
            if (fc * fp <= 0.0)
            {
                if (nbb >= static_cast<int>(xb1.size()) || nbb >= static_cast<int>(xb2.size()))
                {
                    break;
                }
                xb1[nbb] = x - dx;
                xb2[nbb] = x;
                ++nbb;
                if (nbb == nbmax)
                {
                    nb = nbb;
                    return;
                }
            }
            fp = fc;
        }

        nb = nbb;
    }

    double zbrent(const std::function<double(double)> &func, double x1, double x2, double tol)
    {
        double a = x1;
        double b = x2;
        double c = x2;
        double d = 0.0;
        double e = 0.0;

        double fa = func(a);
        double fb = func(b);
        if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
        {
            throw std::runtime_error("Root must be bracketed in zbrent");
        }

        double fc = fb;

        for (int iter = 0; iter < ITMAX_ZBRENT; ++iter)
        {
            if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
            {
                c = a;
                fc = fa;
                d = b - a;
                e = d;
            }

            if (std::abs(fc) < std::abs(fb))
            {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }

            const double tol1 = 2.0 * EPS_ZBRENT * std::abs(b) + 0.5 * tol;
            const double xm = 0.5 * (c - b);

            if (std::abs(xm) <= tol1 || fb == 0.0)
            {
                return b;
            }

            if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb))
            {
                const double s = fb / fa;
                double p, q;
                if (a == c)
                {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                }
                else
                {
                    const double q1 = fa / fc;
                    const double r = fb / fc;
                    p = s * (2.0 * xm * q1 * (q1 - r) - (b - a) * (r - 1.0));
                    q = (q1 - 1.0) * (r - 1.0) * (s - 1.0);
                }

                if (p > 0.0)
                    q = -q;
                p = std::abs(p);
                const double min1 = 3.0 * xm * q - std::abs(tol1 * q);
                const double min2 = std::abs(e * q);

                if (2.0 * p < std::min(min1, min2))
                {
                    e = d;
                    d = p / q;
                }
                else
                {
                    d = xm;
                    e = d;
                }
            }
            else
            {
                d = xm;
                e = d;
            }

            a = b;
            fa = fb;

            if (std::abs(d) > tol1)
            {
                b += d;
            }
            else
            {
                b += std::copysign(tol1, xm);
            }

            fb = func(b);
        }

        throw std::runtime_error("Maximum iterations in zbrent");
    }

    void tridag(const std::vector<double> &a, const std::vector<double> &b, const std::vector<double> &c,
                const std::vector<double> &r, std::vector<double> &u, int n)
    {
        std::vector<double> gam(n, 0.0);
        if (b[0] == 0.0)
        {
            throw std::runtime_error("Error 1 in tridag");
        }

        double bet = b[0];
        u[0] = r[0] / bet;

        for (int j = 1; j < n; ++j)
        {
            gam[j] = c[j - 1] / bet;
            bet = b[j] - a[j] * gam[j];
            if (bet == 0.0)
            {
                throw std::runtime_error("Error 2 in tridag");
            }
            u[j] = (r[j] - a[j] * u[j - 1]) / bet;
        }

        for (int j = n - 2; j >= 0; --j)
        {
            u[j] -= gam[j + 1] * u[j + 1];
        }
    }

    void brent_min(const std::function<double(double)> &f, double ax, double bx, double cx,
                   double tol, double &xmin, double &fmin)
    {
        constexpr int itmax = 200;
        constexpr double cgold = 0.3819660;
        constexpr double zeps = 1.0e-12;

        double a = std::min(ax, cx);
        double b = std::max(ax, cx);
        double x = bx;
        double w = bx;
        double v = bx;
        double fw = f(w);
        double fv = fw;
        double fx = fw;
        double d = 0.0;
        double e = 0.0;

        for (int iter = 0; iter < itmax; ++iter)
        {
            const double xm = 0.5 * (a + b);
            const double tol1 = tol * std::abs(x) + zeps;
            const double tol2 = 2.0 * tol1;

            if (std::abs(x - xm) <= (tol2 - 0.5 * (b - a)))
            {
                xmin = x;
                fmin = fx;
                return;
            }

            if (std::abs(e) > tol1)
            {
                const double r = (x - w) * (fx - fv);
                const double q = (x - v) * (fx - fw);
                double p = (x - v) * q - (x - w) * r;
                double q2 = 2.0 * (q - r);
                if (q2 > 0.0)
                    p = -p;
                q2 = std::abs(q2);
                const double etemp = e;
                e = d;

                if (std::abs(p) >= std::abs(0.5 * q2 * etemp) || p <= q2 * (a - x) || p >= q2 * (b - x))
                {
                    e = (x >= xm) ? (a - x) : (b - x);
                    d = cgold * e;
                }
                else
                {
                    d = p / q2;
                    const double u = x + d;
                    if (u - a < tol2 || b - u < tol2)
                    {
                        d = std::copysign(tol1, xm - x);
                    }
                }
            }
            else
            {
                e = (x >= xm) ? (a - x) : (b - x);
                d = cgold * e;
            }

            const double u = (std::abs(d) >= tol1) ? (x + d) : (x + std::copysign(tol1, d));
            const double fu = f(u);

            if (fu <= fx)
            {
                if (u >= x)
                    a = x;
                else
                    b = x;
                v = w;
                fv = fw;
                w = x;
                fw = fx;
                x = u;
                fx = fu;
            }
            else
            {
                if (u < x)
                    a = u;
                else
                    b = u;
                if (fu <= fw || w == x)
                {
                    v = w;
                    fv = fw;
                    w = u;
                    fw = fu;
                }
                else if (fu <= fv || v == x || v == w)
                {
                    v = u;
                    fv = fu;
                }
            }
        }

        xmin = x;
        fmin = fx;
    }

    void linmin(std::vector<double> &p, std::vector<double> &xi,
                const std::function<double(const std::vector<double> &)> &func, double &fret)
    {
        const std::vector<double> p0 = p;
        const std::vector<double> dir = xi;

        auto f1d = [&](double alpha)
        {
            std::vector<double> xt(p.size());
            for (size_t i = 0; i < p.size(); ++i)
            {
                xt[i] = p0[i] + alpha * dir[i];
            }
            return func(xt);
        };

        constexpr double gold = 1.618034;
        double ax = 0.0;
        double bx = 1.0;
        double fa = f1d(ax);
        double fb = f1d(bx);

        if (fb > fa)
        {
            std::swap(ax, bx);
            std::swap(fa, fb);
        }

        double cx = bx + gold * (bx - ax);
        double fc = f1d(cx);

        for (int i = 0; i < 100; ++i)
        {
            if (fb <= fc)
                break;
            ax = bx;
            fa = fb;
            bx = cx;
            fb = fc;
            cx = bx + gold * (bx - ax);
            fc = f1d(cx);
        }

        double xmin, fmin;
        brent_min(f1d, ax, bx, cx, 1.0e-10, xmin, fmin);

        for (size_t i = 0; i < p.size(); ++i)
        {
            p[i] = p0[i] + xmin * dir[i];
            xi[i] = xmin * dir[i];
        }

        fret = fmin;
    }

    void powell(std::vector<double> &p, std::vector<std::vector<double>> &xi, int n,
                double ftol, int &iter, double &fret,
                const std::function<double(const std::vector<double> &)> &func)
    {
        constexpr int itmax = 200;
        constexpr double tiny = 1.0e-25;

        std::vector<double> pt = p;
        std::vector<double> ptt(n);
        std::vector<double> xit(n);

        fret = func(p);

        for (iter = 1; iter <= itmax; ++iter)
        {
            const double fp = fret;
            int ibig = 0;
            double del = 0.0;

            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    xit[j] = xi[j][i];
                }
                const double fptt = fret;
                linmin(p, xit, func, fret);
                if (fptt - fret > del)
                {
                    del = fptt - fret;
                    ibig = i;
                }
            }

            if (2.0 * (fp - fret) <= ftol * (std::abs(fp) + std::abs(fret)) + tiny)
            {
                return;
            }

            for (int j = 0; j < n; ++j)
            {
                ptt[j] = 2.0 * p[j] - pt[j];
                xit[j] = p[j] - pt[j];
                pt[j] = p[j];
            }

            const double fptt = func(ptt);

            if (fptt < fp)
            {
                const double t = 2.0 * (fp - 2.0 * fret + fptt) * std::pow(fp - fret - del, 2) -
                                 del * std::pow(fp - fptt, 2);
                if (t < 0.0)
                {
                    linmin(p, xit, func, fret);
                    for (int j = 0; j < n; ++j)
                    {
                        xi[j][ibig] = xi[j][n - 1];
                        xi[j][n - 1] = xit[j];
                    }
                }
            }
        }

        throw std::runtime_error("powell exceeding maximum iterations");
    }

    void build_potential_table(double z, double alpha, int n, std::vector<double> &r,
                               std::vector<double> &vc, std::vector<double> &veff)
    {
        const double rmin = 1.0e-4;
        const double rmax = 120.0;
        const double h = (rmax - rmin) / static_cast<double>(n - 1);

        for (int i = 0; i < n; ++i)
        {
            const double rr = rmin + static_cast<double>(i) * h;
            r[i] = rr;
            vc[i] = -(z - 2.0) / rr;
            veff[i] = -z / rr + 2.0 / rr * (1.0 - (1.0 + alpha * rr) * std::exp(-2.0 * alpha * rr));
        }
    }

    void solve_l_spectrum(double z, double alpha, double ell, int n, double x1, double x2,
                          int nscan, std::vector<double> &energies, int &nfound)
    {
        set_spec_params(z, alpha, ell, n);

        std::vector<double> xb1(NBMAX, 0.0), xb2(NBMAX, 0.0);
        int nb = NBMAX;
        zbrak(determ_l, x1, x2, nscan, xb1, xb2, nb);

        nfound = nb;
        std::fill(energies.begin(), energies.end(), 0.0);

        for (int i = 0; i < nb; ++i)
        {
            double tol = 1.0e-12 * std::abs(0.5 * (xb1[i] + xb2[i]));
            if (tol == 0.0)
                tol = 1.0e-12;
            energies[i] = zbrent(determ_l, xb1[i], xb2[i], tol);
        }
    }

    void wave_for_energy(double z, double alpha, double ell, int n, double energy,
                         std::vector<double> &wave)
    {
        set_spec_params(z, alpha, ell, n);

        const double h = g_spec.h;
        std::vector<double> c(n), diag(n), superd(n), subd(n), vk(n, 0.0), vkm1(n, 0.0);

        for (int k = 0; k < n; ++k)
        {
            c[k] = pot_l(h * static_cast<double>(k + 1));
            diag[k] = 2.0 + 2.0 * h * h * c[k] - 2.0 * h * h * energy;
        }

        for (int k = 0; k < n - 1; ++k)
        {
            superd[k] = -1.0;
        }
        superd[n - 1] = 0.0;

        subd[0] = 0.0;
        for (int k = 1; k < n; ++k)
        {
            subd[k] = -1.0;
        }

        vk[0] = 1.0;

        for (int kk = 0; kk < 10; ++kk)
        {
            tridag(subd, diag, superd, vk, vkm1, n);
            double s = 0.0;
            for (int i = 0; i < n; ++i)
                s += vkm1[i] * vkm1[i];
            const double nrm = std::sqrt(s * h);
            if (nrm > 0.0)
            {
                for (int i = 0; i < n; ++i)
                    vk[i] = vkm1[i] / nrm;
            }
            else
            {
                vk = vkm1;
            }
        }

        wave = vk;
    }

    double rel_error(double value, double reference)
    {
        return std::abs((value - reference) / reference);
    }

} // namespace

int main(int argc, char **argv)
{
    double z = 3.0;
    if (argc > 1)
    {
        z = std::stod(argv[1]);
    }

    auto t0 = std::chrono::high_resolution_clock::now();

    const double alpha0 = 1.5;
    const double beta0 = 1.5;

    std::vector<double> p = {alpha0, beta0};
    std::vector<std::vector<double>> xi(2, std::vector<double>(2, 0.0));
    xi[0][0] = 1.0;
    xi[1][1] = 1.0;

    int it_var = 0;
    double e0 = 0.0;
    powell(p, xi, 2, 1.0e-8, it_var, e0, [&](const std::vector<double> &x)
           { return func_vec(x, z); });

    const double alpha_opt = p[0];
    const double beta_opt = p[1];
    const double ei = ionization_energy_li(alpha_opt, beta_opt, z);

    std::vector<double> r(NPT), vc(NPT), veff(NPT);
    std::vector<double> rs2(NPT), rs3(NPT), rs4(NPT), rp2(NPT), rp3(NPT);
    std::vector<double> e_s(NBMAX, 0.0), e_p(NBMAX, 0.0);
    int nfound_s = 0, nfound_p = 0;

    build_potential_table(z, 2.535930, NPT, r, vc, veff);
    solve_l_spectrum(z, 2.535930, 0.0, NPT, -10.0, 0.0, 2500, e_s, nfound_s);
    solve_l_spectrum(z, 2.535930, 1.0, NPT, -10.0, 0.0, 2500, e_p, nfound_p);

    if (nfound_s < 4 || nfound_p < 2)
    {
        std::cerr << "No se encontraron estados suficientes" << std::endl;
        return 1;
    }

    wave_for_energy(z, 2.535930, 0.0, NPT, e_s[1], rs2);
    wave_for_energy(z, 2.535930, 0.0, NPT, e_s[2], rs3);
    wave_for_energy(z, 2.535930, 0.0, NPT, e_s[3], rs4);
    wave_for_energy(z, 2.535930, 1.0, NPT, e_p[0], rp2);
    wave_for_energy(z, 2.535930, 1.0, NPT, e_p[1], rp3);

    auto t1 = std::chrono::high_resolution_clock::now();
    const double elapsed = std::chrono::duration<double>(t1 - t0).count();

    {
        std::ofstream out("potential_litio_cpp.dat");
        out << "# r(a.u.)  Vcoulomb(Ha)  Veff(Ha)\n";
        out << std::setprecision(14) << std::scientific;
        for (int i = 0; i < NPT; ++i)
        {
            out << std::setw(22) << r[i] << ' ' << std::setw(22) << vc[i] << ' ' << std::setw(22) << veff[i] << '\n';
        }
    }

    {
        std::ofstream out("funciones_radiales_litio_cpp.dat");
        out << "# r(a.u.)  R_2s  R_3s  R_4s  R_2p  R_3p\n";
        out << std::setprecision(14) << std::scientific;
        for (int i = 0; i < NPT; ++i)
        {
            out << std::setw(22) << r[i] << ' ' << std::setw(22) << rs2[i] << ' ' << std::setw(22) << rs3[i] << ' '
                << std::setw(22) << rs4[i] << ' ' << std::setw(22) << rp2[i] << ' ' << std::setw(22) << rp3[i] << '\n';
        }
    }

    const double e2s_ev = e_s[1] * HARTREE_TO_EV;
    const double e3s_ev = e_s[2] * HARTREE_TO_EV;
    const double e4s_ev = e_s[3] * HARTREE_TO_EV;
    const double e2p_ev = e_p[0] * HARTREE_TO_EV;
    const double e3p_ev = e_p[1] * HARTREE_TO_EV;

    {
        std::ofstream out("energias_litio_cpp.dat");
        out << "# estado  l  E(Ha)  E(eV)\n";
        out << std::setprecision(14) << std::scientific;
        out << "2s 0 " << e_s[1] << ' ' << e2s_ev << '\n';
        out << "3s 0 " << e_s[2] << ' ' << e3s_ev << '\n';
        out << "4s 0 " << e_s[3] << ' ' << e4s_ev << '\n';
        out << "2p 1 " << e_p[0] << ' ' << e2p_ev << '\n';
        out << "3p 1 " << e_p[1] << ' ' << e3p_ev << '\n';
    }

    const double ei_ev = ei * HARTREE_TO_EV;
    const double e0_ev = e0 * HARTREE_TO_EV;

    const double err_rel_2s = rel_error(e2s_ev, -5.39);
    const double err_rel_ei = rel_error(ei_ev, -4.58064);

    {
        std::ofstream out("metrics_litio_cpp.dat");
        out << "# language time_s var_iter e0_ev ion_ev e2s_ev err_rel_ion err_rel_2s nfound_s nfound_p\n";
        out << std::setprecision(10) << std::fixed;
        out << "cpp " << elapsed << ' ' << it_var << ' ' << e0_ev << ' ' << ei_ev << ' ' << e2s_ev << ' '
            << err_rel_ei << ' ' << err_rel_2s << ' ' << nfound_s << ' ' << nfound_p << '\n';
    }

    std::cout << "Iteraciones (variacional): " << it_var << "\n";
    std::cout << "E0 [eV]: " << e0_ev << "\n";
    std::cout << "Eion [eV]: " << ei_ev << "\n";
    std::cout << "E2s [eV]: " << e2s_ev << "\n";
    std::cout << "Tiempo [s]: " << elapsed << "\n";

    return 0;
}
