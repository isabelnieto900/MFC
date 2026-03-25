#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>

using namespace std;

// Tipos para manejar funciones y vectores
typedef function<double(const vector<double>&)> FuncND;

namespace Optimization {
    // Variables globales para la comunicación interna de linmin
    vector<double> pcom, xicom;
    FuncND nrfunc;

    double f1dim(double x) {
        int n = pcom.size();
        vector<double> xt(n);
        for (int j = 0; j < n; j++) xt[j] = pcom[j] + x * xicom[j];
        return nrfunc(xt);
    }

    void mnbrak(double &ax, double &bx, double &cx, double &fa, double &fb, double &fc) {
        const double GOLD = 1.618034, GLIMIT = 100.0, TINY = 1.0e-20;
        fa = f1dim(ax); fb = f1dim(bx);
        if (fb > fa) {
            swap(ax, bx); swap(fa, fb);
        }
        cx = bx + GOLD * (bx - ax);
        fc = f1dim(cx);
        while (fb > fc) {
            double r = (bx - ax) * (fb - fc);
            double q = (bx - cx) * (fb - fa);
            double u = bx - ((bx - cx) * q - (bx - ax) * r) / (2.0 * max(abs(q - r), TINY) * (q - r >= 0 ? 1 : -1));
            double ulim = bx + GLIMIT * (cx - bx);
            double fu;
            if ((bx - u) * (u - cx) > 0.0) {
                fu = f1dim(u);
                if (fu < fc) {
                    ax = bx; bx = u; fa = fb; fb = fu; return;
                } else if (fu > fb) {
                    cx = u; fc = fu; return;
                }
                u = cx + GOLD * (cx - bx); fu = f1dim(u);
            } else if ((cx - u) * (u - ulim) > 0.0) {
                fu = f1dim(u);
                if (fu < fc) {
                    bx = cx; cx = u; u = cx + GOLD * (cx - bx);
                    fb = fc; fc = fu; fu = f1dim(u);
                }
            } else if ((u - ulim) * (ulim - cx) >= 0.0) {
                u = ulim; fu = f1dim(u);
            } else {
                u = cx + GOLD * (cx - bx); fu = f1dim(u);
            }
            ax = bx; bx = cx; cx = u; fa = fb; fb = fc; fc = fu;
        }
    }

    double brent(double ax, double bx, double cx, double tol, double &xmin) {
        const int ITMAX = 100;
        const double CGOLD = 0.3819660, ZEPS = 1.0e-10;
        double a = (ax < cx ? ax : cx), b = (ax > cx ? ax : cx);
        double x = bx, w = bx, v = bx, e = 0.0, d = 0.0;
        double fx = f1dim(x), fw = fx, fv = fx;

        for (int iter = 1; iter <= ITMAX; iter++) {
            double xm = 0.5 * (a + b);
            double tol1 = tol * abs(x) + ZEPS, tol2 = 2.0 * tol1;
            if (abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
                xmin = x; return fx;
            }
            if (abs(e) > tol1) {
                double r = (x - w) * (fx - fv), q = (x - v) * (fx - fw);
                double p = (x - v) * q - (x - w) * r;
                q = 2.0 * (q - r);
                if (q > 0.0) p = -p;
                q = abs(q);
                double etemp = e; e = d;
                if (abs(p) >= abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
                    e = (x >= xm ? a - x : b - x); d = CGOLD * e;
                } else {
                    d = p / q; double u = x + d;
                    if (u - a < tol2 || b - u < tol2) d = (xm - x >= 0 ? abs(tol1) : -abs(tol1));
                }
            } else {
                e = (x >= xm ? a - x : b - x); d = CGOLD * e;
            }
            double u = (abs(d) >= tol1 ? x + d : x + (d >= 0 ? abs(tol1) : -abs(tol1)));
            double fu = f1dim(u);
            if (fu <= fx) {
                if (u >= x) a = x; else b = x;
                v = w; w = x; x = u; fv = fw; fw = fx; fx = fu;
            } else {
                if (u < x) a = u; else b = u;
                if (fu <= fw || w == x) {
                    v = w; w = u; fv = fw; fw = fu;
                } else if (fu <= fv || v == x || v == w) {
                    v = u; fv = fu;
                }
            }
        }
        xmin = x; return fx;
    }

    void linmin(vector<double> &p, vector<double> &xi, double &fret, FuncND func) {
        int n = p.size();
        nrfunc = func; pcom = p; xicom = xi;
        double ax = 0.0, xx = 1.0, cx, fa, fx, fc, xmin;
        mnbrak(ax, xx, cx, fa, fx, fc);
        fret = brent(ax, xx, cx, 1.0e-8, xmin);
        for (int j = 0; j < n; j++) {
            xi[j] *= xmin;
            p[j] += xi[j];
        }
    }

    void powell(vector<double> &p, vector<vector<double>> &xi, double ftol, int &iter, double &fret, FuncND func) {
        const int ITMAX = 200;
        int n = p.size();
        fret = func(p);
        vector<double> pt = p, ptt(n), xit(n);
        iter = 0;
        while (iter < ITMAX) {
            iter++;
            double fp = fret; int ibig = 0; double del = 0.0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) xit[j] = xi[j][i];
                double fptt = fret;
                linmin(p, xit, fret, func);
                if (abs(fptt - fret) > del) {
                    del = abs(fptt - fret); ibig = i + 1;
                }
            }
            if (2.0 * abs(fp - fret) <= ftol * (abs(fp) + abs(fret))) return;
            for (int j = 0; j < n; j++) {
                ptt[j] = 2.0 * p[j] - pt[j];
                xit[j] = p[j] - pt[j];
                pt[j] = p[j];
            }
            double fptt = func(ptt);
            if (fptt < fp) {
                double t = 2.0 * (fp - 2.0 * fret + fptt) * pow(fp - fret - del, 2) - del * pow(fp - fptt, 2);
                if (t < 0.0) {
                    linmin(p, xit, fret, func);
                    for (int j = 0; j < n; j++) {
                        xi[j][ibig - 1] = xi[j][n - 1];
                        xi[j][n - 1] = xit[j];
                    }
                }
            }
        }
    }
}