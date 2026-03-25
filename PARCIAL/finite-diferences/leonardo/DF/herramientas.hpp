#ifndef HERRAMIENTAS_HPP
#define HERRAMIENTAS_HPP

#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>

namespace Herramientas {
    using namespace std;

    // Prototipo de zbrak: busca intervalos donde existen raíces
    void zbrak(function<double(double)> fx, double x1, double x2, int n, vector<double>& xb1, vector<double>& xb2) {
        double dx = (x2 - x1) / n;
        double x = x1;
        double fp = fx(x);
        
        for (int i = 1; i <= n; i++) {
            x += dx;
            double fc = fx(x);
            if (fc * fp < 0.0) {
                xb1.push_back(x - dx);
                xb2.push_back(x);
            }
            fp = fc;
        }
    }

    // Prototipo de zbrent: encuentra la raíz exacta en un intervalo dado
    double zbrent(function<double(double)> fx, double x1, double x2, double tol) {
        const int ITMAX = 100;
        const double EPS = 3.0e-15;
        double a = x1, b = x2, c = x2, d, e, min1, min2;
        double fa = fx(a), fb = fx(b), fc, p, q, r, s, tol1, xm;

        if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) return 0.0;
        fc = fb;
        d = b - a;
        e = d;

        for (int iter = 1; iter <= ITMAX; iter++) {
            if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
                c = a; fc = fa; d = b - a; e = d;
            }
            if (abs(fc) < abs(fb)) {
                a = b; b = c; c = a;
                fa = fb; fb = fc; fc = fa;
            }
            tol1 = 2.0 * EPS * abs(b) + 0.5 * tol;
            xm = 0.5 * (c - b);
            if (abs(xm) <= tol1 || fb == 0.0) return b;
            if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
                s = fb / fa;
                if (a == c) {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                } else {
                    q = fa / fc;
                    r = fb / fc;
                    p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }
                if (p > 0.0) q = -q;
                p = abs(p);
                min1 = 3.0 * xm * q - abs(tol1 * q);
                min2 = abs(e * q);
                if (2.0 * p < (min1 < min2 ? min1 : min2)) {
                    e = d; d = p / q;
                } else {
                    d = xm; e = d;
                }
            } else {
                d = xm; e = d;
            }
            a = b; fa = fb;
            if (abs(d) > tol1) b += d;
            else b += (xm >= 0 ? abs(tol1) : -abs(tol1));
            fb = fx(b);
        }
        return b;
    }

    // Algoritmo Tridiagonal para resolver sistemas de ecuaciones lineales
    void tridag(const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& r, vector<double>& u) {
        int n = r.size();
        vector<double> gam(n);
        double bet = b[0];
        u[0] = r[0] / bet;
        for (int j = 1; j < n; j++) {
            gam[j] = c[j - 1] / bet;
            bet = b[j] - a[j] * gam[j];
            u[j] = (r[j] - a[j] * u[j - 1]) / bet;
        }
        for (int j = n - 2; j >= 0; j--) {
            u[j] -= gam[j + 1] * u[j + 1];
        }
    }

    // Normalización de vectores (función de onda)
    void norma(vector<double>& v) {
        double suma = 0.0;
        for (double val : v) suma += val * val;
        double modulo = sqrt(suma);
        for (double& val : v) val /= modulo;
    }
}

#endif