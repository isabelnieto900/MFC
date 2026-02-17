#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <utility>
using namespace std;
using namespace std::chrono;

const double b = 0.02, d = 0.015, r = 0.1, p0 = 0.01, t_max = 50.0;

inline double f(double p) { return r * b * (1 - p); }
inline double f_prime(double p) { return -r * b * f(p); }
inline double exacta(double t) { return 1 - (1 - p0) * exp(-r * b * t); }

pair<double, double> euler(double p0, double h, long n)
{
    double p = p0;
    auto start = high_resolution_clock::now();
    for (long i = 0; i < n; i++)
        p = p + h * f(p);
    auto end = high_resolution_clock::now();
    return {p, duration<double>(end - start).count()};
}

vector<double> euler_solucion(double p0, double h, long n)
{
    vector<double> p(n + 1);
    p[0] = p0;
    for (long i = 0; i < n; i++)
        p[i + 1] = p[i] + h * f(p[i]);
    return p;
}

pair<double, double> taylor2(double p0, double h, long n)
{
    double p = p0;
    auto start = high_resolution_clock::now();
    for (long i = 0; i < n; i++)
        p = p + h * f(p) + (h * h / 2.0) * f_prime(p);
    auto end = high_resolution_clock::now();
    return {p, duration<double>(end - start).count()};
}

vector<double> taylor2_solucion(double p0, double h, long n)
{
    vector<double> p(n + 1);
    p[0] = p0;
    for (long i = 0; i < n; i++)
        p[i + 1] = p[i] + h * f(p[i]) + (h * h / 2.0) * f_prime(p[i]);
    return p;
}

pair<double, double> trapecio(double p0, double h, long n)
{
    double p = p0;
    double factor = h * r * b / 2.0;
    auto start = high_resolution_clock::now();
    for (long i = 0; i < n; i++)
        p = (p + (h / 2.0) * f(p) + factor) / (1.0 + factor);
    auto end = high_resolution_clock::now();
    return {p, duration<double>(end - start).count()};
}

vector<double> trapecio_solucion(double p0, double h, long n)
{
    vector<double> p(n + 1);
    p[0] = p0;
    double factor = h * r * b / 2.0;
    for (long i = 0; i < n; i++)
        p[i + 1] = (p[i] + (h / 2.0) * f(p[i]) + factor) / (1.0 + factor);
    return p;
}

int main()
{
    long n_steps = 10000000;
    double h = t_max / n_steps;

    cout << "Benchmark C++ - 10^7 pasos" << endl;

    auto [p1, t_euler] = euler(p0, h, n_steps);
    cout << fixed << setprecision(4) << "Euler:    " << t_euler << " s" << endl;

    auto [p2, t_taylor] = taylor2(p0, h, n_steps);
    cout << "Taylor2:  " << t_taylor << " s" << endl;

    auto [p3, t_trapecio] = trapecio(p0, h, n_steps);
    cout << "Trapecio: " << t_trapecio << " s" << endl;

    ofstream f("benchmark_cpp.dat");
    f << setprecision(10);
    f << "Euler " << t_euler << endl;
    f << "Taylor2 " << t_taylor << endl;
    f << "Trapecio " << t_trapecio << endl;
    f.close();

    // Calcular soluciones con h=0.5 para graficar
    long n_vis = 100;
    double h_vis = t_max / n_vis;

    auto p_euler = euler_solucion(p0, h_vis, n_vis);
    auto p_taylor = taylor2_solucion(p0, h_vis, n_vis);
    auto p_trapecio_sol = trapecio_solucion(p0, h_vis, n_vis);

    ofstream sol("soluciones_cpp.dat");
    sol << setprecision(10) << fixed;
    sol << "# t Exacta Euler Taylor2 Trapecio" << endl;
    for (long i = 0; i <= n_vis; i++)
    {
        double t = i * h_vis;
        sol << t << " " << exacta(t) << " " << p_euler[i] << " "
            << p_taylor[i] << " " << p_trapecio_sol[i] << endl;
    }
    sol.close();

    if (p1 < 0 || p2 < 0 || p3 < 0)
        return 1;
    return 0;
}
