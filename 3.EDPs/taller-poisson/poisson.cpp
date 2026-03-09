/*
 * poisson.cpp — Solución numérica de la ecuación de Poisson
 *
 *   ∇²V = (x² + y²)·exp(x·y)   en Ω = [0,2] × [0,1]
 *
 *   Condiciones de Dirichlet:
 *     V(0,y) = 1          V(2,y) = exp(2y)
 *     V(x,0) = 1          V(x,1) = exp(x)
 *
 *   Solución analítica: V(x,y) = exp(x·y)
 *
 * Compilación:
 *   g++ -O2 -std=c++17 -I/usr/include/eigen3 poisson.cpp -o poisson
 *
 * Uso:
 *   ./poisson N
 *     N : número de puntos interiores en cada dirección (e.g. 32)
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <chrono>
#include <string>
#include <vector>
#include <iomanip>

// ── Eigen (sparse) ──────────────────────────────────────────────────────────
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
// ── Memoria real del proceso (VmRSS de /proc/self/status) ───────────────────
static double get_rss_mb()
{
    std::ifstream f("/proc/self/status");
    std::string line;
    while (std::getline(f, line))
    {
        if (line.rfind("VmRSS:", 0) == 0)
        {
            std::istringstream ss(line.substr(6));
            long kb;
            ss >> kb;
            return kb / 1024.0;
        }
    }
    return 0.0;
}
// ════════════════════════════════════════════════════════════════════════════
//  Funciones del problema
// ════════════════════════════════════════════════════════════════════════════

inline double f_rhs(double x, double y)
{
    return (x * x + y * y) * std::exp(x * y);
}

inline double V_exact(double x, double y)
{
    return std::exp(x * y);
}

// Condiciones de frontera
inline double bc_left(double y) { return 1.0; }
inline double bc_right(double y) { return std::exp(2.0 * y); }
inline double bc_bot(double x) { return 1.0; }
inline double bc_top(double x) { return std::exp(x); }

// ════════════════════════════════════════════════════════════════════════════
//  Construcción del sistema (compartida: mismos A e índices para ambos solvers)
// ════════════════════════════════════════════════════════════════════════════

struct GridData
{
    int N, M;
    double hx, hy;
};

inline int idx(int i, int j, int N) { return j * N + i; }

// ════════════════════════════════════════════════════════════════════════════
//  Eigen: construye y resuelve
// ════════════════════════════════════════════════════════════════════════════

std::vector<double> solve_eigen(const GridData &g,
                                double &elapsed_ms, double &mem_MB)
{
    const int N = g.N, M = g.M, sz = N * M;
    const double hx = g.hx, hy = g.hy;
    const double ax = 1.0 / (hx * hx), ay = 1.0 / (hy * hy);
    const double ac = -2.0 * ax - 2.0 * ay;

    using SpMat = Eigen::SparseMatrix<double>;
    using T = Eigen::Triplet<double>;

    auto t0 = std::chrono::high_resolution_clock::now();
    double mem_before = get_rss_mb(); // RSS antes de construir la matriz

    std::vector<T> triplets;
    triplets.reserve(5 * sz);
    Eigen::VectorXd b(sz);
    b.setZero();

    for (int j = 0; j < M; ++j)
    {
        double y = (j + 1) * hy;
        for (int i = 0; i < N; ++i)
        {
            double x = (i + 1) * hx;
            int k = idx(i, j, N);

            // diagonal
            triplets.emplace_back(k, k, ac);

            // vecino izquierdo
            if (i > 0)
                triplets.emplace_back(k, idx(i - 1, j, N), ax);
            else
                b[k] -= ax * bc_left(y);

            // vecino derecho
            if (i < N - 1)
                triplets.emplace_back(k, idx(i + 1, j, N), ax);
            else
                b[k] -= ax * bc_right(y);

            // vecino inferior
            if (j > 0)
                triplets.emplace_back(k, idx(i, j - 1, N), ay);
            else
                b[k] -= ay * bc_bot(x);

            // vecino superior
            if (j < M - 1)
                triplets.emplace_back(k, idx(i, j + 1, N), ay);
            else
                b[k] -= ay * bc_top(x);

            // término fuente
            b[k] += f_rhs(x, y);
        }
    }

    SpMat A(sz, sz);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();

    Eigen::SparseLU<SpMat> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    if (solver.info() != Eigen::Success)
    {
        std::cerr << "[Eigen] Error en factorización.\n";
        return {};
    }
    Eigen::VectorXd u = solver.solve(b);

    auto t1 = std::chrono::high_resolution_clock::now();
    elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    mem_MB = get_rss_mb() - mem_before; // delta RSS atribuible al solver

    return std::vector<double>(u.data(), u.data() + sz);
}

// ════════════════════════════════════════════════════════════════════════════
//  Métricas de error
// ════════════════════════════════════════════════════════════════════════════

void compute_errors(const std::vector<double> &u, const GridData &g,
                    double &errL2, double &errMax)
{
    const int N = g.N, M = g.M;
    errL2 = 0.0;
    errMax = 0.0;
    for (int j = 0; j < M; ++j)
    {
        double y = (j + 1) * g.hy;
        for (int i = 0; i < N; ++i)
        {
            double x = (i + 1) * g.hx;
            double diff = std::abs(u[idx(i, j, N)] - V_exact(x, y));
            errL2 += diff * diff;
            errMax = std::max(errMax, diff);
        }
    }
    errL2 = std::sqrt(errL2 / (N * M));
}

// ════════════════════════════════════════════════════════════════════════════
//  Guardar solución en CSV para visualización con Python
// ════════════════════════════════════════════════════════════════════════════

void save_csv(const std::string &fname,
              const std::vector<double> &u,
              const GridData &g)
{
    std::ofstream f(fname);
    f << std::scientific << std::setprecision(10);
    f << "x,y,V_num,V_exact\n";
    for (int j = 0; j < g.M; ++j)
    {
        double y = (j + 1) * g.hy;
        for (int i = 0; i < g.N; ++i)
        {
            double x = (i + 1) * g.hx;
            f << x << "," << y << ","
              << u[idx(i, j, g.N)] << ","
              << V_exact(x, y) << "\n";
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
//  Main
// ════════════════════════════════════════════════════════════════════════════

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cerr << "Uso: " << argv[0] << " N\n"
                  << "  N : puntos interiores (e.g. 32)\n";
        return 1;
    }

    int N = std::atoi(argv[1]);
    if (N < 2)
    {
        std::cerr << "N debe ser >= 2\n";
        return 1;
    }

    GridData g;
    g.N = N;
    g.M = N;
    g.hx = 2.0 / (N + 1);
    g.hy = 1.0 / (N + 1);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Poisson FD (Eigen)  malla " << N << "x" << N
              << "  (" << (long long)N * N << " incognitas)\n";

    double ms = 0, mem = 0;
    auto u = solve_eigen(g, ms, mem);
    if (u.empty())
        return 1;

    double eL2, eMax;
    compute_errors(u, g, eL2, eMax);

    std::cout << "  Tiempo    : " << ms << " ms\n"
              << "  Mem. RSS  : " << mem << " MB\n"
              << "  Error L2  : " << eL2 << "\n"
              << "  Error Max : " << eMax << "\n";

    // CSV de la solución
    std::string csvname = "sol_" + std::to_string(N) + "_eigen.csv";
    save_csv(csvname, u, g);
    std::cout << "  CSV: " << csvname << "\n";

    // Métricas acumulativas
    std::ofstream mf("metrics.csv", std::ios::app);
    if (mf.is_open())
    {
        mf.seekp(0, std::ios::end);
        if (mf.tellp() == 0)
            mf << "N,unknowns,time_ms,mem_MB,errL2,errMax\n";
        mf << N << "," << (long long)N * N << ","
           << ms << "," << mem << "," << eL2 << "," << eMax << "\n";
    }

    return 0;
}
