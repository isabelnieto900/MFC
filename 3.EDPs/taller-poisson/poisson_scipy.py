"""
poisson_scipy.py — Solución de la ecuación de Poisson con scipy.sparse
  ∇²V = (x²+y²)·exp(xy)   en Ω = [0,2]×[0,1]
  Solución analítica: V(x,y) = exp(xy)

Uso:
    python poisson_scipy.py N
    python poisson_scipy.py N --all     # corre 32 64 128 256
"""

import sys, os, time, csv, math, tracemalloc
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

# ── Condiciones de frontera ──────────────────────────────────────────────────
bc_left  = lambda y: 1.0
bc_right = lambda y: math.exp(2.0 * y)
bc_bot   = lambda x: 1.0
bc_top   = lambda x: math.exp(x)
V_exact  = lambda x, y: math.exp(x * y)
f_rhs    = lambda x, y: (x**2 + y**2) * math.exp(x * y)

# ── Solver ───────────────────────────────────────────────────────────────────
def solve(N):
    M  = N
    hx = 2.0 / (N + 1)
    hy = 1.0 / (N + 1)
    ax = 1.0 / hx**2
    ay = 1.0 / hy**2
    ac = -2.0 * ax - 2.0 * ay
    sz = N * M

    t0 = time.perf_counter()
    tracemalloc.start()

    A = lil_matrix((sz, sz), dtype=np.float64)
    b = np.zeros(sz)

    for j in range(M):
        y = (j + 1) * hy
        for i in range(N):
            x  = (i + 1) * hx
            k  = j * N + i

            A[k, k] = ac

            if i > 0:       A[k, k-1]   = ax
            else:           b[k] -= ax * bc_left(y)

            if i < N-1:     A[k, k+1]   = ax
            else:           b[k] -= ax * bc_right(y)

            if j > 0:       A[k, k-N]   = ay
            else:           b[k] -= ay * bc_bot(x)

            if j < M-1:     A[k, k+N]   = ay
            else:           b[k] -= ay * bc_top(x)

            b[k] += f_rhs(x, y)

    A_csr = csr_matrix(A)

    u = spsolve(A_csr, b)

    elapsed_ms = (time.perf_counter() - t0) * 1e3
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    mem_MB = peak / 1024**2

    # Errores
    errL2 = 0.0; errMax = 0.0
    for j in range(M):
        y = (j + 1) * hy
        for i in range(N):
            x   = (i + 1) * hx
            diff = abs(u[j*N+i] - V_exact(x, y))
            errL2  += diff**2
            errMax  = max(errMax, diff)
    errL2 = math.sqrt(errL2 / (N * M))

    return elapsed_ms, mem_MB, errL2, errMax, u, hx, hy

# ── CSV de solución ───────────────────────────────────────────────────────────
def save_csv(fname, u, N):
    hx = 2.0 / (N + 1); hy = 1.0 / (N + 1)
    with open(fname, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["x", "y", "V_num", "V_exact"])
        for j in range(N):
            y = (j+1)*hy
            for i in range(N):
                x = (i+1)*hx
                w.writerow([x, y, u[j*N+i], V_exact(x, y)])

# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    args = sys.argv[1:]
    if not args:
        print("Uso: python poisson_scipy.py N [--all]"); sys.exit(1)

    run_all = "--all" in args
    sizes   = [32, 64, 128, 256] if run_all else [int(args[0])]

    metrics_file = "metrics_python.csv"
    write_header = not os.path.exists(metrics_file)
    mf = open(metrics_file, "a", newline="")
    mw = csv.writer(mf)
    if write_header:
        mw.writerow(["N", "unknowns", "time_ms", "mem_MB", "errL2", "errMax"])

    for N in sizes:
        ms, mem, eL2, eMax, u, hx, hy = solve(N)
        print(f"Python/scipy  malla {N}x{N}  ({N*N} incógnitas)")
        print(f"  Tiempo    : {ms:.2f} ms")
        print(f"  Mem. RSS  : {mem:.4f} MB")
        print(f"  Error L2  : {eL2:.6e}")
        print(f"  Error Max : {eMax:.6e}\n")
        mw.writerow([N, N*N, f"{ms:.4f}", f"{mem:.6f}", f"{eL2:.8e}", f"{eMax:.8e}"])
        csv_name = f"sol_{N}_python.csv"
        save_csv(csv_name, u, N)
        print(f"  CSV: {csv_name}\n")

    mf.close()

if __name__ == "__main__":
    main()
