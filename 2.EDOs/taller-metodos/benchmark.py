import numpy as np
import time

# Parámetros
b, d, r, p0, t_max = 0.02, 0.015, 0.1, 0.01, 50.0

def f(p): return r * b * (1 - p)
def f_prime(p): return -r * b * f(p)
def exacta(t): return 1 - (1 - p0) * np.exp(-r * b * t)

def euler(p0, h, n):
    p = np.zeros(n + 1)
    p[0] = p0
    start = time.time()
    for i in range(n):
        p[i + 1] = p[i] + h * f(p[i])
    return p, time.time() - start

def taylor2(p0, h, n):
    p = np.zeros(n + 1)
    p[0] = p0
    start = time.time()
    for i in range(n):
        p[i + 1] = p[i] + h * f(p[i]) + (h**2 / 2) * f_prime(p[i])
    return p, time.time() - start

def trapecio(p0, h, n):
    p = np.zeros(n + 1)
    p[0] = p0
    factor = h * r * b / 2
    start = time.time()
    for i in range(n):
        p[i + 1] = (p[i] + (h / 2) * f(p[i]) + factor) / (1 + factor)
    return p, time.time() - start

if __name__ == '__main__':
    # Benchmark con 10^7 pasos
    n_steps_bench = int(1e7)
    h_bench = t_max / n_steps_bench
    
    print(f"Benchmark Python - 10^7 pasos")
    
    _, t_euler = euler(p0, h_bench, n_steps_bench)
    print(f"Euler:    {t_euler:.4f} s")
    
    _, t_taylor = taylor2(p0, h_bench, n_steps_bench)
    print(f"Taylor2:  {t_taylor:.4f} s")
    
    _, t_trapecio = trapecio(p0, h_bench, n_steps_bench)
    print(f"Trapecio: {t_trapecio:.4f} s")
    
    # Guardar benchmark
    with open('benchmark_python.dat', 'w') as archivo:
        archivo.write(f"Euler {t_euler:.10f}\n")
        archivo.write(f"Taylor2 {t_taylor:.10f}\n")
        archivo.write(f"Trapecio {t_trapecio:.10f}\n")
    
    # Calcular soluciones con h=0.5 para graficar
    n_steps_vis = 100
    h_vis = t_max / n_steps_vis
    
    p_euler, _ = euler(p0, h_vis, n_steps_vis)
    p_taylor, _ = taylor2(p0, h_vis, n_steps_vis)
    p_trapecio, _ = trapecio(p0, h_vis, n_steps_vis)
    
    t = np.linspace(0, t_max, n_steps_vis + 1)
    p_exact = exacta(t)
    
    # Guardar soluciones
    with open('soluciones_python.dat', 'w') as archivo:
        archivo.write("# t Exacta Euler Taylor2 Trapecio\n")
        for i in range(len(t)):
            archivo.write(f"{t[i]:.6f} {p_exact[i]:.10f} {p_euler[i]:.10f} {p_taylor[i]:.10f} {p_trapecio[i]:.10f}\n")
