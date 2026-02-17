#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os
import glob

# Borrar gráficas anteriores para mantener orden
for archivo in glob.glob('*.png'):
    try:
        os.remove(archivo)
        print(f"Eliminado: {archivo}")
    except:
        pass


# Leer los datos de soluciones
data_py = np.loadtxt('soluciones_python.dat')
data_cpp = np.loadtxt('soluciones_cpp.dat')
data_f90 = np.loadtxt('soluciones_fortran.dat')

t_py = data_py[:, 0]
exacta_py = data_py[:, 1]
euler_py = data_py[:, 2]
taylor_py = data_py[:, 3]
trapecio_py = data_py[:, 4]

t_cpp = data_cpp[:, 0]
exacta_cpp = data_cpp[:, 1]
euler_cpp = data_cpp[:, 2]
taylor_cpp = data_cpp[:, 3]
trapecio_cpp = data_cpp[:, 4]

t_f90 = data_f90[:, 0]
exacta_f90 = data_f90[:, 1]
euler_f90 = data_f90[:, 2]
taylor_f90 = data_f90[:, 3]
trapecio_f90 = data_f90[:, 4]

# Leer datos de benchmark
def leer_benchmark(nombre_archivo):
    tiempos = {}
    with open(nombre_archivo, 'r') as ff:
        for linea in ff:
            partes = linea.strip().split()
            if len(partes) == 2:
                tiempos[partes[0]] = float(partes[1])
    return tiempos

t_py_bench = leer_benchmark('benchmark_python.dat')
t_cpp_bench = leer_benchmark('benchmark_cpp.dat')
t_f90_bench = leer_benchmark('benchmark_fortran.dat')

# =============================================================================
# FIGURA 1: Comparación de tiempos (gráfico de barras)
# =============================================================================
metodos = ['Euler', 'Taylor2', 'Trapecio']
py_times = [t_py_bench[m] for m in metodos]
cpp_times = [t_cpp_bench[m] for m in metodos]
f90_times = [t_f90_bench[m] for m in metodos]

x = np.arange(len(metodos))
ancho = 0.25

fig1, ax = plt.subplots(figsize=(12, 6))
bars1 = ax.bar(x - ancho, py_times, ancho, label='Python', color='#3776ab', alpha=0.8, edgecolor='black')
bars2 = ax.bar(x, cpp_times, ancho, label='C++', color='#00599c', alpha=0.8, edgecolor='black')
bars3 = ax.bar(x + ancho, f90_times, ancho, label='Fortran', color='#734f96', alpha=0.8, edgecolor='black')

# Poner valores arriba de cada barra
for bars in [bars1, bars2, bars3]:
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.4f}s', ha='center', va='bottom', fontsize=9, fontweight='bold')

ax.set_xlabel('Método', fontsize=12, fontweight='bold')
ax.set_ylabel('Tiempo (s)', fontsize=12, fontweight='bold')
ax.set_title('Comparación de Tiempos - Benchmark 10^7 pasos', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(metodos)
ax.legend(fontsize=11)
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('comparacion_tiempos.png', dpi=150)
print("Gráfica guardada: comparacion_tiempos.png")

# =============================================================================
# FIGURA 2: Soluciones y errores (2 filas x 3 columnas)
# =============================================================================
fig2 = plt.figure(figsize=(18, 10))

# --- Fila 1: Soluciones ---
# Python
ax1 = plt.subplot(2, 3, 1)
ax1.plot(t_py, exacta_py, 'k-', linewidth=2, label='Exacta')
ax1.plot(t_py, euler_py, 'r--', label='Euler')
ax1.plot(t_py, taylor_py, 'g--', label='Taylor2')
ax1.plot(t_py, trapecio_py, 'b--', label='Trapecio')
ax1.set_xlabel('t')
ax1.set_ylabel('p(t)')
ax1.set_title(f'Python (Euler: {t_py_bench["Euler"]:.4f}s, Taylor: {t_py_bench["Taylor2"]:.4f}s, Trap: {t_py_bench["Trapecio"]:.4f}s)')
ax1.legend()
ax1.grid(True, alpha=0.3)

# C++
ax2 = plt.subplot(2, 3, 2)
ax2.plot(t_cpp, exacta_cpp, 'k-', linewidth=2, label='Exacta')
ax2.plot(t_cpp, euler_cpp, 'r--', label='Euler')
ax2.plot(t_cpp, taylor_cpp, 'g--', label='Taylor2')
ax2.plot(t_cpp, trapecio_cpp, 'b--', label='Trapecio')
ax2.set_xlabel('t')
ax2.set_ylabel('p(t)')
ax2.set_title(f'C++ (Euler: {t_cpp_bench["Euler"]:.4f}s, Taylor: {t_cpp_bench["Taylor2"]:.4f}s, Trap: {t_cpp_bench["Trapecio"]:.4f}s)')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Fortran
ax3 = plt.subplot(2, 3, 3)
ax3.plot(t_f90, exacta_f90, 'k-', linewidth=2, label='Exacta')
ax3.plot(t_f90, euler_f90, 'r--', label='Euler')
ax3.plot(t_f90, taylor_f90, 'g--', label='Taylor2')
ax3.plot(t_f90, trapecio_f90, 'b--', label='Trapecio')
ax3.set_xlabel('t')
ax3.set_ylabel('p(t)')
ax3.set_title(f'Fortran (Euler: {t_f90_bench["Euler"]:.4f}s, Taylor: {t_f90_bench["Taylor2"]:.4f}s, Trap: {t_f90_bench["Trapecio"]:.4f}s)')
ax3.legend()
ax3.grid(True, alpha=0.3)

# --- Fila 2: Errores ---
# Python
ax4 = plt.subplot(2, 3, 4)
ax4.semilogy(t_py, np.abs(euler_py - exacta_py), 'r-', label='Euler')
ax4.semilogy(t_py, np.abs(taylor_py - exacta_py), 'g-', label='Taylor2')
ax4.semilogy(t_py, np.abs(trapecio_py - exacta_py), 'b-', label='Trapecio')
ax4.set_xlabel('t')
ax4.set_ylabel('|Error|')
ax4.set_title('Error Python')
ax4.legend()
ax4.grid(True, alpha=0.3)

# C++
ax5 = plt.subplot(2, 3, 5)
ax5.semilogy(t_cpp, np.abs(euler_cpp - exacta_cpp), 'r-', label='Euler')
ax5.semilogy(t_cpp, np.abs(taylor_cpp - exacta_cpp), 'g-', label='Taylor2')
ax5.semilogy(t_cpp, np.abs(trapecio_cpp - exacta_cpp), 'b-', label='Trapecio')
ax5.set_xlabel('t')
ax5.set_ylabel('|Error|')
ax5.set_title('Error C++')
ax5.legend()
ax5.grid(True, alpha=0.3)

# Fortran
ax6 = plt.subplot(2, 3, 6)
ax6.semilogy(t_f90, np.abs(euler_f90 - exacta_f90), 'r-', label='Euler')
ax6.semilogy(t_f90, np.abs(taylor_f90 - exacta_f90), 'g-', label='Taylor2')
ax6.semilogy(t_f90, np.abs(trapecio_f90 - exacta_f90), 'b-', label='Trapecio')
ax6.set_xlabel('t')
ax6.set_ylabel('|Error|')
ax6.set_title('Error Fortran')
ax6.legend()
ax6.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('comparacion_soluciones_errores.png', dpi=150)
print("Gráfica guardada: comparacion_soluciones_errores.png")

print("\n=== Resumen de Tiempos ===")
print(f"Python   - Euler: {py_times[0]:.4f}s, Taylor2: {py_times[1]:.4f}s, Trapecio: {py_times[2]:.4f}s")
print(f"C++      - Euler: {cpp_times[0]:.4f}s, Taylor2: {cpp_times[1]:.4f}s, Trapecio: {cpp_times[2]:.4f}s")
print(f"Fortran  - Euler: {f90_times[0]:.4f}s, Taylor2: {f90_times[1]:.4f}s, Trapecio: {f90_times[2]:.4f}s")

# =============================
# NUEVA FIGURA: Comparación por método
# =============================
fig3, axes = plt.subplots(2, 3, figsize=(18, 10))
metodos = ['Euler', 'Taylor2', 'Trapecio']
colores = ['#3776ab', '#00599c', '#734f96']
labels = ['Python', 'C++', 'Fortran']
# Soluciones
for j, (metodo, idx) in enumerate(zip(metodos, [2,3,4])):
    ax = axes[0, j]
    ax.plot(t_py, exacta_py, 'k-', linewidth=2, label='Exacta')
    ax.plot(t_py, data_py[:, idx], color=colores[0], linestyle='--', label='Python')
    ax.plot(t_cpp, data_cpp[:, idx], color=colores[1], linestyle='--', label='C++')
    ax.plot(t_f90, data_f90[:, idx], color=colores[2], linestyle='--', label='Fortran')
    ax.set_xlabel('t')
    ax.set_ylabel('p(t)')
    ax.set_title(f'Solución {metodo}')
    ax.legend()
    ax.grid(True, alpha=0.3)
# Errores
for j, (metodo, idx) in enumerate(zip(metodos, [2,3,4])):
    ax = axes[1, j]
    ax.semilogy(t_py, np.abs(data_py[:, idx] - exacta_py), color=colores[0], label='Python')
    ax.semilogy(t_cpp, np.abs(data_cpp[:, idx] - exacta_cpp), color=colores[1], label='C++')
    ax.semilogy(t_f90, np.abs(data_f90[:, idx] - exacta_f90), color=colores[2], label='Fortran')
    ax.set_xlabel('t')
    ax.set_ylabel('|Error|')
    ax.set_title(f'Error {metodo}')
    ax.legend()
    ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('comparacion_metodo_vs_lenguaje.png', dpi=150)
print("Gráfica guardada: comparacion_metodo_vs_lenguaje.png")

