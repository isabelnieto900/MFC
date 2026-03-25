#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "sturm_lib.hpp"
#include <chrono>
using namespace std;

// Función del potencial del Litio (usando el alfa de Fortran para consistencia)
double pot_litio(double r, double l) {
    const double alfa = 2.535930;
    double v_c = -(1.0 / r) - (2.0 / r) * (1.0 + alfa * r) * exp(-2.0 * alfa * r);
    double v_l = l * (l + 1.0) / (2.0 * r * r);
    return v_c + v_l;
}

int main() {
    const int npt = 6000;
    const double d_max = 25.0;
    const double h = d_max / npt;
    const double tol = 1e-12;
    double l_angular;

    cout << "MOMENTUM ANGULAR (l) = ";
    cin >> l_angular;
    auto start = std::chrono::high_resolution_clock::now();
    // 1. Construcción de la Matriz Hamiltoniana (Tridiagonal)
    vector<double> diag(npt);
    vector<double> subdiag(npt - 1);
    
    for (int i = 0; i < npt; i++) {
        double r = (i + 1) * h;
        // H_ii = 1/h^2 + V(r)
        diag[i] = (1.0 / (h * h)) + pot_litio(r, l_angular);
        if (i < npt - 1) {
            // H_i,i+1 = -1/(2h^2)
            subdiag[i] = -1.0 / (2.0 * h * h);
        }
    }

    ofstream archivo("LITIO_STURM_CPP.DAT");
    cout << fixed << setprecision(8);
    cout << "\nCalculando estados con Secuencia de Sturm...\n";
    cout << "----------------------------------------------\n";

    // 2. Encontrar los primeros 4 estados
    for (int n = 1; n <= 4; n++) {
        // Buscamos la energía mediante bisección controlada por Sturm
        double e_au = SturmEngine::buscar_energia(npt, diag, subdiag, n, -50.0, 0.0, tol);
        
        cout << "N=" << n << " | E(au)=" << setw(12) << e_au 
             << " | E(eV)=" << setw(12) << e_au * 27.211 << endl;

        // 3. Iteración Inversa para la Función de Onda
        vector<double> vk(npt, 1.0);
        vector<double> v_next(npt);
        vector<double> diag_shift(npt);
        vector<double> sup(npt - 1, -1.0 / (2.0 * h * h));
        vector<double> sub(npt - 1, -1.0 / (2.0 * h * h));

        for (int i = 0; i < npt; i++) diag_shift[i] = diag[i] - e_au;

        for (int iter = 0; iter < 15; iter++) {
            SturmEngine::solver_tridiag(sub, diag_shift, sup, vk, v_next);
            SturmEngine::normalizar(v_next);
            vk = v_next;
        }
        // 1. Escala Física: Dividir por la raíz de h para integral(psi^2 dr) = 1
        double factor_h = sqrt(h);
        for (double &val : vk) {
            val /= factor_h;
        }

        // 2. Corrección de Fase: Forzar que el inicio sea positivo
        // Buscamos el primer valor que supere el umbral de ruido numérico
        for (int k = 0; k < npt; k++) {
            if (abs(vk[k]) > 1e-5) {
                if (vk[k] < 0.0) {
                    for (double &val : vk) val = -val; // Invertimos toda la función
                }
                break; // Ya corregimos la fase, salimos del buscador
            }
        }
        // Guardar en archivo (1 de cada 10 puntos)
        for (int i = 0; i < npt; i += 10) {
            archivo << (i + 1) * h << " " << vk[i] << "\n";
        }
        archivo << "\n";
    }

    archivo.close();
    cout << "----------------------------------------------\n";
    cout << "Proceso completado. Datos en LITIO_STURM_CPP.DAT\n";
    auto end = std::chrono::high_resolution_clock::now();   // Fin
    std::chrono::duration<double> duration = end - start;

    cout << "Tiempo de ejecución: " << duration.count() << " segundos" << endl;
    return 0;
}