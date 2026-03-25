#ifndef STURM_LIB_HPP
#define STURM_LIB_HPP

#include <vector>
#include <cmath>
#include <iostream>

namespace SturmEngine {
    using namespace std;

    // Cuenta cuántos autovalores son menores que e_test
    int contar_sturm(int n, const vector<double>& diag, const vector<double>& subdiag, double e_test) {
        int count = 0;
        double q = diag[0] - e_test;
        
        if (q < 0.0) count++;

        for (int i = 1; i < n; i++) {
            if (abs(q) < 1.0e-30) q = 1.0e-30; // Evitar división por cero
            
            // Recurrencia de Sturm: q_i = (d_i - E) - b_{i-1}^2 / q_{i-1}
            q = (diag[i] - e_test) - (subdiag[i-1] * subdiag[i-1]) / q;
            
            if (q < 0.0) count++;
        }
        return count;
    }

    // Bisección para encontrar exactamente el nivel de energía deseado
    double buscar_energia(int n, const vector<double>& diag, const vector<double>& subdiag, 
                          int n_estado, double e_min, double e_max, double tol) {
        double low = e_min;
        double high = e_max;
        
        while ((high - low) > tol) {
            double mid = 0.5 * (low + high);
            if (contar_sturm(n, diag, subdiag, mid) >= n_estado) {
                high = mid;
            } else {
                low = mid;
            }
        }
        return 0.5 * (low + high);
    }

    // Algoritmo de Thomas para resolver sistemas tridiagonales (Inverse Iteration)
    void solver_tridiag(const vector<double>& sub, const vector<double>& diag, 
                        const vector<double>& sup, const vector<double>& r, vector<double>& u) {
        int n = r.size();
        vector<double> gam(n);
        double bet = diag[0];
        u[0] = r[0] / bet;
        
        for (int j = 1; j < n; j++) {
            gam[j] = sup[j-1] / bet;
            bet = diag[j] - sub[j-1] * gam[j];
            u[j] = (r[j] - sub[j-1] * u[j-1]) / bet;
        }
        for (int j = n - 2; j >= 0; j--) {
            u[j] -= gam[j+1] * u[j+1];
        }
    }

    // Normalización estándar
    void normalizar(vector<double>& v) {
        double norma = 0.0;
        for (double val : v) norma += val * val;
        norma = sqrt(norma);
        for (double& val : v) val /= norma;
    }
}

#endif