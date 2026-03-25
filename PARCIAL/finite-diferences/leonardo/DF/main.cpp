#include "optimization.hpp" // Asumiendo que guardaste lo anterior en este header
#include <iomanip>
#include <chrono>
double z_global;

// Función de energía (func)
double energy_function(const vector<double>& x) {
    double a = x[0];
    double b = x[1];
    double y = 2.0 * a * b / pow(2.0 * a + b, 5);
    
    return pow(a, 2) - 2.0 * z_global * a + (5.0 / 8.0) * a + (1.0 / 8.0) * pow(b, 2) - (z_global * b / 4.0) +
           y * (8.0 * pow(a, 4) + 20.0 * pow(a, 3) * b + 12.0 * pow(a, 2) * pow(b, 2) + 10.0 * a * pow(b, 3) + pow(b, 4));
}

// Energía de Ionización (ip_func)
double ip_function(double a, double b) {
    double a2 = a * a;
    double b2 = b * b;
    double y = 2.0 * a * b / pow(2.0 * a + b, 5);
    
    return (1.0 / 8.0) * b2 - (z_global * b / 4.0) +
           y * (8.0 * a2 * a2 + 20.0 * a2 * a * b + 12.0 * a2 * b2 + 10.0 * a * b * b2 + b2 * b2);
}

int main() {
    
    int n = 2;
    double ftol = 1.0e-8;
    int iter;
    double fret;
    
    vector<double> p = {1.5, 1.5};
    vector<vector<double>> xi(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) xi[i][i] = 1.0;

    cout << "CUAL ES EL VALOR DE Z: ";
    cin >> z_global;
    auto start = std::chrono::high_resolution_clock::now();
    Optimization::powell(p, xi, ftol, iter, fret, energy_function);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    cout << "Tiempo de ejecución: " << duration.count() << " segundos" << endl;

    cout << fixed << setprecision(8);
    cout << "--------------------------------------" << endl;
    cout << "ITERACIONES:        = " << iter << endl;
    cout << "ALFA, BETA:         = " << p[0] << ", " << p[1] << endl;
    cout << "ENER. ESTADO BASE (eV) = " << 27.211 * fret << endl;
    cout << "ENER. IONIZACION  (eV) = " << 27.211 * ip_function(p[0], p[1]) << endl;
    cout << "--------------------------------------" << endl;

    return 0;
    }