#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <chrono>
#include "herramientas.hpp"

using namespace std;

// Variables globales para simular el COMMON BLOCK /params/
double h_global, l_global;
int npt_global;

// Definición del potencial local para el Litio
double pot_local(double r, double l) {
    const double alfa = 2.535930; // CORRECCIÓN: Unificado con Fortran
    double v_central = -(1.0 / r) - (2.0 / r) * (1.0 + alfa * r) * exp(-2.0 * alfa * r);
    double v_centrifugo = l * (l + 1.0) / (2.0 * r * r);
    return v_central + v_centrifugo;
}

// Función determinante basada en la integración numérica
double determ(double e) {
    double h2 = h_global * h_global;
    double p0 = 0.0; // psi(0)
    double p1 = 1.0; // psi(h) - Punto de partida
    double p2;

    for (int i = 1; i <= npt_global; i++) {
        double r = i * h_global;
        double pot = pot_local(r, l_global);
        
        if (i == 1) {
            p1 = (2.0 + 2.0 * h2 * pot - 2.0 * e * h2);
        } else {
            p2 = (2.0 + 2.0 * h2 * pot - 2.0 * h2 * e) * p1 - p0;
            p0 = p1;
            p1 = p2;
        }

        if (std::abs(p1) > 1.0e20) {
            p1 *= 1.0e-20;
            p0 *= 1.0e-20;
        }
    }
    return p1;
}

int main() {
    
    npt_global = 6000;
    double r_max = 25.0;
    h_global = r_max / npt_global;

    cout << "MOMENTUM ANGULAR (l) = ";
    cin >> l_global;
    
    vector<double> xb1, xb2;
    cout << "Buscando intervalos de energia..." << endl;
    auto start = std::chrono::high_resolution_clock::now();
    Herramientas::zbrak(determ, -10.0, 0.0, 20000, xb1, xb2);

    ofstream salida("LITIO_CPP.DAT");
    salida << fixed << setprecision(8);

    cout << "\nNiveles de energia encontrados:" << endl;
    cout << "--------------------------------------" << endl;

    for (size_t i = 0; i < xb1.size(); i++) {
        double raiz = Herramientas::zbrent(determ, xb1[i], xb2[i], 1.0e-12);
        
        if (raiz < -1.0e-8) {
            cout << "N=" << i + 1 << " | E(au)=" << setw(12) << raiz 
                 << " | E(eV)=" << setw(12) << raiz * 27.211 << endl;

            vector<double> diag(npt_global);
            vector<double> sub(npt_global, -1.0);
            vector<double> sup(npt_global, -1.0);
            vector<double> vk(npt_global, 0.0);
            vector<double> v_next(npt_global);

            for (int k = 0; k < npt_global; k++) {
                double r = (k + 1) * h_global;
                diag[k] = 2.0 + 2.0 * h_global * h_global * (pot_local(r, l_global) - raiz);
            }

            vk[0] = 1.0;
            for (int iter = 0; iter < 15; iter++) {
                Herramientas::tridag(sub, diag, sup, vk, v_next);
                Herramientas::norma(v_next);
                vk = v_next;
            }

            // --- CORRECCIÓN DE FASE (Reflexión) ---
            if (vk[0] < 0.0) {
                for (int k = 0; k < npt_global; k++) {
                    vk[k] = -vk[k];
                }
            }
            // --------------------------------------

            for (int k = 0; k < npt_global; k += 10) {
                salida << (k + 1) * h_global << " " << vk[k] << endl;
            }
            salida << endl; 
        }
    }
    
    salida.close();
    cout << "--------------------------------------" << endl;
    cout << "Datos de funciones de onda guardados en LITIO_CPP.DAT" << endl;
    auto end = std::chrono::high_resolution_clock::now();   // Fin
	  std::chrono::duration<double> duration = end - start;

	  cout << "Tiempo de ejecución: " << duration.count() << " segundos" << endl;	


    return 0;
	
}