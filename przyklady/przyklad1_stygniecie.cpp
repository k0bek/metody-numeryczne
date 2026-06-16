/**
 * @file przyklad1_stygniecie.cpp
 * @brief Przyklad 1: Stygniecie ciala (RRZ)
 *
 * Model: T'(t) = -alpha * T^4,  T(0) = T0
 * Rozwiazanie numeryczne: Euler, Heun, Punkt srodkowy, RK4.
 * Porownanie z rozwiazaniem analitycznym.
 *
 * Kompilacja: g++ -std=c++17 -I../include przyklad1_stygniecie.cpp ../src/*.cpp -o przyklad1
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "../include/rownania_rozniczkowe.h"

using namespace std;
using namespace biblioteka_mn;

// Dokladne rozwiazanie analityczne (Wolfram)
double dokladne(double t) {
    return 76010000.0 / pow(1317447908403.0 * t + 1e12, 1.0/3.0);
}

int main() {
    const double T0      = 7601.0;
    const double alpha   = 1e-12;
    const double t_koniec = 7601.0;
    const int    N        = 2000;

    auto f = [alpha](double /*t*/, double T) {
        return -alpha * T * T * T * T;
    };

    auto eu = euler(f, 0.0, T0, t_koniec, N);
    auto he = heun (f, 0.0, T0, t_koniec, N);
    auto mi = punkt_srodkowy(f, 0.0, T0, t_koniec, N);
    auto rk = rk4  (f, 0.0, T0, t_koniec, N);

    cout << fixed << setprecision(4);
    cout << "=== Stygniecie ciala - porownanie metod ===\n\n";
    cout << setw(10) << "t [s]"
         << setw(14) << "Dokladny [K]"
         << setw(14) << "Euler [K]"
         << setw(12) << "Blad%"
         << setw(14) << "RK4 [K]"
         << setw(12) << "Blad%"
         << "\n";
    cout << string(80, '-') << "\n";

    int krok = N / 10;
    for (int i = 0; i <= N; i += krok) {
        double t   = eu[i].first;
        double dok = dokladne(t);
        double blad_eu = fabs(eu[i].second - dok) / dok * 100.0;
        double blad_rk = fabs(rk[i].second - dok) / dok * 100.0;
        cout << setw(10) << t
             << setw(14) << dok
             << setw(14) << eu[i].second
             << setw(11) << blad_eu << "%"
             << setw(14) << rk[i].second
             << setw(11) << blad_rk << "%"
             << "\n";
    }

    cout << "\nBlad bezwzgledny w t_koniec:\n";
    cout << "  Euler:           " << fabs(eu.back().second - dokladne(t_koniec)) << " K\n";
    cout << "  Heun:            " << fabs(he.back().second - dokladne(t_koniec)) << " K\n";
    cout << "  Punkt srodkowy:  " << fabs(mi.back().second - dokladne(t_koniec)) << " K\n";
    cout << "  RK4:             " << fabs(rk.back().second - dokladne(t_koniec)) << " K\n";
    return 0;
}
