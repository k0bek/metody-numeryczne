/**
 * @file przyklad2_pierwiastki_calki.cpp
 * @brief Przyklad 2: Miejsca zerowe i calkowanie
 *
 * Demonstruje:
 *  - Wszystkie cztery metody wyznaczania miejsc zerowych
 *  - Porownanie dokladnosci metod calkowania (trapez, Simpson, G-L)
 *
 * Kompilacja: g++ -std=c++17 -I../include przyklad2_pierwiastki_calki.cpp ../src/*.cpp -o przyklad2
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "../include/rownania_nieliniowe.h"
#include "../include/calkowanie.h"

using namespace std;
using namespace biblioteka_mn;

int main() {
    cout << fixed << setprecision(10);

    // -- Miejsca zerowe -----------------------------------------------------
    cout << "=== MIEJSCA ZEROWE ===\n\n";

    auto f1  = [](double x) { return cosh(x) - sqrt(x) - 1.0; };
    auto df1 = [](double x) { return sinh(x) - 0.5/sqrt(x); };

    cout << "f(x) = cosh(x) - sqrt(x) - 1,  pierwiastek na [0.01, 4]\n";
    cout << "  Bisekcja:          " << bisekcja(f1, 0.01, 4.0, 1e-8) << "\n";
    cout << "  Regula falsi:      " << regula_falsi(f1, 0.01, 4.0, 1e-8) << "\n";
    cout << "  Newton analitycz.: " << newton(f1, df1, 1.0, 1e-10) << "\n";
    cout << "  Newton numerycz.:  " << newton_numeryczny(f1, 1.0, 1e-10) << "\n";
    cout << "  Sieczne:           " << sieczne(f1, 0.5, 1.5, 1e-10) << "\n\n";

    auto f2  = [](double x) { return pow(x,3) + x*x - 3*x - 3; };
    auto df2 = [](double x) { return 3*x*x + 2*x - 3; };

    cout << "f(x) = x^3 + x^2 - 3x - 3, pierwiastek ~1.732\n";
    cout << "  Newton analitycz.: " << newton(f2, df2, 2.0, 1e-10) << "\n";
    cout << "  Sieczne:           " << sieczne(f2, 2.0, 3.0, 1e-10) << "\n\n";

    cout << "f(x) = x^2 + 1 na [-2, 2]  (brak pierwiastka):\n";
    try {
        bisekcja([](double x){ return x*x+1; }, -2.0, 2.0);
    } catch (const exception& e) {
        cout << "  Blad (oczekiwany): " << e.what() << "\n";
    }

    // -- Calkowanie -----------------------------------------------------------
    cout << "\n=== CALKOWANIE NUMERYCZNE ===\n\n";

    auto g       = [](double x){ return (1.0 + x) * exp(cos(x)); };
    double dokladna = -116.659583194914;
    double a = -14.0, b = 1.0;

    cout << "Calka -14 do 1: (1+x)*exp(cos(x)) dx\n";
    cout << "  Wartosc dokladna:  " << dokladna << "\n\n";

    auto pokaz = [&](const string& nazwa, double wartosc) {
        cout << "  " << setw(28) << left << nazwa
             << wartosc << "   blad: " << fabs(wartosc - dokladna) << "\n";
    };

    pokaz("Metoda trapezow (n=1000):",  metoda_trapezow(g, a, b, 1000));
    pokaz("Metoda Simpsona (n=1000):", metoda_simpsona(g, a, b, 1000));
    pokaz("Gauss-Legendre 2w (n=50):", gauss_legendre(g, a, b, 2, 50));
    pokaz("Gauss-Legendre 3w (n=50):", gauss_legendre(g, a, b, 3, 50));
    pokaz("Gauss-Legendre 4w (n=50):", gauss_legendre(g, a, b, 4, 50));

    return 0;
}
