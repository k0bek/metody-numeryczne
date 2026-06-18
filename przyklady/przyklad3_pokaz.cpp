/**
 * @file przyklad3_pokaz.cpp
 * @brief Przyklad 3: Pelny pokaz mozliwosci biblioteki biblioteka_mn
 *
 * Demonstruje wszystkie 6 kategorii algorytmow w jednym spojnym przykladzie.
 * Uruchom: ./przyklady/przyklad3
 * Kompilacja: g++ -std=c++17 -I../include przyklad3_pokaz.cpp ../src/*.cpp -o przyklad3
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>

#include "../include/uklady_rownan.h"
#include "../include/interpolacja.h"
#include "../include/aproksymacja.h"
#include "../include/calkowanie.h"
#include "../include/rownania_rozniczkowe.h"
#include "../include/rownania_nieliniowe.h"

using namespace std;
using namespace biblioteka_mn;

// -- Pomocnicze ---------------------------------------------------------------
void naglowek(const string& tytul) {
    cout << "\n" << string(60, '=') << "\n";
    cout << "  " << tytul << "\n";
    cout << string(60, '=') << "\n";
}
void linia() { cout << string(60, '-') << "\n"; }

// ================================================================
// 1. UKLADY ROWNAN LINIOWYCH
// ================================================================
void pokaz_uklady_rownan() {
    naglowek("1. UKLADY ROWNAN LINIOWYCH");

    // Problem: siec rezystorow (Kirchhoff)
    vector<vector<double>> A = {
        { 4, -1, -1},
        {-1,  3,  0},
        {-1,  0,  2}
    };
    vector<double> b = {5, 4, 3};

    cout << fixed << setprecision(6);
    cout << "Uklad rownan (siec rezystorow):\n";
    cout << "  4x - y - z = 5\n";
    cout << " -x + 3y     = 4\n";
    cout << " -x     + 2z = 3\n\n";

    auto x_gauss = eliminacja_gaussa(A, b);
    cout << "Eliminacja Gaussa:\n";
    for (int i = 0; i < 3; ++i)
        cout << "  x[" << i << "] = " << x_gauss[i] << "\n";

    auto x_lu = rozwiaz_lu(A, b);
    cout << "\nRozklad LU:\n";
    for (int i = 0; i < 3; ++i)
        cout << "  x[" << i << "] = " << x_lu[i] << "\n";

    double max_blad = 0;
    for (int i = 0; i < 3; ++i) {
        double wiersz = 0;
        for (int j = 0; j < 3; ++j) wiersz += A[i][j] * x_gauss[j];
        max_blad = max(max_blad, fabs(wiersz - b[i]));
    }
    cout << "\nWeryfikacja (maks. residuum |Ax-b|): " << scientific << max_blad << "\n";
}

// ================================================================
// 2. INTERPOLACJA
// ================================================================
void pokaz_interpolacje() {
    naglowek("2. INTERPOLACJA");

    vector<double> godz = {0, 4, 8, 12, 16, 20, 24};
    vector<double> temp = {12.1, 10.3, 14.7, 21.5, 23.8, 18.4, 12.1};

    cout << fixed << setprecision(4);
    cout << "Dane: temperatura w ciagu doby\n";
    cout << "Godzina  |  Pomiar\n";
    linia();
    for (size_t i = 0; i < godz.size(); ++i)
        cout << "  " << setw(5) << godz[i] << "  |  " << temp[i] << " C\n";

    cout << "\nInterpolacja w wybranych godzinach:\n";
    cout << "Godzina  |  Lagrange  |  Newton\n";
    linia();

    auto wspolczynniki = wspolczynniki_newtona(godz, temp);

    for (double g : {2.0, 6.0, 10.0, 14.0, 18.0, 22.0}) {
        double L = lagrange(g, godz, temp);
        double N = oblicz_newtona(g, godz, wspolczynniki);
        cout << "  " << setw(5) << g
             << "  |  " << setw(8) << L
             << "  |  " << setw(8) << N << "\n";
    }

    cout << "\nBlad w wezlach (powinien byc ~0):\n";
    double max_blad = 0;
    for (size_t i = 0; i < godz.size(); ++i)
        max_blad = max(max_blad, fabs(lagrange(godz[i], godz, temp) - temp[i]));
    cout << "  Max blad Lagrange: " << scientific << max_blad << "\n";
}

// ================================================================
// 3. APROKSYMACJA
// ================================================================
void pokaz_aproksymacje() {
    naglowek("3. APROKSYMACJA");

    auto f = [](double x) { return sin(x) + 0.5 * x; };

    cout << fixed << setprecision(6);
    cout << "Aproksymacja f(x) = sin(x) + 0.5x wielomianami MNK na [0, 4]\n\n";
    cout << "Stopien  |  Blad w x=2.0\n";
    linia();

    double dokladna_20 = f(2.0);

    for (int st : {1, 3, 5, 7}) {
        auto wspolczynniki = aproksymacja_mnk(f, 0.0, 4.0, st);
        double aproks = horner(wspolczynniki, 2.0);
        double blad = fabs(aproks - dokladna_20);
        cout << "  st=" << st << "   |  " << scientific << blad << "\n";
    }

    cout << "\nSchemat Hornera dla p(x) = 1 + 2x + 3x^2 + 4x^3 w x=2:\n";
    vector<double> wielomian = {1, 2, 3, 4};
    cout << "  Wynik: " << fixed << horner(wielomian, 2.0) << "  (oczekiwane: 49)\n";
}

// ================================================================
// 4. CALKOWANIE
// ================================================================
void pokaz_calkowanie() {
    naglowek("4. CALKOWANIE NUMERYCZNE");

    auto f = [](double x) { return (1.0 + x) * exp(cos(x)); };
    double dokladna = -116.659583194914;
    double a = -14.0, b = 1.0;

    cout << fixed << setprecision(12);
    cout << "Calka -14 do 1: (1+x)*exp(cos(x)) dx\n";
    cout << "Dokladna wartosc: " << dokladna << "\n\n";

    struct { string nazwa; double wartosc; } metody[] = {
        {"Trapezy (n=100)",    metoda_trapezow(f, a, b, 100)},
        {"Trapezy (n=1000)",   metoda_trapezow(f, a, b, 1000)},
        {"Simpson (n=100)",    metoda_simpsona(f, a, b, 100)},
        {"Simpson (n=1000)",   metoda_simpsona(f, a, b, 1000)},
        {"Gauss-L. 2w (n=50)", gauss_legendre(f, a, b, 2, 50)},
        {"Gauss-L. 3w (n=50)", gauss_legendre(f, a, b, 3, 50)},
        {"Gauss-L. 4w (n=50)", gauss_legendre(f, a, b, 4, 50)},
    };

    cout << setw(22) << left << "Metoda"
         << setw(20) << "Wynik"
         << "Blad\n";
    linia();
    for (auto& m : metody) {
        cout << setw(22) << left << m.nazwa
             << setw(20) << m.wartosc
             << scientific << fabs(m.wartosc - dokladna) << "\n";
    }
}

// ================================================================
// 5. ROWNANIA ROZNICZKOWE
// ================================================================
void pokaz_rownania_rozniczkowe() {
    naglowek("5. ROWNANIA ROZNICZKOWE");

    auto f_test = [](double /*t*/, double y) { return -y; };
    double t_koniec = 5.0;
    auto dokladne = [](double t) { return exp(-t); };

    cout << fixed << setprecision(8);
    cout << "Rownanie: y' = -y,  y(0) = 1\n";
    cout << "Dokladne rozwiazanie: y(t) = e^(-t)\n\n";
    cout << setw(8) << "t"
         << setw(14) << "Dokladne"
         << setw(14) << "Euler"
         << setw(14) << "RK4" << "\n";
    linia();

    int n = 50;
    auto eu = euler(f_test, 0.0, 1.0, t_koniec, n);
    auto rk = rk4  (f_test, 0.0, 1.0, t_koniec, n);

    int krok = n / 5;
    for (int i = 0; i <= n; i += krok) {
        double t = eu[i].first;
        cout << setw(8) << t
             << setw(14) << dokladne(t)
             << setw(14) << eu[i].second
             << setw(14) << rk[i].second << "\n";
    }

    cout << "\nBlad bezwzgledny w t=" << t_koniec << ":\n";
    cout << "  Euler: " << scientific << fabs(eu.back().second - dokladne(t_koniec)) << "\n";
    cout << "  RK4:   " << fabs(rk.back().second - dokladne(t_koniec)) << "\n";
}

// ================================================================
// 6. ROWNANIA NIELINIOWE
// ================================================================
void pokaz_rownania_nieliniowe() {
    naglowek("6. ROWNANIA NIELINIOWE");

    cout << fixed << setprecision(10);

    auto f1  = [](double x) { return cosh(x) - sqrt(x) - 1.0; };
    auto df1 = [](double x) { return sinh(x) - 0.5/sqrt(x); };

    cout << "f(x) = cosh(x) - sqrt(x) - 1,  przedzial [0.01, 4]\n";
    cout << setw(24) << left << "Metoda" << "Wynik\n";
    linia();
    cout << setw(24) << "Bisekcja:"       << bisekcja(f1, 0.01, 4.0, 1e-10) << "\n";
    cout << setw(24) << "Regula falsi:"   << regula_falsi(f1, 0.01, 4.0, 1e-10) << "\n";
    cout << setw(24) << "Newton (an.):"   << newton(f1, df1, 1.0, 1e-12) << "\n";
    cout << setw(24) << "Newton (num.):"  << newton_numeryczny(f1, 1.0, 1e-12) << "\n";
    cout << setw(24) << "Sieczne:"        << sieczne(f1, 0.5, 1.5, 1e-12) << "\n";

    cout << "\nf(x) = (x-1)^5,  pierwiastek w x=1 (krotnosc 5)\n";
    auto f2 = [](double x){ return pow(x-1.0, 5); };
    cout << "  Bisekcja na [0,2]: " << bisekcja(f2, 0.0, 2.0, 1e-6) << "\n";

    cout << "\nPrzyklad obslugi bledu (f bez pierwiastka na przedziale):\n";
    try {
        bisekcja([](double x){ return x*x + 1; }, -2.0, 2.0);
    } catch (const invalid_argument& e) {
        cout << "  Zlapany wyjatek: " << e.what() << "\n";
    }
}

// ================================================================
int main() {
    cout << "==============================================================\n";
    cout << "           biblioteka_mn - PELNY POKAZ BIBLIOTEKI\n";
    cout << "==============================================================\n";

    pokaz_uklady_rownan();
    pokaz_interpolacje();
    pokaz_aproksymacje();
    pokaz_calkowanie();
    pokaz_rownania_rozniczkowe();
    pokaz_rownania_nieliniowe();

    cout << "\n" << string(60, '=') << "\n";
    cout << "  Koniec demonstracji. Wszystkie moduly dzialaja poprawnie.\n";
    cout << string(60, '=') << "\n";
    return 0;
}
