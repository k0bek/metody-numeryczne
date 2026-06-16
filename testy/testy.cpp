/**
 * @file testy.cpp
 * @brief Testy jednostkowe biblioteki biblioteka_mn
 *
 * Kazda funkcja posiada co najmniej 2 testy: przypadek poprawny i brzegowy/bledny.
 * Kompilacja: g++ -std=c++17 -I../include testy.cpp ../src/*.cpp -o testy
 */

#include <iostream>
#include <cmath>
#include <stdexcept>
#include <iomanip>

#include "../include/uklady_rownan.h"
#include "../include/interpolacja.h"
#include "../include/aproksymacja.h"
#include "../include/calkowanie.h"
#include "../include/rownania_rozniczkowe.h"
#include "../include/rownania_nieliniowe.h"

using namespace std;
using namespace biblioteka_mn;

// -- Narzedzia pomocnicze ---------------------------------------------------

static int zaliczone = 0, niezaliczone = 0;

void sprawdz(bool warunek, const string& nazwa) {
    if (warunek) {
        cout << "  [OK]   " << nazwa << "\n";
        zaliczone++;
    } else {
        cout << "  [BLAD] " << nazwa << "\n";
        niezaliczone++;
    }
}

bool bliskie(double a, double b, double tol = 1e-5) {
    return fabs(a - b) < tol;
}

// -- 1. Uklady rownan liniowych ----------------------------------------------

void test_uklady_rownan() {
    cout << "\n=== 1. Uklady rownan liniowych ===\n";

    // Test 1: klasyczny uklad 3x3  (wynik dokladny: x={2,3,-1})
    {
        vector<vector<double>> A = {{2,1,-1},{-3,-1,2},{-2,1,2}};
        vector<double> b = {8,-11,-3};
        auto x = eliminacja_gaussa(A, b);
        sprawdz(bliskie(x[0],2) && bliskie(x[1],3) && bliskie(x[2],-1),
              "eliminacja_gaussa: uklad 3x3");
    }

    // Test 2: macierz 1x1
    {
        vector<vector<double>> A = {{5.0}};
        vector<double> b = {15.0};
        auto x = eliminacja_gaussa(A, b);
        sprawdz(bliskie(x[0], 3.0), "eliminacja_gaussa: macierz 1x1");
    }

    // Test 3: blad - uklad sprzeczny
    {
        bool zlapano = false;
        try {
            vector<vector<double>> A = {{1,1},{1,1}};
            vector<double> b = {2, 3};
            eliminacja_gaussa(A, b);
        } catch (const runtime_error&) { zlapano = true; }
        sprawdz(zlapano, "eliminacja_gaussa: wyjatek dla ukladu sprzecznego");
    }

    // Test 4: rozwiaz_lu ten sam uklad 3x3
    {
        vector<vector<double>> A = {{2,1,-1},{-3,-1,2},{-2,1,2}};
        vector<double> b = {8,-11,-3};
        auto x = rozwiaz_lu(A, b);
        sprawdz(bliskie(x[0],2) && bliskie(x[1],3) && bliskie(x[2],-1),
              "rozwiaz_lu: uklad 3x3");
    }

    // Test 5: rozklad_lu - weryfikacja L*U = A
    {
        vector<vector<double>> A = {{4,3},{6,3}};
        vector<vector<double>> L, U;
        rozklad_lu(A, L, U);
        double lu00 = L[0][0]*U[0][0];
        double lu01 = L[0][0]*U[0][1];
        double lu10 = L[1][0]*U[0][0] + L[1][1]*U[1][0];
        double lu11 = L[1][0]*U[0][1] + L[1][1]*U[1][1];
        sprawdz(bliskie(lu00,4) && bliskie(lu01,3) && bliskie(lu10,6) && bliskie(lu11,3),
              "rozklad_lu: L*U = A dla macierzy 2x2");
    }
}

// -- 2. Interpolacja ----------------------------------------------------------

void test_interpolacja() {
    cout << "\n=== 2. Interpolacja ===\n";

    vector<double> xs = {0, 1, 2, 3};
    vector<double> ys = {1, 3, 7, 13}; // f(x) = x^2 + x + 1

    sprawdz(bliskie(lagrange(1.0, xs, ys), 3.0), "lagrange: wartosc w wezle");
    sprawdz(bliskie(lagrange(1.5, xs, ys), 4.75), "lagrange: wartosc posrednia (1.5)");

    auto wspolczynniki = wspolczynniki_newtona(xs, ys);
    sprawdz(bliskie(oblicz_newtona(2.0, xs, wspolczynniki), 7.0),
          "oblicz_newtona: wartosc w wezle");

    double xtest = 2.5;
    sprawdz(bliskie(lagrange(xtest, xs, ys), oblicz_newtona(xtest, xs, wspolczynniki)),
          "lagrange == oblicz_newtona poza wezlami");
}

// -- 3. Aproksymacja ------------------------------------------------------------

void test_aproksymacja() {
    cout << "\n=== 3. Aproksymacja ===\n";

    vector<double> c = {1, 2, 3};
    sprawdz(bliskie(horner(c, 2.0), 17.0), "horner: 1 + 2x + 3x^2 w x=2");

    vector<double> cz = {0, 0, 0};
    sprawdz(bliskie(horner(cz, 5.0), 0.0), "horner: wielomian zerowy");

    auto f_const = [](double) { return 3.0; };
    auto cf = aproksymacja_mnk(f_const, 0.0, 1.0, 0);
    sprawdz(bliskie(cf[0], 3.0, 1e-4), "aproksymacja_mnk: stala funkcja stopnia 0");

    auto f_lin = [](double x) { return x + 1.0; };
    auto cl = aproksymacja_mnk(f_lin, 0.0, 2.0, 1);
    sprawdz(bliskie(cl[0], 1.0, 1e-3) && bliskie(cl[1], 1.0, 1e-3),
          "aproksymacja_mnk: f(x)=x+1 stopnia 1");
}

// -- 4. Calkowanie numeryczne ---------------------------------------------------

void test_calkowanie() {
    cout << "\n=== 4. Calkowanie numeryczne ===\n";

    auto f2 = [](double x) { return x * x; };

    sprawdz(bliskie(metoda_trapezow(f2, 0.0, 1.0, 10000), 1.0/3.0, 1e-4),
          "metoda_trapezow: calka x^2 na [0,1]");

    sprawdz(bliskie(metoda_simpsona(f2, 0.0, 1.0, 100), 1.0/3.0, 1e-8),
          "metoda_simpsona: calka x^2 na [0,1]");

    sprawdz(bliskie(gauss_legendre(f2, 0.0, 1.0, 2, 10), 1.0/3.0, 1e-10),
          "gauss_legendre (2 wezly): calka x^2 na [0,1]");

    sprawdz(bliskie(gauss_legendre(f2, 0.0, 1.0, 3, 10), 1.0/3.0, 1e-12),
          "gauss_legendre (3 wezly): calka x^2 na [0,1]");

    auto fsin = [](double x) { return sin(x); };
    sprawdz(bliskie(metoda_simpsona(fsin, 0.0, M_PI, 1000), 2.0, 1e-6),
          "metoda_simpsona: calka sin(x) na [0,pi]");

    {
        bool zlapano = false;
        try { gauss_legendre(f2, 0, 1, 5); }
        catch (const invalid_argument&) { zlapano = true; }
        sprawdz(zlapano, "gauss_legendre: wyjatek dla blednej liczby wezlow");
    }
}

// -- 5. Rownania rozniczkowe ------------------------------------------------------

void test_rownania_rozniczkowe() {
    cout << "\n=== 5. Rownania rozniczkowe (RRZ) ===\n";

    auto f = [](double /*t*/, double y) { return -y; };
    double t_koniec = 1.0, y_dokladne = exp(-1.0);

    auto traj_euler = euler(f, 0.0, 1.0, t_koniec, 1000);
    sprawdz(bliskie(traj_euler.back().second, y_dokladne, 1e-3),
          "euler: y'=-y, blad < 1e-3");

    auto traj_heun = heun(f, 0.0, 1.0, t_koniec, 100);
    sprawdz(bliskie(traj_heun.back().second, y_dokladne, 1e-5),
          "heun: y'=-y, blad < 1e-5");

    auto traj_mid = punkt_srodkowy(f, 0.0, 1.0, t_koniec, 100);
    sprawdz(bliskie(traj_mid.back().second, y_dokladne, 1e-5),
          "punkt_srodkowy: y'=-y, blad < 1e-5");

    auto traj_rk4 = rk4(f, 0.0, 1.0, t_koniec, 50);
    sprawdz(bliskie(traj_rk4.back().second, y_dokladne, 1e-8),
          "rk4: y'=-y, blad < 1e-8");

    sprawdz((int)traj_rk4.size() == 51, "rk4: dlugosc trajektorii = n+1");
}

// -- 6. Rownania nieliniowe ----------------------------------------------------

void test_rownania_nieliniowe() {
    cout << "\n=== 6. Rownania nieliniowe ===\n";

    auto f1  = [](double x) { return x*x - 2.0; };
    auto df1 = [](double x) { return 2.0*x; };

    sprawdz(bliskie(bisekcja(f1, 1.0, 2.0, 1e-8), sqrt(2.0), 1e-7),
          "bisekcja: x^2 - 2 = 0");

    sprawdz(bliskie(regula_falsi(f1, 1.0, 2.0, 1e-8), sqrt(2.0), 1e-7),
          "regula_falsi: x^2 - 2 = 0");

    sprawdz(bliskie(newton(f1, df1, 1.5, 1e-10), sqrt(2.0), 1e-9),
          "newton (analityczny): x^2 - 2 = 0");

    sprawdz(bliskie(newton_numeryczny(f1, 1.5, 1e-10), sqrt(2.0), 1e-8),
          "newton_numeryczny: x^2 - 2 = 0");

    sprawdz(bliskie(sieczne(f1, 1.0, 2.0, 1e-10), sqrt(2.0), 1e-8),
          "sieczne: x^2 - 2 = 0");

    {
        bool zlapano = false;
        try { bisekcja([](double x){ return x*x + 1; }, -2.0, 2.0); }
        catch (const invalid_argument&) { zlapano = true; }
        sprawdz(zlapano, "bisekcja: wyjatek dla f bez pierwiastka");
    }

    auto f3 = [](double x){ return pow(x-1,5); };
    sprawdz(bliskie(bisekcja(f3, 0.0, 2.0, 1e-6), 1.0, 1e-5),
          "bisekcja: pierwiastek krotny (x-1)^5 = 0");
}

// -- main -----------------------------------------------------------------------

int main() {
    cout << "====================================\n";
    cout << "  TESTY JEDNOSTKOWE - biblioteka_mn\n";
    cout << "====================================\n";

    test_uklady_rownan();
    test_interpolacja();
    test_aproksymacja();
    test_calkowanie();
    test_rownania_rozniczkowe();
    test_rownania_nieliniowe();

    cout << "\n====================================\n";
    cout << "  Wynik: " << zaliczone << " OK, " << niezaliczone << " BLAD\n";
    cout << "====================================\n";
    return niezaliczone == 0 ? 0 : 1;
}
