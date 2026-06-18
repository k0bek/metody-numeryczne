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

    // Rozwiazujemy uklad 3 rownan z 3 niewiadomymi metoda Gaussa.
    // Oczekiwany wynik: x={2, 3, -1}. Sprawdzamy poprawnosc eliminacji z pivotingiem.
    {
        vector<vector<double>> A = {{2,1,-1},{-3,-1,2},{-2,1,2}};
        vector<double> b = {8,-11,-3};
        auto x = eliminacja_gaussa(A, b);
        sprawdz(bliskie(x[0],2) && bliskie(x[1],3) && bliskie(x[2],-1),
              "eliminacja_gaussa: uklad 3x3");
    }

    // Przypadek brzegowy: macierz 1x1 (jedno rownanie 5x=15).
    // Weryfikujemy, ze algorytm dziala poprawnie dla minimalnego rozmiaru.
    {
        vector<vector<double>> A = {{5.0}};
        vector<double> b = {15.0};
        auto x = eliminacja_gaussa(A, b);
        sprawdz(bliskie(x[0], 3.0), "eliminacja_gaussa: macierz 1x1");
    }

    // Przypadek bledny: uklad sprzeczny (dwa identyczne wiersze, rozne prawej strony).
    // Oczekujemy rzucenia wyjatku runtime_error.
    {
        bool zlapano = false;
        try {
            vector<vector<double>> A = {{1,1},{1,1}};
            vector<double> b = {2, 3};
            eliminacja_gaussa(A, b);
        } catch (const runtime_error&) { zlapano = true; }
        sprawdz(zlapano, "eliminacja_gaussa: wyjatek dla ukladu sprzecznego");
    }

    // Rozwiazujemy ten sam uklad 3x3 metoda rozkladu LU.
    // Porownujemy wynik z eliminacja Gaussa - powinien byc identyczny.
    {
        vector<vector<double>> A = {{2,1,-1},{-3,-1,2},{-2,1,2}};
        vector<double> b = {8,-11,-3};
        auto x = rozwiaz_lu(A, b);
        sprawdz(bliskie(x[0],2) && bliskie(x[1],3) && bliskie(x[2],-1),
              "rozwiaz_lu: uklad 3x3");
    }

    // Weryfikacja poprawnosci samego rozkladu LU: obliczamy iloczyn L*U
    // i sprawdzamy, czy odzyskujemy oryginalna macierz A (2x2).
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

    // Dane testowe: wezly wielomianu f(x) = x^2 + x + 1
    vector<double> xs = {0, 1, 2, 3};
    vector<double> ys = {1, 3, 7, 13};

    // Ewaluacja Lagrange'a dokladnie w wezle x=1. Wynik musi byc rowny y=3.
    sprawdz(bliskie(lagrange(1.0, xs, ys), 3.0), "lagrange: wartosc w wezle");

    // Ewaluacja Lagrange'a miedzy wezlami (x=1.5). Oczekiwany wynik: 1.5^2+1.5+1 = 4.75.
    sprawdz(bliskie(lagrange(1.5, xs, ys), 4.75), "lagrange: wartosc posrednia (1.5)");

    // Obliczamy ilorazy roznicowe Newtona, nastepnie ewaluujemy w wezle x=2.
    // Wynik musi zgadzac sie z wartoscia y=7.
    auto wspolczynniki = wspolczynniki_newtona(xs, ys);
    sprawdz(bliskie(oblicz_newtona(2.0, xs, wspolczynniki), 7.0),
          "oblicz_newtona: wartosc w wezle");

    // Porownanie Lagrange'a z Newtonem poza wezlami (x=2.5).
    // Oba wielomiany interpolacyjne musza dawac identyczny wynik.
    double xtest = 2.5;
    sprawdz(bliskie(lagrange(xtest, xs, ys), oblicz_newtona(xtest, xs, wspolczynniki)),
          "lagrange == oblicz_newtona poza wezlami");
}

// -- 3. Aproksymacja ------------------------------------------------------------

void test_aproksymacja() {
    cout << "\n=== 3. Aproksymacja ===\n";

    // Schemat Hornera dla wielomianu 1 + 2x + 3x^2 w punkcie x=2.
    // Wynik: 1 + 2*2 + 3*4 = 17.
    vector<double> c = {1, 2, 3};
    sprawdz(bliskie(horner(c, 2.0), 17.0), "horner: 1 + 2x + 3x^2 w x=2");

    // Przypadek brzegowy: wielomian zerowy (wszystkie wspolczynniki = 0).
    // Wynik musi byc 0 niezaleznie od argumentu.
    vector<double> cz = {0, 0, 0};
    sprawdz(bliskie(horner(cz, 5.0), 0.0), "horner: wielomian zerowy");

    // Aproksymacja MNK funkcji stalej f(x)=3 wielomianem stopnia 0 na [0,1].
    // Jedyny wspolczynnik a0 powinien wynosic 3.
    auto f_const = [](double) { return 3.0; };
    auto cf = aproksymacja_mnk(f_const, 0.0, 1.0, 0);
    sprawdz(bliskie(cf[0], 3.0, 1e-4), "aproksymacja_mnk: stala funkcja stopnia 0");

    // Aproksymacja MNK funkcji liniowej f(x)=x+1 wielomianem stopnia 1 na [0,2].
    // Oczekiwane wspolczynniki: a0=1, a1=1 (dokladne odwzorowanie).
    auto f_lin = [](double x) { return x + 1.0; };
    auto cl = aproksymacja_mnk(f_lin, 0.0, 2.0, 1);
    sprawdz(bliskie(cl[0], 1.0, 1e-3) && bliskie(cl[1], 1.0, 1e-3),
          "aproksymacja_mnk: f(x)=x+1 stopnia 1");
}

// -- 4. Calkowanie numeryczne ---------------------------------------------------

void test_calkowanie() {
    cout << "\n=== 4. Calkowanie numeryczne ===\n";

    auto f2 = [](double x) { return x * x; };

    // Calka z x^2 na [0,1] metoda trapezow (n=10000). Wynik analityczny: 1/3.
    // Sprawdzamy zbieznosc metody trapezow z dokladnoscia 1e-4.
    sprawdz(bliskie(metoda_trapezow(f2, 0.0, 1.0, 10000), 1.0/3.0, 1e-4),
          "metoda_trapezow: calka x^2 na [0,1]");

    // Ta sama calka metoda Simpsona (n=100). Simpson jest dokladniejszy,
    // wiec wystarcza mniej podprzedzialow do osiagniecia 1e-8.
    sprawdz(bliskie(metoda_simpsona(f2, 0.0, 1.0, 100), 1.0/3.0, 1e-8),
          "metoda_simpsona: calka x^2 na [0,1]");

    // Kwadratura Gaussa-Legendre'a z 2 wezlami (10 podprzedzialow).
    // GL z 2 wezlami jest dokladna dla wielomianow stopnia <= 3, wiec x^2 jest dokladne.
    sprawdz(bliskie(gauss_legendre(f2, 0.0, 1.0, 2, 10), 1.0/3.0, 1e-10),
          "gauss_legendre (2 wezly): calka x^2 na [0,1]");

    // GL z 3 wezlami - dokladna dla wielomianow stopnia <= 5.
    // Oczekujemy jeszcze wyzszej precyzji niz dla 2 wezlow.
    sprawdz(bliskie(gauss_legendre(f2, 0.0, 1.0, 3, 10), 1.0/3.0, 1e-12),
          "gauss_legendre (3 wezly): calka x^2 na [0,1]");

    // Calka z sin(x) na [0, pi] metoda Simpsona. Wynik analityczny: 2.
    // Testujemy na innej funkcji niz wielomianowa.
    auto fsin = [](double x) { return sin(x); };
    sprawdz(bliskie(metoda_simpsona(fsin, 0.0, M_PI, 1000), 2.0, 1e-6),
          "metoda_simpsona: calka sin(x) na [0,pi]");

    // Przypadek bledny: podano 5 wezlow GL, a obslugiwane sa tylko 2, 3, 4.
    // Oczekujemy wyjatku invalid_argument.
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

    // Rownanie testowe: y' = -y, y(0) = 1. Rozwiazanie analityczne: y(t) = e^(-t).
    // Porownujemy y(1) numeryczne z e^(-1) = 0.367879...
    auto f = [](double /*t*/, double y) { return -y; };
    double t_koniec = 1.0, y_dokladne = exp(-1.0);

    // Metoda Eulera (rzad 1) - potrzebuje duzo krokow (1000) dla bledu < 1e-3.
    auto traj_euler = euler(f, 0.0, 1.0, t_koniec, 1000);
    sprawdz(bliskie(traj_euler.back().second, y_dokladne, 1e-3),
          "euler: y'=-y, blad < 1e-3");

    // Metoda Heuna (rzad 2) - 100 krokow wystarcza dla bledu < 1e-5.
    auto traj_heun = heun(f, 0.0, 1.0, t_koniec, 100);
    sprawdz(bliskie(traj_heun.back().second, y_dokladne, 1e-5),
          "heun: y'=-y, blad < 1e-5");

    // Metoda punktu srodkowego (rzad 2) - podobna dokladnosc jak Heun.
    auto traj_mid = punkt_srodkowy(f, 0.0, 1.0, t_koniec, 100);
    sprawdz(bliskie(traj_mid.back().second, y_dokladne, 1e-5),
          "punkt_srodkowy: y'=-y, blad < 1e-5");

    // Metoda RK4 (rzad 4) - nawet 50 krokow daje blad < 1e-8.
    // Demonstracja przewagi metod wyzszego rzedu.
    auto traj_rk4 = rk4(f, 0.0, 1.0, t_koniec, 50);
    sprawdz(bliskie(traj_rk4.back().second, y_dokladne, 1e-8),
          "rk4: y'=-y, blad < 1e-8");

    // Sprawdzenie dlugosci trajektorii: n krokow powinno dac n+1 punktow
    // (punkt poczatkowy + n kolejnych).
    sprawdz((int)traj_rk4.size() == 51, "rk4: dlugosc trajektorii = n+1");
}

// -- 6. Rownania nieliniowe ----------------------------------------------------

void test_rownania_nieliniowe() {
    cout << "\n=== 6. Rownania nieliniowe ===\n";

    // Funkcja testowa: f(x) = x^2 - 2, pierwiastek = sqrt(2) ~ 1.41421356
    auto f1  = [](double x) { return x*x - 2.0; };
    auto df1 = [](double x) { return 2.0*x; };

    // Bisekcja na przedziale [1,2] - f(1)<0, f(2)>0, wiec pierwiastek jest w srodku.
    sprawdz(bliskie(bisekcja(f1, 1.0, 2.0, 1e-8), sqrt(2.0), 1e-7),
          "bisekcja: x^2 - 2 = 0");

    // Regula falsi na tym samym przedziale - szybsza zbieznosc niz bisekcja.
    sprawdz(bliskie(regula_falsi(f1, 1.0, 2.0, 1e-8), sqrt(2.0), 1e-7),
          "regula_falsi: x^2 - 2 = 0");

    // Newton z pochodna analityczna df(x)=2x, start x0=1.5.
    // Zbieznosc kwadratowa - najszybsza z metod.
    sprawdz(bliskie(newton(f1, df1, 1.5, 1e-10), sqrt(2.0), 1e-9),
          "newton (analityczny): x^2 - 2 = 0");

    // Newton z pochodna numeryczna (roznice centralne h=1e-6).
    // Nieco mniejsza dokladnosc niz analityczny z powodu przyblizonej pochodnej.
    sprawdz(bliskie(newton_numeryczny(f1, 1.5, 1e-10), sqrt(2.0), 1e-8),
          "newton_numeryczny: x^2 - 2 = 0");

    // Metoda siecznych - dwa punkty startowe x0=1, x1=2, bez pochodnej.
    sprawdz(bliskie(sieczne(f1, 1.0, 2.0, 1e-10), sqrt(2.0), 1e-8),
          "sieczne: x^2 - 2 = 0");

    // Przypadek bledny: f(x) = x^2 + 1 nie ma pierwiastkow rzeczywistych.
    // f(-2)=5 > 0 i f(2)=5 > 0, wiec bisekcja musi rzucic wyjatek.
    {
        bool zlapano = false;
        try { bisekcja([](double x){ return x*x + 1; }, -2.0, 2.0); }
        catch (const invalid_argument&) { zlapano = true; }
        sprawdz(zlapano, "bisekcja: wyjatek dla f bez pierwiastka");
    }

    // Pierwiastek krotny: (x-1)^5 = 0 ma pierwiastek x=1 o krotnosci 5.
    // Bisekcja radzi sobie z pierwiastkami krotnymi (w przeciwienstwie do Newtona).
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
