#include "aproksymacja.h"
#include "uklady_rownan.h"
#include <cmath>
#include <stdexcept>

namespace biblioteka_mn {

// Pomocnicza: calka z x^n na [a,b] - wzor analityczny
static double calka_xn(double a, double b, int n) {
    return (std::pow(b, n + 1) - std::pow(a, n + 1)) / (n + 1);
}

// Pomocnicza: calka z f(x)*x^k na [a,b] - regula Simpsona
static double calka_f_xk(std::function<double(double)> f,
                         double a, double b, int k,
                         int kroki = 5000) {
    if (kroki % 2 != 0) kroki++;
    double h = (b - a) / kroki;
    auto g = [&](double x) {
        return f(x) * std::pow(x, k);
    };
    double suma = g(a) + g(b);
    for (int i = 1; i < kroki; ++i) {
        double x = a + i * h;
        suma += (i % 2 == 0 ? 2.0 : 4.0) * g(x);
    }
    return (h / 3.0) * suma;
}

std::vector<double> aproksymacja_mnk(std::function<double(double)> f,
                                      double a, double b, int stopien) {
    int rozmiar = stopien + 1;
    std::vector<std::vector<double>> A(rozmiar, std::vector<double>(rozmiar));
    std::vector<double> prawa_strona(rozmiar);

    // Budowanie macierzy Grama: A[i][j] = calka(x^(i+j), a, b)
    for (int i = 0; i <= stopien; ++i) {
        for (int j = 0; j <= stopien; ++j)
            A[i][j] = calka_xn(a, b, i + j);
        prawa_strona[i] = calka_f_xk(f, a, b, i);
    }

    return eliminacja_gaussa(A, prawa_strona); // [a0, a1, ..., a_stopien]
}

double horner(const std::vector<double>& wspolczynniki, double x) {
    // wspolczynniki[0] = wyraz wolny, wspolczynniki[n] = wspolczynnik przy x^n
    int n = (int)wspolczynniki.size() - 1;
    double wynik = wspolczynniki[n];
    for (int i = n - 1; i >= 0; --i)
        wynik = wspolczynniki[i] + x * wynik;
    return wynik;
}

} // namespace biblioteka_mn
