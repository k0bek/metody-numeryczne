#include "calkowanie.h"
#include <cmath>
#include <stdexcept>
#include <vector>

namespace biblioteka_mn {

double metoda_trapezow(std::function<double(double)> f, double a, double b, int n) {
    double h = (b - a) / n;
    double suma = 0.5 * (f(a) + f(b));
    for (int i = 1; i < n; ++i)
        suma += f(a + i * h);
    return suma * h;
}

double metoda_simpsona(std::function<double(double)> f, double a, double b, int n) {
    if (n % 2 != 0) n++;
    double h = (b - a) / n;
    double suma = f(a) + f(b);
    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        suma += (i % 2 != 0 ? 4.0 : 2.0) * f(x);
    }
    return suma * (h / 3.0);
}

double gauss_legendre(std::function<double(double)> f,
                      double a, double b,
                      int liczba_wezlow,
                      int podprzedzialy) {
    // Wezly i wagi kwadratury GL na [-1, 1]
    std::vector<double> wezly, wagi;

    if (liczba_wezlow == 2) {
        double t = 1.0 / std::sqrt(3.0);
        wezly = {-t, t};
        wagi  = {1.0, 1.0};
    } else if (liczba_wezlow == 3) {
        double t = std::sqrt(3.0 / 5.0);
        wezly = {-t, 0.0, t};
        wagi  = {5.0/9.0, 8.0/9.0, 5.0/9.0};
    } else if (liczba_wezlow == 4) {
        double t1 = std::sqrt((3.0/7.0) + (2.0/7.0)*std::sqrt(6.0/5.0));
        double t2 = std::sqrt((3.0/7.0) - (2.0/7.0)*std::sqrt(6.0/5.0));
        double w1 = (18.0 - std::sqrt(30.0)) / 36.0;
        double w2 = (18.0 + std::sqrt(30.0)) / 36.0;
        wezly = {-t1, -t2, t2, t1};
        wagi  = {w1, w2, w2, w1};
    } else {
        throw std::invalid_argument("gauss_legendre: liczba_wezlow musi byc 2, 3 lub 4");
    }

    double calka_total = 0.0;
    double szerokosc_podprzedzialu = (b - a) / podprzedzialy;

    for (int i = 0; i < podprzedzialy; ++i) {
        double m_a = a + i * szerokosc_podprzedzialu;
        double m_b = m_a + szerokosc_podprzedzialu;
        double wspolczynnik = (m_b - m_a) / 2.0;
        double przesuniecie = (m_a + m_b) / 2.0;
        double suma_lokalna = 0.0;
        for (std::size_t j = 0; j < wezly.size(); ++j)
            suma_lokalna += wagi[j] * f(wspolczynnik * wezly[j] + przesuniecie);
        calka_total += wspolczynnik * suma_lokalna;
    }
    return calka_total;
}

} // namespace biblioteka_mn
