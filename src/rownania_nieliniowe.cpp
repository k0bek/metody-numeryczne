#include "rownania_nieliniowe.h"
#include <cmath>
#include <stdexcept>

namespace biblioteka_mn {

double bisekcja(std::function<double(double)> f,
                double a, double b, double eps, int max_iter) {
    if (f(a) * f(b) > 0)
        throw std::invalid_argument(
            "bisekcja: f(a) i f(b) maja ten sam znak - brak gwarantowanego pierwiastka");

    double c = a;
    for (int i = 0; i < max_iter; ++i) {
        c = (a + b) / 2.0;
        double fc = f(c);
        if (std::abs(fc) < eps) break;
        if (f(a) * fc < 0)
            b = c;
        else
            a = c;
    }
    return c;
}

double regula_falsi(std::function<double(double)> f,
                    double a, double b, double eps, int max_iter) {
    if (f(a) * f(b) > 0)
        throw std::invalid_argument(
            "regula_falsi: f(a) i f(b) maja ten sam znak - brak gwarantowanego pierwiastka");

    double x = a;
    for (int i = 0; i < max_iter; ++i) {
        double fa = f(a), fb = f(b);
        x = a - fa * (b - a) / (fb - fa);
        double fx = f(x);
        if (std::abs(fx) < eps) break;
        if (fa * fx < 0)
            b = x;
        else
            a = x;
    }
    return x;
}

double newton(std::function<double(double)> f,
              std::function<double(double)> df,
              double x0, double eps, int max_iter) {
    double x = x0;
    for (int i = 0; i < max_iter; ++i) {
        double dfx = df(x);
        if (std::abs(dfx) < 1e-12)
            throw std::runtime_error("newton: pochodna bliska zeru");
        double x_nowe = x - f(x) / dfx;
        if (std::abs(x_nowe - x) < eps) { x = x_nowe; break; }
        x = x_nowe;
    }
    return x;
}

double newton_numeryczny(std::function<double(double)> f,
                         double x0, double eps, int max_iter) {
    // Pochodna centralna: f'(x) ~ (f(x+h) - f(x-h)) / (2h)
    auto df = [&](double x) {
        const double h = 1e-6;
        return (f(x + h) - f(x - h)) / (2.0 * h);
    };
    return newton(f, df, x0, eps, max_iter);
}

double sieczne(std::function<double(double)> f,
              double x0, double x1, double eps, int max_iter) {
    double x_prev = x0, x_aktualne = x1, x_nowe;
    for (int i = 0; i < max_iter; ++i) {
        double f_prev = f(x_prev), f_akt = f(x_aktualne);
        double mianownik = f_akt - f_prev;
        if (std::abs(mianownik) < 1e-12)
            throw std::runtime_error("sieczne: dzielenie przez zero");
        x_nowe = (f_akt * x_prev - f_prev * x_aktualne) / mianownik;
        if (std::abs(x_nowe - x_aktualne) < eps) { x_aktualne = x_nowe; break; }
        x_prev = x_aktualne;
        x_aktualne = x_nowe;
    }
    return x_aktualne;
}

} // namespace biblioteka_mn
