#include "interpolacja.h"
#include <vector>

namespace biblioteka_mn {

double lagrange(double x,
                const std::vector<double>& wezly_x,
                const std::vector<double>& wezly_y) {
    double wynik = 0.0;
    int n = (int)wezly_x.size();
    for (int i = 0; i < n; ++i) {
        double li = 1.0;
        for (int j = 0; j < n; ++j) {
            if (i != j)
                li *= (x - wezly_x[j]) / (wezly_x[i] - wezly_x[j]);
        }
        wynik += wezly_y[i] * li;
    }
    return wynik;
}

std::vector<double> wspolczynniki_newtona(const std::vector<double>& wezly_x,
                                           const std::vector<double>& wezly_y) {
    std::vector<double> a = wezly_y; // kopia wartosci, modyfikowana w miejscu
    int n = (int)wezly_x.size();
    for (int j = 1; j < n; ++j) {
        for (int i = n - 1; i >= j; --i) {
            // iloraz roznicowy: [x_{i-j}, ..., x_i]
            a[i] = (a[i] - a[i - 1]) / (wezly_x[i] - wezly_x[i - j]);
        }
    }
    return a; // a[k] = iloraz roznicowy k-tego rzedu
}

double oblicz_newtona(double x,
                      const std::vector<double>& wezly_x,
                      const std::vector<double>& wspolczynniki) {
    double wynik = wspolczynniki[0];
    double iloczyn = 1.0;
    for (std::size_t i = 1; i < wspolczynniki.size(); ++i) {
        iloczyn *= (x - wezly_x[i - 1]);
        wynik += wspolczynniki[i] * iloczyn;
    }
    return wynik;
}

} // namespace biblioteka_mn
