#pragma once
/**
 * @file calkowanie.h
 * @brief Numeryczne calkowanie
 *
 * Zawiera implementacje:
 *  - Metody trapezow
 *  - Metody Simpsona
 *  - Zlozonej kwadratury Gaussa-Legendre'a (2, 3, 4 wezly)
 *
 * @example
 *   auto f = [](double x){ return x*x; };
 *   double I = metoda_simpsona(f, 0.0, 1.0, 100);   // ~0.333333
 *   double J = gauss_legendre(f, 0.0, 1.0, 3, 10);  // ~0.333333
 */

#include <functional>
#include <vector>

namespace biblioteka_mn {

/**
 * @brief Metoda trapezow.
 *
 * @param f  Funkcja podcalkowa
 * @param a  Dolna granica calkowania
 * @param b  Gorna granica calkowania
 * @param n  Liczba podprzedzialow (im wieksza, tym dokladniej)
 * @return   Przyblizona wartosc calki
 */
double metoda_trapezow(std::function<double(double)> f, double a, double b, int n);

/**
 * @brief Metoda Simpsona (zlozona regula 1/3 Simpsona).
 *
 * @param f  Funkcja podcalkowa
 * @param a  Dolna granica calkowania
 * @param b  Gorna granica calkowania
 * @param n  Liczba podprzedzialow (zaokraglana do parzystej)
 * @return   Przyblizona wartosc calki
 */
double metoda_simpsona(std::function<double(double)> f, double a, double b, int n);

/**
 * @brief Zlozona kwadratura Gaussa-Legendre'a.
 *
 * Calkuje funkcje f na [a,b] dzielac przedzial na podprzedzialy podprzedzialy,
 * a kazdy calkuje kwadratura GL z liczba_wezlow wezlami (2, 3 lub 4).
 *
 * @param f              Funkcja podcalkowa
 * @param a              Dolna granica
 * @param b              Gorna granica
 * @param liczba_wezlow  Liczba wezlow GL: 2, 3 lub 4
 * @param podprzedzialy  Liczba podprzedzialow zlozonej kwadratury
 * @return               Przyblizona wartosc calki
 * @throws std::invalid_argument gdy liczba_wezlow nie nalezy do {2,3,4}
 */
double gauss_legendre(std::function<double(double)> f,
                      double a, double b,
                      int liczba_wezlow,
                      int podprzedzialy = 50);

} // namespace biblioteka_mn
