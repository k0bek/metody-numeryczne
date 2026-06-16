#pragma once
/**
 * @file interpolacja.h
 * @brief Interpolacja wielomianowa
 *
 * Zawiera implementacje:
 *  - Interpolacji Lagrange'a
 *  - Interpolacji Newtona (ilorazy roznicowe)
 *
 * @example
 *   vector<double> wezly_x = {0, 1, 2};
 *   vector<double> wezly_y = {1, 3, 7};
 *   double wynik = lagrange(1.5, wezly_x, wezly_y);  // ~4.75
 */

#include <vector>

namespace biblioteka_mn {

/**
 * @brief Oblicza wartosc wielomianu interpolacyjnego Lagrange'a w punkcie x.
 *
 * @param x        Punkt, w ktorym obliczamy wartosc
 * @param wezly_x  Wezly interpolacji (wartosci argumentu)
 * @param wezly_y  Wartosci funkcji w wezlach
 * @return         Wartosc wielomianu L(x)
 */
double lagrange(double x,
                const std::vector<double>& wezly_x,
                const std::vector<double>& wezly_y);

/**
 * @brief Oblicza wspolczynniki (ilorazy roznicowe) wielomianu Newtona.
 *
 * Wynik mozna przekazac do oblicz_newtona().
 *
 * @param wezly_x  Wezly interpolacji
 * @param wezly_y  Wartosci funkcji w wezlach
 * @return         Wektor wspolczynnikow a[0..n]
 */
std::vector<double> wspolczynniki_newtona(const std::vector<double>& wezly_x,
                                           const std::vector<double>& wezly_y);

/**
 * @brief Oblicza wartosc wielomianu Newtona W(x).
 *
 * @param x          Punkt ewaluacji
 * @param wezly_x    Wezly interpolacji (te same co przy obliczaniu wspolczynnikow)
 * @param wspolczynniki  Wspolczynniki zwrocone przez wspolczynniki_newtona()
 * @return           Wartosc W(x)
 */
double oblicz_newtona(double x,
                      const std::vector<double>& wezly_x,
                      const std::vector<double>& wspolczynniki);

} // namespace biblioteka_mn
