#pragma once
/**
 * @file aproksymacja.h
 * @brief Aproksymacja srednio-kwadratowa i schemat Hornera
 *
 * Zawiera implementacje:
 *  - Aproksymacji wielomianowej metoda najmniejszych kwadratow
 *  - Schematu Hornera do ewaluacji wielomianow
 *
 * @example
 *   // Aproksymacja funkcji f na [-2.5, 1.0] wielomianem stopnia 5
 *   auto wspolczynniki = aproksymacja_mnk(f, -2.5, 1.0, 5);
 *   double wartosc = horner(wspolczynniki, 0.0);
 */

#include <vector>
#include <functional>

namespace biblioteka_mn {

/**
 * @brief Aproksymacja funkcji wielomianem stopnia n metoda najmniejszych kwadratow.
 *
 * Buduje uklad rownan normalnych (macierz Grama) i rozwiazuje go eliminacja
 * Gaussa. Calki wymagane do budowy macierzy i wektora prawej strony obliczane
 * sa regula Simpsona.
 *
 * @param f       Aproksymowana funkcja
 * @param a       Lewy koniec przedzialu
 * @param b       Prawy koniec przedzialu
 * @param stopien Stopien wielomianu aproksymujacego
 * @return        Wektor wspolczynnikow [a0, a1, ..., a_n] (od wyrazu wolnego)
 */
std::vector<double> aproksymacja_mnk(std::function<double(double)> f,
                                      double a, double b, int stopien);

/**
 * @brief Schemat Hornera - ewaluacja wielomianu w punkcie x.
 *
 * Wspolczynniki podawane od wyrazu wolnego: wspolczynniki[0] + wspolczynniki[1]*x + ...
 *
 * @param wspolczynniki  Wektor wspolczynnikow (indeks 0 = wyraz wolny)
 * @param x              Argument
 * @return               Wartosc wielomianu
 */
double horner(const std::vector<double>& wspolczynniki, double x);

} // namespace biblioteka_mn
