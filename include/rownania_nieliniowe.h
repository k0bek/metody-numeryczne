#pragma once
/**
 * @file rownania_nieliniowe.h
 * @brief Rozwiazywanie rownan nieliniowych f(x) = 0
 *
 * Zawiera implementacje:
 *  - Metody bisekcji
 *  - Metody regula falsi (falszywej linii)
 *  - Metody Newtona-Raphsona (z analityczna lub numeryczna pochodna)
 *  - Metody siecznych
 *
 * @example
 *   auto f  = [](double x){ return x*x - 2.0; };
 *   auto df = [](double x){ return 2.0*x; };
 *   double pierwiastek  = bisekcja(f, 1.0, 2.0, 1e-8);   // ~1.41421356
 *   double pierwiastek2 = newton(f, df, 1.5, 1e-8);      // ~1.41421356
 */

#include <functional>

namespace biblioteka_mn {

/**
 * @brief Metoda bisekcji.
 *
 * Wymaga f(a)*f(b) < 0 (zmiana znaku na przedziale [a,b]).
 *
 * @param f          Funkcja ciagla
 * @param a          Lewy koniec przedzialu
 * @param b          Prawy koniec przedzialu
 * @param eps        Dokladnosc (kryterium stopu: |f(c)| < eps)
 * @param max_iter   Limit iteracji (domyslnie 1000)
 * @return           Przyblizone miejsce zerowe
 * @throws std::invalid_argument gdy f(a)*f(b) > 0
 */
double bisekcja(std::function<double(double)> f,
                double a, double b, double eps = 1e-6,
                int max_iter = 1000);

/**
 * @brief Metoda regula falsi (falszywa linia).
 *
 * Podobna do bisekcji, ale punkt podzialu wyznaczany jest jako przeciecie
 * prostej laczacej (a, f(a)) i (b, f(b)) z osia OX.
 *
 * @param f          Funkcja ciagla
 * @param a          Lewy koniec przedzialu (f(a)*f(b) < 0)
 * @param b          Prawy koniec przedzialu
 * @param eps        Dokladnosc
 * @param max_iter   Limit iteracji
 * @return           Przyblizone miejsce zerowe
 * @throws std::invalid_argument gdy f(a)*f(b) > 0
 */
double regula_falsi(std::function<double(double)> f,
                    double a, double b, double eps = 1e-6,
                    int max_iter = 1000);

/**
 * @brief Metoda Newtona-Raphsona z analityczna pochodna.
 *
 * @param f          Funkcja
 * @param df         Pochodna funkcji f
 * @param x0         Punkt startowy
 * @param eps        Dokladnosc (|x_nowe - x_stare| < eps)
 * @param max_iter   Limit iteracji
 * @return           Przyblizone miejsce zerowe
 */
double newton(std::function<double(double)> f,
              std::function<double(double)> df,
              double x0, double eps = 1e-6,
              int max_iter = 100);

/**
 * @brief Metoda Newtona-Raphsona z numeryczna pochodna (roznice centralne).
 *
 * Przydatna gdy pochodna analityczna jest trudna do wyznaczenia.
 *
 * @param f          Funkcja
 * @param x0         Punkt startowy
 * @param eps        Dokladnosc
 * @param max_iter   Limit iteracji
 * @return           Przyblizone miejsce zerowe
 */
double newton_numeryczny(std::function<double(double)> f,
                         double x0, double eps = 1e-6,
                         int max_iter = 100);

/**
 * @brief Metoda siecznych.
 *
 * Nie wymaga obliczania pochodnej; uzywa dwoch poprzednich przyblizen.
 *
 * @param f          Funkcja
 * @param x0         Pierwsze przyblizenie startowe
 * @param x1         Drugie przyblizenie startowe
 * @param eps        Dokladnosc
 * @param max_iter   Limit iteracji
 * @return           Przyblizone miejsce zerowe
 */
double sieczne(std::function<double(double)> f,
               double x0, double x1, double eps = 1e-6,
               int max_iter = 100);

} // namespace biblioteka_mn
