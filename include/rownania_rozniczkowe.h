#pragma once
/**
 * @file rownania_rozniczkowe.h
 * @brief Numeryczne rozwiazywanie rownan rozniczkowych zwyczajnych (RRZ)
 *
 * Zawiera implementacje metod jednokrokowych dla zagadnienia poczatkowego
 * y' = f(t, y),  y(t0) = y0:
 *  - Metoda Eulera
 *  - Metoda Heuna (Euler ulepszony / RK2)
 *  - Metoda punktu srodkowego (RK2 midpoint)
 *  - Klasyczna metoda Runge-Kutty 4. rzedu (RK4)
 *
 * Kazda funkcja zwraca wektor par {t, y} reprezentujacy trajektorie.
 *
 * @example
 *   // y' = -y,  y(0) = 1
 *   auto f = [](double t, double y){ return -y; };
 *   auto traj = rk4(f, 0.0, 1.0, 1.0, 100);
 *   // traj.back().second ~ exp(-1) ~ 0.3679
 */

#include <vector>
#include <functional>
#include <utility>

namespace biblioteka_mn {

/// Typ trajektorii: ciag par (t_i, y_i)
using Trajektoria = std::vector<std::pair<double, double>>;

/**
 * @brief Metoda Eulera (rzad 1).
 *
 * @param f      Prawa strona RRZ: f(t, y)
 * @param t0     Punkt startowy
 * @param y0     Warunek poczatkowy y(t0)
 * @param t_koniec  Koniec przedzialu calkowania
 * @param n      Liczba krokow
 * @return       Trajektoria {(t0,y0), (t1,y1), ..., (t_koniec, y_n)}
 */
Trajektoria euler(std::function<double(double, double)> f,
                  double t0, double y0, double t_koniec, int n);

/**
 * @brief Metoda Heuna (RK2 - korektor trapezowy, rzad 2).
 *
 * @param f      Prawa strona RRZ
 * @param t0     Punkt startowy
 * @param y0     Warunek poczatkowy
 * @param t_koniec  Koniec przedzialu
 * @param n      Liczba krokow
 * @return       Trajektoria
 */
Trajektoria heun(std::function<double(double, double)> f,
                 double t0, double y0, double t_koniec, int n);

/**
 * @brief Metoda punktu srodkowego (RK2 midpoint, rzad 2).
 *
 * @param f      Prawa strona RRZ
 * @param t0     Punkt startowy
 * @param y0     Warunek poczatkowy
 * @param t_koniec  Koniec przedzialu
 * @param n      Liczba krokow
 * @return       Trajektoria
 */
Trajektoria punkt_srodkowy(std::function<double(double, double)> f,
                           double t0, double y0, double t_koniec, int n);

/**
 * @brief Klasyczna metoda Runge-Kutty 4. rzedu (RK4).
 *
 * @param f      Prawa strona RRZ
 * @param t0     Punkt startowy
 * @param y0     Warunek poczatkowy
 * @param t_koniec  Koniec przedzialu
 * @param n      Liczba krokow
 * @return       Trajektoria
 */
Trajektoria rk4(std::function<double(double, double)> f,
                double t0, double y0, double t_koniec, int n);

} // namespace biblioteka_mn
