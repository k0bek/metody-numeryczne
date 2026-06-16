#pragma once
/**
 * @file uklady_rownan.h
 * @brief Rozwiazywanie ukladow rownan liniowych
 *
 * Zawiera implementacje:
 *  - Eliminacji Gaussa z czesciowym wyborem elementu glownego
 *  - Rozkladu LU
 *
 * @example
 *   vector<vector<double>> A = {{2,1,-1},{-3,-1,2},{-2,1,2}};
 *   vector<double> b = {8,-11,-3};
 *   auto x = eliminacja_gaussa(A, b);
 *   // x = {2, 3, -1}
 */

#include <vector>
#include <stdexcept>
#include <cmath>

namespace biblioteka_mn {

/**
 * @brief Eliminacja Gaussa z czesciowym wyborem elementu glownego.
 *
 * Rozwiazuje uklad Ax = b metoda eliminacji Gaussa z pivotingiem.
 *
 * @param A  Macierz wspolczynnikow (n x n)
 * @param b  Wektor prawej strony (n)
 * @return   Wektor rozwiazania x (n)
 * @throws   std::runtime_error gdy uklad jest sprzeczny lub nieoznaczony
 */
std::vector<double> eliminacja_gaussa(std::vector<std::vector<double>> A,
                                      std::vector<double> b);

/**
 * @brief Rozklad LU macierzy kwadratowej.
 *
 * Rozklada macierz A na iloczyn L*U, gdzie L - dolnotrojkatna (jedynki na
 * przekatnej), U - gornotrojkatna.
 *
 * @param A  Macierz wejsciowa (n x n)
 * @param L  Wynikowa macierz L (n x n)
 * @param U  Wynikowa macierz U (n x n)
 * @throws   std::runtime_error gdy pivot wynosi zero
 */
void rozklad_lu(const std::vector<std::vector<double>>& A,
                std::vector<std::vector<double>>& L,
                std::vector<std::vector<double>>& U);

/**
 * @brief Rozwiazuje Ax = b uzywajac rozkladu LU.
 *
 * Kolejno: rozklada A = LU, rozwiazuje Lz = b (podstawianie w przod),
 * nastepnie Ux = z (podstawianie wstecz).
 *
 * @param A  Macierz wspolczynnikow (n x n)
 * @param b  Wektor prawej strony (n)
 * @return   Wektor rozwiazania x (n)
 */
std::vector<double> rozwiaz_lu(const std::vector<std::vector<double>>& A,
                               const std::vector<double>& b);

} // namespace biblioteka_mn
