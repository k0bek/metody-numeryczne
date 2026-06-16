#include "uklady_rownan.h"
#include <stdexcept>
#include <cmath>
#include <vector>

namespace biblioteka_mn {

std::vector<double> eliminacja_gaussa(std::vector<std::vector<double>> A,
                                      std::vector<double> b) {
    int n = (int)A.size();

    // Budowanie macierzy rozszerzonej [A|b]
    for (int i = 0; i < n; ++i)
        A[i].push_back(b[i]);

    // Eliminacja z czesciowym wyborem elementu glownego
    for (int i = 0; i < n - 1; ++i) {
        // Pivoting - szukamy wiersza z najwiekszym elementem w kolumnie i
        int pivot = i;
        for (int k = i + 1; k < n; ++k)
            if (std::abs(A[k][i]) > std::abs(A[pivot][i]))
                pivot = k;
        std::swap(A[i], A[pivot]);

        if (std::abs(A[i][i]) < 1e-12)
            continue; // kolumna juz wyzerowana

        for (int k = i + 1; k < n; ++k) {
            double mnoznik = A[k][i] / A[i][i];
            for (int j = i; j <= n; ++j)
                A[k][j] -= mnoznik * A[i][j];
        }
    }

    // Sprawdzenie rozwiazywalnosci
    for (int i = 0; i < n; ++i) {
        bool wiersz_zerowy = true;
        for (int j = 0; j < n; ++j)
            if (std::abs(A[i][j]) > 1e-12) { wiersz_zerowy = false; break; }
        if (wiersz_zerowy && std::abs(A[i][n]) > 1e-12)
            throw std::runtime_error("eliminacja_gaussa: uklad sprzeczny");
    }

    // Podstawianie wstecz
    std::vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        double suma = A[i][n];
        for (int j = i + 1; j < n; ++j)
            suma -= A[i][j] * x[j];
        if (std::abs(A[i][i]) < 1e-12)
            throw std::runtime_error("eliminacja_gaussa: uklad nieoznaczony (pivot = 0)");
        x[i] = suma / A[i][i];
    }
    return x;
}

void rozklad_lu(const std::vector<std::vector<double>>& A,
                std::vector<std::vector<double>>& L,
                std::vector<std::vector<double>>& U) {
    int n = (int)A.size();
    L.assign(n, std::vector<double>(n, 0.0));
    U.assign(n, std::vector<double>(n, 0.0));

    for (int k = 0; k < n; ++k) {
        L[k][k] = 1.0;

        // Wiersz U
        for (int j = k; j < n; ++j) {
            double suma = 0.0;
            for (int s = 0; s < k; ++s)
                suma += L[k][s] * U[s][j];
            U[k][j] = A[k][j] - suma;
        }

        if (std::abs(U[k][k]) < 1e-12)
            throw std::runtime_error("rozklad_lu: pivot zerowy - rozklad niemozliwy");

        // Kolumna L
        for (int i = k + 1; i < n; ++i) {
            double suma = 0.0;
            for (int s = 0; s < k; ++s)
                suma += L[i][s] * U[s][k];
            L[i][k] = (A[i][k] - suma) / U[k][k];
        }
    }
}

std::vector<double> rozwiaz_lu(const std::vector<std::vector<double>>& A,
                               const std::vector<double>& b) {
    int n = (int)A.size();
    std::vector<std::vector<double>> L, U;
    rozklad_lu(A, L, U);

    // Podstawianie w przod: Lz = b
    std::vector<double> z(n);
    for (int i = 0; i < n; ++i) {
        double suma = 0.0;
        for (int j = 0; j < i; ++j)
            suma += L[i][j] * z[j];
        z[i] = b[i] - suma;
    }

    // Podstawianie wstecz: Ux = z
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        double suma = 0.0;
        for (int j = i + 1; j < n; ++j)
            suma += U[i][j] * x[j];
        x[i] = (z[i] - suma) / U[i][i];
    }
    return x;
}

} // namespace biblioteka_mn
