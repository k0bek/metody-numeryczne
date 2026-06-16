# biblioteka_mn – Biblioteka Metod Numerycznych w C++

Biblioteka zaimplementowana w ramach kursu **Metody Numeryczne**.
Zawiera szesc kategorii algorytmow numerycznych w jezyku C++17.
Wszystkie nazwy plikow, funkcji i zmiennych sa w jezyku polskim.

---

## Struktura projektu

```
biblioteka_mn/
├── include/          # Pliki naglowkowe z dokumentacja API
│   ├── uklady_rownan.h
│   ├── interpolacja.h
│   ├── aproksymacja.h
│   ├── calkowanie.h
│   ├── rownania_rozniczkowe.h
│   └── rownania_nieliniowe.h
├── src/              # Implementacje
│   ├── uklady_rownan.cpp
│   ├── interpolacja.cpp
│   ├── aproksymacja.cpp
│   ├── calkowanie.cpp
│   ├── rownania_rozniczkowe.cpp
│   └── rownania_nieliniowe.cpp
├── testy/
│   └── testy.cpp     # Testy jednostkowe (≥2 testy na funkcje)
├── przyklady/
│   ├── przyklad1_stygniecie.cpp           # Stygniecie ciala – RRZ
│   ├── przyklad2_pierwiastki_calki.cpp    # Miejsca zerowe i calkowanie
│   └── przyklad3_pokaz.cpp                # Pelny pokaz wszystkich modulow
└── README.md
```

---

## Zakres funkcjonalny

### 1. Rozwiazywanie ukladow rownan liniowych (`uklady_rownan.h`)
| Funkcja | Opis |
|---|---|
| `eliminacja_gaussa(A, b)` | Eliminacja Gaussa z czesciowym pivotingiem |
| `rozklad_lu(A, L, U)` | Rozklad LU macierzy |
| `rozwiaz_lu(A, b)` | Rozwiazanie przez rozklad LU |

### 2. Interpolacja (`interpolacja.h`)
| Funkcja | Opis |
|---|---|
| `lagrange(x, wezly_x, wezly_y)` | Wielomian interpolacyjny Lagrange'a |
| `wspolczynniki_newtona(wezly_x, wezly_y)` | Ilorazy roznicowe (wspolczynniki Newtona) |
| `oblicz_newtona(x, wezly_x, wspolczynniki)` | Ewaluacja wielomianu Newtona |

### 3. Aproksymacja (`aproksymacja.h`)
| Funkcja | Opis |
|---|---|
| `aproksymacja_mnk(f, a, b, stopien)` | Aproksymacja MNK wielomianem stopnia `stopien` |
| `horner(wspolczynniki, x)` | Schemat Hornera – ewaluacja wielomianu |

### 4. Calkowanie numeryczne (`calkowanie.h`)
| Funkcja | Opis |
|---|---|
| `metoda_trapezow(f, a, b, n)` | Metoda trapezow |
| `metoda_simpsona(f, a, b, n)` | Regula Simpsona 1/3 |
| `gauss_legendre(f, a, b, liczba_wezlow, podprzedzialy)` | Kwadratura Gaussa-Legendre'a (2/3/4 wezly) |

### 5. Rownania rozniczkowe (`rownania_rozniczkowe.h`)
| Funkcja | Opis |
|---|---|
| `euler(f, t0, y0, t_koniec, n)` | Metoda Eulera (rzad 1) |
| `heun(f, t0, y0, t_koniec, n)` | Metoda Heuna / RK2 |
| `punkt_srodkowy(f, t0, y0, t_koniec, n)` | Metoda punktu srodkowego |
| `rk4(f, t0, y0, t_koniec, n)` | Klasyczna metoda RK4 |

### 6. Rownania nieliniowe (`rownania_nieliniowe.h`)
| Funkcja | Opis |
|---|---|
| `bisekcja(f, a, b, eps)` | Metoda bisekcji |
| `regula_falsi(f, a, b, eps)` | Metoda falszywej linii |
| `newton(f, df, x0, eps)` | Metoda Newtona (analityczna pochodna) |
| `newton_numeryczny(f, x0, eps)` | Metoda Newtona (numeryczna pochodna) |
| `sieczne(f, x0, x1, eps)` | Metoda siecznych |

## Instalacja i kompilacja

Wymagany kompilator obslugujacy **C++17** (np. `g++ >= 7`).

### Budowanie wszystkiego (testy + przyklady)

```bash
make wszystko
```

### Tylko testy jednostkowe

```bash
make testy
```

### Tylko przyklady

```bash
make przyklad1
make przyklad2
make przyklad3
```


## Testy
Plik `testy/testy.cpp` zawiera **≥ 2 testy na kazda funkcje**, obejmujace:
- przypadki poprawne (weryfikacja wartosci wzgledem rozwiazania analitycznego),
- przypadki bledne (sprawdzenie rzucania wyjatkow).

Po uruchomieniu `make testy` na ekranie pojawi sie podsumowanie `OK / BLAD`.
