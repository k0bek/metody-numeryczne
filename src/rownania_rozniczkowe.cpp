#include "rownania_rozniczkowe.h"

namespace biblioteka_mn {

Trajektoria euler(std::function<double(double, double)> f,
                  double t0, double y0, double t_koniec, int n) {
    Trajektoria traj;
    traj.reserve(n + 1);
    double h = (t_koniec - t0) / n;
    double t = t0, y = y0;
    traj.push_back({t, y});
    for (int i = 0; i < n; ++i) {
        y = y + h * f(t, y);
        t = t + h;
        traj.push_back({t, y});
    }
    return traj;
}

Trajektoria heun(std::function<double(double, double)> f,
                 double t0, double y0, double t_koniec, int n) {
    Trajektoria traj;
    traj.reserve(n + 1);
    double h = (t_koniec - t0) / n;
    double t = t0, y = y0;
    traj.push_back({t, y});
    for (int i = 0; i < n; ++i) {
        double k1 = f(t, y);
        double y_tilde = y + h * k1;          // predyktor Eulera
        double k2 = f(t + h, y_tilde);
        y = y + (h / 2.0) * (k1 + k2);       // korektor trapezowy
        t = t + h;
        traj.push_back({t, y});
    }
    return traj;
}

Trajektoria punkt_srodkowy(std::function<double(double, double)> f,
                          double t0, double y0, double t_koniec, int n) {
    Trajektoria traj;
    traj.reserve(n + 1);
    double h = (t_koniec - t0) / n;
    double t = t0, y = y0;
    traj.push_back({t, y});
    for (int i = 0; i < n; ++i) {
        double k1 = f(t, y);
        double y_mid = y + (h / 2.0) * k1;
        y = y + h * f(t + h / 2.0, y_mid);
        t = t + h;
        traj.push_back({t, y});
    }
    return traj;
}

Trajektoria rk4(std::function<double(double, double)> f,
               double t0, double y0, double t_koniec, int n) {
    Trajektoria traj;
    traj.reserve(n + 1);
    double h = (t_koniec - t0) / n;
    double t = t0, y = y0;
    traj.push_back({t, y});
    for (int i = 0; i < n; ++i) {
        double k1 = h * f(t,           y);
        double k2 = h * f(t + h/2.0,   y + k1/2.0);
        double k3 = h * f(t + h/2.0,   y + k2/2.0);
        double k4 = h * f(t + h,        y + k3);
        y = y + (1.0/6.0) * (k1 + 2*k2 + 2*k3 + k4);
        t = t + h;
        traj.push_back({t, y});
    }
    return traj;
}

} // namespace biblioteka_mn
