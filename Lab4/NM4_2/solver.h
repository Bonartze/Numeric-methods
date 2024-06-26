#ifndef NM4_2_SOLVER_H
#define NM4_2_SOLVER_H

#include <cmath>
#include <vector>
#include <tuple>
#include <functional>
#include "tridiag.cpp"
#include "desolve.h"

using vec = std::vector<double>;
using vect = std::vector<std::tuple<double, double, double>>;

class shooting {
private:
    double a, b;
    func f, g;
    double alpha, beta, y0;
    double delta, gamma, y1;

public:
    shooting(const double _a, const double _b,
             const func _f, const func _g,
             const double _alpha, const double _beta, const double _y0,
             const double _delta, const double _gamma, const double _y1)
            : a(_a), b(_b), f(_f), g(_g),
              alpha(_alpha), beta(_beta), y0(_y0),
              delta(_delta), gamma(_gamma), y1(_y1) {}

    double get_start_cond(double eta) {
        return (y0 - alpha * eta) / beta;
    }

    double get_eta_next(double eta_prev, double eta, const vect sol_prev, const vect sol) {
        double yb_prev = std::get<1>(sol_prev.back());
        double zb_prev = std::get<2>(sol_prev.back());
        double phi_prev = delta * yb_prev + gamma * zb_prev - y1;
        double yb = std::get<1>(sol.back());
        double zb = std::get<2>(sol.back());
        double phi = delta * yb + gamma * zb - y1;
        return eta - (eta - eta_prev) / (phi - phi_prev) * phi;
    }

    vect solve(double h, double eps) {
        double eta_prev = 1.0;
        double eta = 0.8;
        while (1) {
            double runge_z0_prev = get_start_cond(eta_prev);
            runge de_solver_prev(a, b, f, g, eta_prev, runge_z0_prev);
            vect sol_prev = de_solver_prev.solve(h);

            double runge_z0 = get_start_cond(eta);
            runge de_solver(a, b, f, g, eta, runge_z0);
            vect sol = de_solver.solve(h);

            double eta_next = get_eta_next(eta_prev, eta, sol_prev, sol);
            if (std::abs(eta_next - eta) < eps) {
                return sol;
            } else {
                eta_prev = eta;
                eta = eta_next;
            }
        }
    }
};

class fin_dif {
private:
    using fx = std::function<double(double)>;
    using tridiag = Tridiag<double>;

    double a, b;
    fx p, q, f;
    double alpha, beta, y0;
    double delta, gamma, y1;

public:
    fin_dif(const double _a, const double _b,
            const fx _p, const fx _q, const fx _f,
            const double _alpha, const double _beta, const double _y0,
            const double _delta, const double _gamma, const double _y1)
            : a(_a), b(_b), p(_p), q(_q), f(_f),
              alpha(_alpha), beta(_beta), y0(_y0),
              delta(_delta), gamma(_gamma), y1(_y1) {}

    vect solve(double h) {
        size_t n = static_cast<size_t>((b - a) / h);
        vec xk(n + 1);
        for (size_t i = 0; i <= n; ++i) {
            xk[i] = a + h * i;
        }
        vec A(n + 1), B(n + 1), C(n + 1), D(n + 1);
        B[0] = alpha - beta / h;
        C[0] = beta / h;
        D[0] = y0;
        A[n] = -gamma / h;
        B[n] = delta + gamma / h;
        D[n] = y1;
        for (size_t i = 1; i < n; ++i) {
            A[i] = 1.0 - p(xk[i]) * h * 0.5;
            B[i] = -2.0 + h * h * q(xk[i]);
            C[i] = 1.0 + p(xk[i]) * h * 0.5;
            D[i] = h * h * f(xk[i]);
        }
        tridiag sys_eq(A, B, C);
        vec yk = sys_eq.Solve(D);
        vect res;
        for (size_t i = 0; i <= n; ++i) {
            res.push_back(std::make_tuple(xk[i], yk[i], NAN));
        }
        return res;
    }
};

#endif // NM4_2_SOLVER_H
