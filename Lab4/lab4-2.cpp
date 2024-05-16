#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <tuple>

using namespace std;

double runge_romberg(double h1, double h2, double y1, double y2, int n = 2) {
    return abs((y1 - y2) / (pow((h2 / h1), n) - 1.0));
}

double exact(double x) {
    return 3 * x + exp(-2 * x);  // The exact solution from the table
}

double f_1(double x, double y, double z) {
    return z;
}

double g_1(double x, double y, double z) {
    return (4 * x * z - 4 * y) / (2 * x + 1);
}

double p(double x) {
    return 0.;
}

double q(double x) {
    return 0.;  // Adjust if necessary
}

double f(double x) {
    return 0.;
}

vector<double>
runge_kutt(double xa, double ya, double eta, double (*f)(double, double, double), double (*g)(double, double, double),
           double h, int n) {
    vector<double> x(n), y(n), z(n);
    y[0] = ya;
    z[0] = eta;

    for (int i = 0; i < n - 1; ++i) {
        double k1 = h * f(x[i], y[i], z[i]);
        double l1 = h * g(x[i], y[i], z[i]);
        double k2 = h * f(x[i] + h / 2, y[i] + k1 / 2, z[i] + l1 / 2);
        double l2 = h * g(x[i] + h / 2, y[i] + k1 / 2, z[i] + l1 / 2);
        double k3 = h * f(x[i] + h / 2, y[i] + k2 / 2, z[i] + l2 / 2);
        double l3 = h * g(x[i] + h / 2, y[i] + k2 / 2, z[i] + l2 / 2);
        double k4 = h * f(x[i] + h, y[i] + k3, z[i] + l3);
        double l4 = h * g(x[i] + h, y[i] + k3, z[i] + l3);

        y[i + 1] = y[i] + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        z[i + 1] = z[i] + (l1 + 2 * l2 + 2 * l3 + l4) / 6;
        x[i + 1] = x[i] + h;
    }
    return y;
}

pair<vector<double>, vector<double>>
shooting(double xa, double xb, double ya, double yb, double h, double (*f)(double, double, double),
         double (*g)(double, double, double)) {
    int n = ceil((xb - xa) / h);
    vector<double> eta = {1, 0.8};
    double eps = 0.000001;
    vector<double> F;

    for (double et: eta) {
        vector<double> y = runge_kutt(xa, ya, et, f, g, h, n);
        F.push_back(y.back() - yb);
    }

    int k = 2;
    while (true) {
        eta.push_back(eta[k - 1] - (eta[k - 1] - eta[k - 2]) / (F[k - 1] - F[k - 2]) * F[k - 1]);

        vector<double> y = runge_kutt(xa, ya, eta[k], f, g, h, n);
        F.push_back(y.back() - yb);

        if (abs(F[k]) < eps) {
            break;
        }
        ++k;
    }

    vector<double> x(n);
    for (int i = 0; i < n; ++i) {
        x[i] = xa + i * h;
    }
    return {x, runge_kutt(xa, ya, eta.back(), f, g, h, n)};
}

vector<double>
tridiagonal(const vector<double> &a, const vector<double> &b, const vector<double> &c, const vector<double> &d, int n) {
    vector<double> y(n);
    vector<double> P(n), Q(n);

    P[0] = -c[0] / b[0];
    Q[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double denom = b[i] + a[i] * P[i - 1];
        P[i] = -c[i] / denom;
        Q[i] = (d[i] - a[i] * Q[i - 1]) / denom;
    }

    y[n - 1] = Q[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        y[i] = P[i] * y[i + 1] + Q[i];
    }

    return y;
}

vector<double> finite_diff(const vector<double> &x, double ya, double yb, double h, int n) {
    vector<double> a(n), b(n), c(n), d(n);
    int last = n - 1;

    b[0] = -2 + pow(h, 2) * q(x[0]);
    c[0] = 1 + p(x[0]) * h / 2;
    d[0] = (pow(h, 2) * f(x[0])) - (1 - p(x[0]) * h / 2) * ya;

    for (int k = 1; k < last; ++k) {
        a[k] = 1 - p(x[k]) * h / 2;
        b[k] = -2 + pow(h, 2) * q(x[k]);
        c[k] = 1 + p(x[k]) * h / 2;
        d[k] = pow(h, 2) * f(x[k]);
    }

    a[last] = 1 - p(x[last]) * h / 2;
    b[last] = -2 + pow(h, 2) * q(x[last]);
    d[last] = (pow(h, 2) * f(x[last])) - (1 + p(x[last]) * h / 2) * yb;

    return tridiagonal(a, b, c, d, n);
}

int main() {
    double h = 0.1;
    double xa = -2, xb = 0;
    double ya = -4.5;
    double yb = 1;
    int n = ceil((xb - xa) / h);

    cout << "Finite Difference Method" << endl;
    vector<double> x(n);
    for (int i = 0; i < n; ++i) {
        x[i] = xa + i * h;
    }
    vector<double> y = finite_diff(x, ya, yb, h, n);
    cout << "x        ";
    for (double xi: x) {
        printf("%5.5f ", xi);
    }
    cout << endl;
    cout << "y        ";
    for (double yi: y) {
        printf("%5.5f ", yi);
    }
    cout << endl;
    cout << "Error    ";
    for (int i = 0; i < n; ++i) {
        double val = abs(exact(x[i]) - y[i]);
        printf("%5.5f ", val);
    }
    cout << endl;
    cout << "Runge Romberg Error ";
    vector<double> x2(2 * n);
    for (int i = 0; i < 2 * n; ++i) {
        x2[i] = xa + i * (h / 2);
    }
    vector<double> y2 = finite_diff(x2, ya, yb, h / 2, 2 * n);
    for (int i = 0; i < n; ++i) {
        printf("%5.5f ", runge_romberg(h, h / 2, y[i], y2[i * 2]));
    }
    cout << endl;

    vector<double> x_shooting, y_shooting;
    tie(x_shooting, y_shooting) = shooting(xa, xb, ya, yb, h, f_1, g_1);
    cout << endl;
    cout << "Shooting Method" << endl;
    cout << "x        ";
    for (double xi: x_shooting) {
        printf("%5.5f ", xi);
    }
    cout << endl;
    cout << "y        ";
    for (double yi: y_shooting) {
        printf("%5.5f ", yi);
    }
    cout << endl;
    cout << "Error    ";
    for (int i = 0; i < n; ++i) {
        double val = abs(exact(x[i]) - y_shooting[i]);
        printf("%5.5f ", val);
    }
    cout << endl;
    cout << "Runge Romberg Error ";
    vector<double> x2_shooting, y2_shooting;
    tie(x2_shooting, y2_shooting) = shooting(xa, xb, ya, yb, h / 2, f_1, g_1);
    for (int i = 0; i < n; ++i) {
        printf("%5.5f ", runge_romberg(h, h / 2, y_shooting[i], y2_shooting[i * 2]));
    }
    cout << endl;

    return 0;
}
