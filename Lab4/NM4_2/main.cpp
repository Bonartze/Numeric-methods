#include "solver.h"

using namespace std;

double g(double x, double y, double z) {
    return (4 * y - 4 * x * z) / (2 * x + 1);
}

double f(double x, double y, double z) {
    (void) x;
    (void) y;
    return z;
}

double px(double x) {
    return (4 * x) / (2 * x + 1);
}

double qx(double x) {
    return -4 / (2 * x + 1);
}

double fx(double x) {
    (void) x;
    return 0.0;
}

using tddd = tuple<double, double, double>;

int main() {
    cout.precision(6);
    cout << fixed;
    double h, eps;
    cin >> h >> eps;
    /*
    Краевые условия 3 рода
    alpha * y(a) + beta * y'(a) = y0
    delta * y(b) + gamma * y'(b) = y1
    */
    double a = -2, b = -1;
    double alpha = 2, beta = 1, y0 = -9;
    double delta = 0, gamma = 1, y1 = 1;

    shooting de_shooting(a, b, f, g, alpha, beta, y0, delta, gamma, y1);
    vector<tddd> sol_shooting = de_shooting.solve(h, eps);
    cout << "Метод стрельбы:" << endl;
    print_data(sol_shooting);
    cout << "Погрешность вычислений:" << endl;
    double shooting_err = runge_romberg(de_shooting.solve(h, eps), de_shooting.solve(h / 2, eps), 4);
    cout << shooting_err << endl;

    fin_dif de_fin_dif(a, b, px, qx, fx, alpha, beta, y0, delta, gamma, y1);
    vector<tddd> sol_fin_dif = de_fin_dif.solve(h);
    cout << "Конечно-разностный метод:" << endl;
    print_data(sol_fin_dif);
    cout << "Погрешность вычислений:" << endl;
    double fin_dif_err = runge_romberg(de_fin_dif.solve(h), de_fin_dif.solve(h / 2), 2);
    cout << fin_dif_err << endl;
}