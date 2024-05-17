#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>
#include "desolve.h"
/* решить задачу Коши для ОДУ 2-го порядка
 * на указанном отрезке.*/
using namespace std;

double g(double x, double y, double z);

double f(double x, double y, double z);

double exact_solution(double x);

void print_data_with_error(const vector<tuple<double, double, double>> &data);

double g(double x, double y, double z) {
    return 4 * x * z - (4 * x * x - 3) * y + exp(x * x);
}

double f(double x, double y, double z) {
    (void) x;
    (void) y;
    return z;
}

double exact_solution(double x) {
    return (exp(x) + exp(-x) - 1) * exp(x * x);
}

using tddd = tuple<double, double, double>;

void print_data_with_error(const vector<tddd> &data) {
    cout << "x        y        Exact    Error" << endl;
    for (const auto &[x, y, z]: data) {
        double exact_y = exact_solution(x);
        double error = abs(exact_y - y);
        cout << x << " " << y << " " << exact_y << " " << error << endl;
    }
}

/*оду своддистя к оду первого проядка с помощью замены*/
int main() {
    cout.precision(6);
    cout << fixed;
    double l, r, y0, z0, h;
    cin >> l >> r;
    cin >> y0 >> z0 >> h;

    euler de_euler(l, r, f, g, y0, z0);
    vector<tddd> sol_euler = de_euler.solve(h);
    cout << "Метод Эйлера:" << endl;
    print_data_with_error(sol_euler);
    cout << "Погрешность вычислений:" << endl;
    double euler_err = runge_romberg(de_euler.solve(h), de_euler.solve(h / 2), 1);
    cout << euler_err << endl;

    runge de_runge(l, r, f, g, y0, z0);
    vector<tddd> sol_runge = de_runge.solve(h);
    cout << "Метод Рунге-Кутты:" << endl;
    print_data_with_error(sol_runge);
    cout << "Погрешность вычислений:" << endl;
    double runge_err = runge_romberg(de_runge.solve(h), de_runge.solve(h / 2), 4);
    cout << runge_err << endl;

    adams de_adams(l, r, f, g, y0, z0);
    vector<tddd> sol_adams = de_adams.solve(h);
    cout << "Метод Адамса:" << endl;
    print_data_with_error(sol_adams);
    cout << "Погрешность вычислений:" << endl;
    double adams_err = runge_romberg(de_adams.solve(h), de_adams.solve(h / 2), 4);
    cout << adams_err << endl;

    return 0;
}
