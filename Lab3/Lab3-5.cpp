#include <iostream>
#include <cmath>
using namespace std;

class IntegralCalculation {
private:
    double x0 = -2;
    double xk = 2;
    double h1 = 1;
    double h2 = 0.5;

    double f(double x) {
        return 1/(x*x+4);
    }

public:
    pair<double, double> SquareMethod() {
        double res1 = 0;
        for (double x_s = x0 + h1; x_s <= xk; x_s += h1)
            res1 += h1 * f((x_s - h1 + x_s) / 2);
        double res2 = 0;
        for (double x_s = x0 + h2; x_s <= xk; x_s += h2)
            res2 += h2 * f((x_s - h2 + x_s) / 2);
        return {res1, res2};
    }

    pair<double, double> TrapMethod() {
        double res1 = h1 * f(x0) / 2;
        for (double x_s = x0 + h1; x_s < xk; x_s += h1)
            res1 += h1 * f(x_s);
        double res2 = h2 * f(x0) / 2;
        for (double x_s = x0 + h2; x_s < xk; x_s += h2)
            res2 += h2 * f(x_s);
        res1 += h1 * f(xk) / 2;
        res2 += h2 * f(xk) / 2;
        return {res1, res2};
    }

    pair<double, double> SimpsonMethod() {
        double res1 = (h1 * f(x0) / 3) + (h1 * f(xk) / 3);
        int k = 4;
        for (double x_s = x0 + h1; x_s < xk; x_s += h1) {
            res1 += k * h1 * f(x_s) / 3;
            k = (k == 4) ? 2 : 4;
        }
        k = 4;
        double res2 = h2 * f(x0) / 3 + h2 * f(xk) / 3;
        for (double x_s = x0 + h2; x_s < xk; x_s += h2) {
            res2 += k * h2 * f(x_s) / 3;
            k = (k == 4) ? 2 : 4;
        }
        return {res1, res2};
    }

    double RungeRombergRichardson() {
        auto [I_h1, I_h2] = SimpsonMethod();
        return I_h1 + (I_h1 - I_h2) / (0.17 - 1);
    }

};

int main() {
    IntegralCalculation it;
    cout << it.RungeRombergRichardson();
}