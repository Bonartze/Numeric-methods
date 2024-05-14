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
        return 1 / (x * x + 4);
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

    double RungeRombergRichardson(pair<double, double> results, int n) {
        double I_h1 = results.first;
        double I_h2 = results.second;
        return I_h2 + (I_h2 - I_h1) / (pow(2, n) - 1);
    }
};

int main() {
    IntegralCalculation it;

    auto squareResults = it.SquareMethod();
    double squareRRR = it.RungeRombergRichardson(squareResults, 2);

    auto trapResults = it.TrapMethod();
    double trapRRR = it.RungeRombergRichardson(trapResults, 2);

    auto simpsonResults = it.SimpsonMethod();
    double simpsonRRR = it.RungeRombergRichardson(simpsonResults, 4);

    cout << "Square Method: " << squareResults.first << ", " << squareResults.second << " | RRR: " << squareRRR << endl;
    cout << "Trap Method: " << trapResults.first << ", " << trapResults.second << " | RRR: " << trapRRR << endl;
    cout << "Simpson Method: " << simpsonResults.first << ", " << simpsonResults.second << " | RRR: " << simpsonRRR
         << endl;

    return 0;
}
