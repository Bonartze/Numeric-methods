#include <iostream>
#include <cmath>

using namespace std;

double epsilon = 0.000000001;

class SimpleIterationMethod {
private:
    double f(double x) {
        return log(x + 1) - 2 * x * x + 1;
    }

    double df(double x) {
        return (f(epsilon + x) - f(x)) / epsilon;
    }

    double g(double x) {
        return x - f(x) / df(x);
    }

public:
    double SimpleIterations() {
        double x0 = 0.8;
        double x_prev = x0;
        double x = x0;
        do {
            x_prev = x;
            x = g(x);
        } while (fabs(x - x_prev) >= epsilon);
        return x;
    }
};

class NewtonMethod {
private:
    double f(double x) {
        return log(x + 1) - 2 * x * x + 1;
    }

    double df(double x) {
        return (f(epsilon + x) - f(x)) / epsilon;
    }

public:
    double Newton_method() {
        double x0 = 0.8;
        double x_prev = x0;
        double x = x0;
        do {
            x_prev = x;
            x = x - f(x) / df(x);
        } while (fabs(x - x_prev) > epsilon);
        return x;
    }
};

int main() {
    SimpleIterationMethod si;
    NewtonMethod nm;
    cout << "Simple iterations: " << si.SimpleIterations() << endl;
    cout << "Newton method: " << nm.Newton_method() << endl;
    return 0;
}
