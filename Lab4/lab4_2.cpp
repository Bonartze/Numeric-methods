#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <map>

// Define lambda functions as std::function
using LambdaFunction = std::function<double(double)>;
using DDFunction = std::function<double(double, double, double)>;

// Function to solve differential equations using the Runge-Kutta method
std::pair<std::vector<double>, std::vector<double>> runge_kutta(DDFunction ddy, const std::vector<double>& borders, double y0, double z0, double h) {
    std::vector<double> x, y, z;
    for (double val = borders[0]; val <= borders[1] + h; val += h) {
        x.push_back(val);
    }
    int N = x.size();
    y.resize(N, 0);
    z.resize(N, 0);
    y[0] = y0;
    z[0] = z0;

    for (int i = 0; i < N - 1; ++i) {
        double K1 = h * z[i];
        double L1 = h * ddy(x[i], y[i], z[i]);
        double K2 = h * (z[i] + 0.5 * L1);
        double L2 = h * ddy(x[i] + 0.5 * h, y[i] + 0.5 * K1, z[i] + 0.5 * L1);
        double K3 = h * (z[i] + 0.5 * L2);
        double L3 = h * ddy(x[i] + 0.5 * h, y[i] + 0.5 * K2, z[i] + 0.5 * L2);
        double K4 = h * (z[i] + L3);
        double L4 = h * ddy(x[i] + h, y[i] + K3, z[i] + L3);
        double delta_y = (K1 + 2 * K2 + 2 * K3 + K4) / 6;
        double delta_z = (L1 + 2 * L2 + 2 * L3 + L4) / 6;
        y[i + 1] = y[i] + delta_y;
        z[i + 1] = z[i] + delta_z;
    }
    return { y, z };
}

double diff_left(const std::map<std::string, double>& bcondition, double h, LambdaFunction f, double x0) {
    return (bcondition.at("c") - (bcondition.at("b") / h) * f(x0)) / (bcondition.at("a") - (bcondition.at("b") / h));
}

double diff_right(const std::map<std::string, double>& bcondition, double h, const std::vector<double>& y) {
    return (bcondition.at("c") + (bcondition.at("b") / h) * y[y.size() - 2]) / (bcondition.at("a") + (bcondition.at("b") / h));
}

std::vector<double> solve_gauss(std::vector<std::vector<double>>& A, std::vector<double>& b) {
    int n = b.size();
    for (int i = 0; i < n; ++i) {
        // Find the maximum element in column i
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(A[k][i]) > std::abs(A[maxRow][i])) {
                maxRow = k;
            }
        }

        // Swap maximum row with current row
        std::swap(A[maxRow], A[i]);
        std::swap(b[maxRow], b[i]);

        // Make all rows below this one 0 in current column
        for (int k = i + 1; k < n; ++k) {
            double c = -A[k][i] / A[i][i];
            for (int j = i; j < n; ++j) {
                if (i == j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
            b[k] += c * b[i];
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i] / A[i][i];
        for (int k = i - 1; k >= 0; --k) {
            b[k] -= A[k][i] * x[i];
        }
    }
    return x;
}

std::vector<double> shooting_method(DDFunction ddy, const std::vector<double>& borders, const std::map<std::string, double>& bcondition1,
                                    const std::map<std::string, double>& bcondition2, double h, LambdaFunction f) {
    double y0 = diff_left(bcondition1, h, f, borders[0]);
    double eta1 = 0.5;
    double eta2 = 2.0;
    auto resolve1 = runge_kutta(ddy, borders, y0, eta1, h).first;
    auto resolve2 = runge_kutta(ddy, borders, y0, eta2, h).first;
    double Phi1 = resolve1.back() - diff_right(bcondition2, h, resolve1);
    double Phi2 = resolve2.back() - diff_right(bcondition2, h, resolve2);
    while (std::abs(Phi2 - Phi1) > h / 10) {
        double temp = eta2;
        eta2 = eta2 - (eta2 - eta1) / (Phi2 - Phi1) * Phi2;
        eta1 = temp;
        resolve1 = runge_kutta(ddy, borders, y0, eta1, h).first;
        resolve2 = runge_kutta(ddy, borders, y0, eta2, h).first;
        Phi1 = resolve1.back() - diff_right(bcondition2, h, resolve1);
        Phi2 = resolve2.back() - diff_right(bcondition2, h, resolve2);
    }
    return runge_kutta(ddy, borders, y0, eta2, h).first;
}

std::vector<double> finite_difference_method(DDFunction ddy, LambdaFunction f, const std::map<std::string, double>& bcondition1,
                                             const std::map<std::string, double>& bcondition2, const std::map<std::string, LambdaFunction>& equation, const std::vector<double>& borders, double h) {
    std::vector<double> x, b;
    for (double val = borders[0]; val <= borders[1] + h; val += h) {
        x.push_back(val);
    }
    int N = x.size();
    std::vector<std::vector<double>> A(N, std::vector<double>(N, 0));
    b.resize(N, 0);

    A[0][0] = bcondition1.at("a") - bcondition1.at("b") / h;
    A[0][1] = bcondition1.at("b") / h;
    b[0] = bcondition1.at("c");

    for (int i = 1; i < N - 1; ++i) {
        A[i][i - 1] = 1 / (h * h) - equation.at("p")(x[i]) / (2 * h);
        A[i][i] = -2 / (h * h) + equation.at("q")(x[i]);
        A[i][i + 1] = 1 / (h * h) + equation.at("p")(x[i]) / (2 * h);
        b[i] = equation.at("f")(x[i]);
    }

    A[N - 1][N - 2] = -bcondition2.at("b") / h;
    A[N - 1][N - 1] = bcondition2.at("a") + bcondition2.at("b") / h;
    b[N - 1] = bcondition2.at("c");

    return solve_gauss(A, b);
}

std::vector<double> runge_rombert(const std::vector<double>& y1, const std::vector<double>& y2, double h1, double h2, int p) {
    int k = h1 > h2 ? static_cast<int>(h1 / h2) : static_cast<int>(h2 / h1);
    int N = std::min(y1.size(), y2.size() * k);
    std::vector<double> y(N);

    if (h1 > h2) {
        for (int i = 0; i < N; ++i) {
            y[i] = y2[i * k] + (y2[i * k] - y1[i]) / (pow(k, p) - 1);
        }
    } else {
        for (int i = 0; i < N; ++i) {
            y[i] = y1[i * k] + (y1[i * k] - y2[i]) / (pow(k, p) - 1);
        }
    }

    return y;
}

double sqr_error(const std::vector<double>& y, const std::vector<double>& y_correct) {
    double error = 0;
    for (size_t i = 0; i < y.size(); ++i) {
        error += pow(y[i] - y_correct[i], 2);
    }
    return sqrt(error);
}

int main() {
    std::map<std::string, double> bcondition1, bcondition2;
    double h = 0.01;

    bcondition1 = {
            {"a", 0},
            {"b", 1},
            {"c", 1}
    };

    bcondition2 = {
            {"a", 1},
            {"b", 2},
            {"c", -9}
    };

    std::vector<double> borders = { -2, 0 };

    DDFunction ddf = [](double x, double y, double dy) { return (-4 * x * dy + 4 * y) / (2 * x + 1); };
    LambdaFunction f = [](double x) { return 3 * x + exp(-2 * x); };
    LambdaFunction p = [](double x) { return 4 * x / (2 * x + 1); };
    LambdaFunction q = [](double x) { return -4 / (2 * x + 1); };
    LambdaFunction right_f = [](double x) { return 0; };

    std::map<std::string, LambdaFunction> equation = { {"p", p}, {"q", q}, {"f", right_f} };

    std::vector<double> x, y_exact;
    for (double val = borders[0]; val <= borders[1] + h; val += h) {
        x.push_back(val);
    }
    for (double val : x) {
        y_exact.push_back(f(val));
    }

    std::vector<double> y_shooting = shooting_method(ddf, borders, bcondition1, bcondition2, h, f);
    std::vector<double> y_finite_diff = finite_difference_method(ddf, f, bcondition1, bcondition2, equation, borders, h);

    double h2 = h / 2;
    std::vector<double> y_shooting_2 = shooting_method(ddf, borders, bcondition1, bcondition2, h2, f);
    std::vector<double> y_finite_diff_2 = finite_difference_method(ddf, f, bcondition1, bcondition2, equation, borders, h2);

    std::cout << "Runge Rombert errors:\n";
    std::cout << "Shooting method: " << sqr_error(y_shooting, runge_rombert(y_shooting, y_shooting_2, h, h2, 1)) << "\n";
    std::cout << "Finite difference method: " << sqr_error(y_finite_diff, runge_rombert(y_finite_diff, y_finite_diff_2, h, h2, 1)) << "\n\n";

    std::cout << "Exact solution errors:\n";
    std::cout << "Shooting method: " << sqr_error(y_shooting_2, y_exact) << "\n";
    std::cout << "Finite difference method: " << sqr_error(y_finite_diff, y_exact) << "\n\n";

    // Print the solutions for comparison
    std::cout << "x\tExact\tShooting\tFinite Difference\n";
    for (size_t i = 0; i < x.size(); ++i) {
        std::cout << x[i] << "\t" << y_exact[i] << "\t" << y_shooting[i] << "\t" << y_finite_diff[i] << "\n";
    }

    std::cout << "\n\nx\tExact\tShooting\tFinite Difference\n";
    for (size_t i = 0; i < x.size(); ++i) {
        std::cout << x[i] << "\t" << y_exact[i] << "\t" << y_shooting_2[i] << "\t" << y_finite_diff[i] << "\n";
    }
    return 0;
}
