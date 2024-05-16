#include <cmath>
#include <iostream>
#include <vector>
#include <utility>

using namespace std;

class du_methods {
private:
    double g(double x) {
        return x;
    }

    double f(double x, double y, double z) {
        return 4 * x * z - (4 * x * x - 3) * y + exp(x * x);
    }

    double exact_solution(double x) {
        return (exp(x) + exp(-x) - 1) * exp(x * x);
    }

    double compute_error(double numerical, double exact) {
        return abs(numerical - exact);
    }

    double runge_romberg_richardson(double yh, double y2h, int p) {
        return yh + (yh - y2h) / (pow(2, p) - 1);
    }

public:
    void EulerMethod() {
        double x0 = 0;
        double h = 0.1;
        double x_last = 1;
        double x = x0;
        double z_s = 0, y = 1;

        while (x < x_last) {
            double last_y = y;
            y += h * g(z_s);
            z_s += h * f(x, last_y, g(z_s));
            x += h;
            double exact_y = exact_solution(x);
            double error = compute_error(y, exact_y);
            if (x <= x_last)
                cout << "x: " << x << ", y: " << y << ", y\' " << z_s << ", exact y: " << exact_y << ", error: "
                     << error << endl;
        }
    }

    pair<vector<double>, vector<double>> RungeKutta(double h) {
        double x0 = 0;
        double x_last = 1;
        double x = x0;
        double z = 0, y = 1;
        vector<double> adams_y, adams_z;
        adams_y.push_back(y);
        adams_z.push_back(z);
        int i = 0;
        while (x < x_last) {
            double K1_y = h * g(z);
            double L1_z = h * f(x, y, z);

            double K2_y = h * g(z + 0.5 * L1_z);
            double L2_z = h * f(x + 0.5 * h, y + 0.5 * K1_y, z + 0.5 * L1_z);

            double K3_y = h * g(z + 0.5 * L2_z);
            double L3_z = h * f(x + 0.5 * h, y + 0.5 * K2_y, z + 0.5 * L2_z);

            double K4_y = h * g(z + L3_z);
            double L4_z = h * f(x + h, y + K3_y, z + L3_z);

            y = y + (K1_y + 2 * K2_y + 2 * K3_y + K4_y) / 6;
            z = z + (L1_z + 2 * L2_z + 2 * L3_z + L4_z) / 6;
            if (i < 3) {
                adams_y.push_back(y);
                adams_z.push_back(z);
                i++;
            }
            x += h;
            double exact_y = exact_solution(x);
            double error = compute_error(y, exact_y);
            if (x <= x_last)

                cout << "x: " << x << ", y: " << y << ", y': " << z << ", exact y: " << exact_y << ", error: " << error
                     << endl;
        }
        cout << "\n\n";
        return {adams_y, adams_z};
    }

    void Adams() {
        double h = 0.1;
        double x0 = 0;
        double x_last = 1;
        auto [y, z] = RungeKutta(h);
        vector<double> adams_y(y.begin(), y.end());
        vector<double> adams_z(z.begin(), z.end());
        double x = x0 + 4 * h;
        while (x < x_last) {
            double f_n = f(x, adams_y[3], adams_z[3]);
            double f_n1 = f(x - h, adams_y[2], adams_z[2]);
            double f_n2 = f(x - 2 * h, adams_y[1], adams_z[1]);
            double f_n3 = f(x - 3 * h, adams_y[0], adams_z[0]);

            double g_n = g(adams_z[3]);
            double g_n1 = g(adams_z[2]);
            double g_n2 = g(adams_z[1]);
            double g_n3 = g(adams_z[0]);

            double y_next = adams_y[3] + (h / 24) * (55 * g_n - 59 * g_n1 + 37 * g_n2 - 9 * g_n3);
            double z_next = adams_z[3] + (h / 24) * (55 * f_n - 59 * f_n1 + 37 * f_n2 - 9 * f_n3);

            adams_y.erase(adams_y.begin());
            adams_z.erase(adams_z.begin());

            adams_y.push_back(y_next);
            adams_z.push_back(z_next);

            double exact_y = exact_solution(x);
            double error = compute_error(y_next, exact_y);
            x += h;
            if (x <= x_last)
                cout << "x: " << x << ", y: " << y_next << ", y': " << z_next << ", exact y: " << exact_y << ", error: "
                     << error << endl;
        }
    }

    void RungeRombergRichardsonComparison() {
        double h = 0.1;
        double h2 = 0.2; // This should be twice the step size h

        auto [yh_vec, _] = RungeKutta(h);
        double yh = yh_vec.back();

        auto [y2h_vec, _2] = RungeKutta(h2);
        double y2h = y2h_vec.back();

        double exact_y = exact_solution(1);

        double rrr_approx = runge_romberg_richardson(yh, y2h, 4);

        double rrr_error = compute_error(rrr_approx, exact_y);

        cout << "Runge-Romberg-Richardson approximation: " << rrr_approx << ", exact y: " << exact_y << ", RRR error: "
             << rrr_error << endl;
    }
};

int main() {
    du_methods du;

    cout << "Euler Method:\n";
    du.EulerMethod();

    cout << "\nRunge-Kutta Method:\n";
    du.RungeKutta(0.1);

    cout << "\nAdams Method:\n";
    du.Adams();

    cout << "\nRunge-Romberg-Richardson Comparison:\n";
    du.RungeRombergRichardsonComparison();
}
