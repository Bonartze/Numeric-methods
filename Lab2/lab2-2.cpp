#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double epsilon = 0.000001;

class NewtonMethod {
private:
    vector<double> F(vector<double> x) {
        return {x[0] * x[0] + x[1] * x[1] - 9, x[0] - exp(x[1]) + 3};
    }


    double get_norm(vector<double> vec) {
        double sum = 0.0;
        for (auto &a: vec)
            sum += a * a;
        return sqrt(sum);
    }

    vector<double> diff(vector<double> v1, vector<double> v2) {
        vector<double> res;
        for (size_t i = 0; i < v1.size(); i++) {
            res.push_back(v1[i] - v2[i]);
        }
        return res;
    }

    void get_inverse_2(vector<vector<double>> &matrix) {
        double deter = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        double temp = matrix[0][0];
        matrix[0][0] = matrix[1][1];
        matrix[1][1] = temp;

        matrix[0][1] *= -1;
        matrix[1][0] *= -1;

        for (size_t i = 0; i < matrix.size(); i++)
            for (size_t j = 0; j < matrix[0].size(); j++)
                matrix[i][j] /= deter;
    }


    vector<vector<double>> get_jacobi_matrix(vector<double> x) {
        return {{2 * x[0], 2 * x[1]},
                {1,        -exp(x[1])}};
    }

    vector<double> mult(vector<vector<double>> jacobian, vector<double> x) {
        vector<double> new_x(x.size());
        for (size_t i = 0; i < jacobian.size(); i++)
            for (size_t j = 0; j < jacobian.size(); j++) {
                new_x[i] += jacobian[i][j] * x[j];
            }
        return new_x;

    }

public:
    vector<double> Newton_method() {
        vector<double> x = {1.7, 2.2};
        auto x_prev = x;
        while (true) {
            x_prev = x;
            auto jacobian = get_jacobi_matrix(x);
            get_inverse_2(jacobian);
            x = diff(x, mult(jacobian, F(x)));
            if (get_norm(x) - get_norm(x_prev) < epsilon)
                break;
        }
        return x;
    }
};

class SimpleIterations {
private:
    vector<double> g(vector<double> x) {
        double x1 = sqrt(9 - x[1] * x[1]);
        double x2 = log(x[0] + 3);
        return {x1, x2};
    }


    double get_norm(vector<double> vec) {
        double sum = 0.0;
        for (auto &a: vec)
            sum += a * a;
        return sqrt(sum);
    }

    vector<double> diff(vector<double> v1, vector<double> v2) {
        vector<double> res;
        for (size_t i = 0; i < v1.size(); i++)
            res.push_back(v1[i] - v2[i]);
        return res;
    }

public:
    vector<double> Simple_Iterations() {
        vector<double> x = {2.4, 1.7};
        auto x_prev = x;
        while (true) {
            x = g(x);
            if (get_norm(x) - get_norm(x_prev) < epsilon)
                break;
            x_prev = x;
        }
        return x;
    }
};

int main() {
    SimpleIterations nm;
    cout << "Simple iterations: ";
    for (auto n: nm.Simple_Iterations())
        cout << n << ' ';
    cout << endl;
    NewtonMethod m;
    cout << "Newton method: ";
    for (auto n: m.Newton_method())
        cout << n << ' ';
    return 0;
}
