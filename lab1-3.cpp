#include <iostream>
#include <vector>
#include <cmath>

double epsilon = 0.00000000001;

using namespace std;

class IterationsAlgorithms {
private:
    vector<vector<double>> matrix = {{-7, -1,  2,  2},
                                     {3,  -20, 0,  -8},
                                     {-9, 1,   18, -6},
                                     {-1, 0,   -1, -6}};
    vector<vector<double>> results = {{-24},
                                      {-47},
                                      {28},
                                      {-50}};
public:
    vector<vector<double>> simple_iterations() {
        vector<vector<double>> x0 = {{1},
                                     {0},
                                     {1},
                                     {2}};
        vector<vector<double>> b(matrix.size(), (vector<double>(matrix.size(), 0)));
        vector<vector<double>> beta(matrix.size(), (vector<double>(1, 0)));
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix.size(); j++) {
                if (i != j)
                    b[i][j] = -matrix[i][j] / matrix[i][i];
            }
        }
        for (size_t i = 0; i < matrix.size(); i++)
            beta[i][0] = results[i][0] / matrix[i][i];
        auto x_prev = x0;
        auto x_curr = x0;
        bool flag = true;
        int n = 0;
        while (flag) {
            n++;
            for (size_t i = 0; i < matrix.size(); i++) {
                x_curr[i][0] = beta[i][0];
                for (size_t j = 0; j < matrix.size(); j++) {
                    x_curr[i][0] += b[i][j] * x_prev[j][0];
                }
                if (fabs(x_curr[i][0] - x_prev[i][0]) < epsilon)
                    flag = false;
                x_prev = x_curr;
            }
        }
        cout << n << endl;
        return x_curr;
    }

    vector<vector<double>> Zaidel_method() {
        vector<vector<double>> x0 = {{1},
                                     {0},
                                     {1},
                                     {2}};
        vector<vector<double>> b(matrix.size(), (vector<double>(matrix.size(), 0)));
        vector<vector<double>> beta(matrix.size(), (vector<double>(1, 0)));
        for (size_t i = 0; i < matrix.size(); i++) {
            for (size_t j = 0; j < matrix.size(); j++) {
                if (i != j)
                    b[i][j] = -matrix[i][j] / matrix[i][i];
            }
        }
        for (size_t i = 0; i < matrix.size(); i++)
            beta[i][0] = results[i][0] / matrix[i][i];
        auto x_prev = x0;
        auto x_curr = x0;
        bool flag = true;
        int n = 0;
        while (flag) {
            n++;
            for (size_t i = 0; i < matrix.size(); i++) {
                x_curr[i][0] = beta[i][0];
                for (size_t j = 0; j < matrix.size(); j++) {
                    if (j >= i)
                        x_curr[i][0] += b[i][j] * x_prev[j][0];
                    else
                        x_curr[i][0] += b[i][j] * x_curr[j][0];

                }
                if (fabs(x_curr[i][0] - x_prev[i][0]) < epsilon)
                    flag = false;
                x_prev = x_curr;
            }
        }
        cout << n << endl;
        return x_curr;
    }
};

int main() {
    IterationsAlgorithms a;
    cin >> epsilon;
    auto b = a.Zaidel_method();
    for (auto &c: b) {
        for (auto &z: c)
            cout << z << ' ';
    }
}
