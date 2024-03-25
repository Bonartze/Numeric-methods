#include <iostream>
#include <vector>

using namespace std;

class RunMethod {
private:
    vector<vector<double>> matrix = {{-11, -8,  0,  0,   0},
                                     {9,   -17, 1,  0,   0},
                                     {0,   -4,  20, 9,   0},
                                     {0,   0,   -4, -14, 3},
                                     {0,   0,   0,  -6,  14}};
    vector<vector<double>> results = {{99},
                                      {-75},
                                      {66},
                                      {54},
                                      {8}};
public:
    vector<double> run_through() {
        // for first row of matrix
        vector<double> y(matrix.size());
        vector<double> a(matrix.size());
        vector<double> b(matrix.size());
        vector<double> x(matrix.size());
        y[0] = matrix[0][0];
        a[0] = -matrix[0][1] / y[0];
        b[0] = results[0][0] / y[0];
        // for the next rows of matrix
        for (size_t i = 1; i < matrix.size(); i++) {
            y[i] = matrix[i][i] + matrix[i][i - 1] * a[i - 1];
            a[i] = -matrix[i][i + 1] / y[i];
            b[i] = (results[i][0] - matrix[i][i - 1] * b[i - 1]) / y[i];
        }
        // for the x[n] element
        x[matrix.size() - 1] = b[matrix.size() - 1];
        // for the other x elements
        for (int i = (int) matrix.size() - 2; i >= 0; i--)
            x[i] = a[i] * x[i + 1] + b[i];
        return x;
    }
};

int main() {
    RunMethod rm;
    auto a = rm.run_through();
    for (auto &x: a) {
        cout << x << ' ';
    }
    return 0;
}