#include "LU_decomposition.hpp"

LU_Decomposition::LU_Decomposition(const vector<vector<double>> &matrix, const vector<vector<double>> &result) : A(
        matrix),
                                                                                                                 B(result) {}

void LU_Decomposition::LU_decompose(vector<vector<double>> &L, vector<vector<double>> &U) {
    int n = A.size();
    L.assign(n, vector<double>(n, 0));
    U.assign(n, vector<double>(n, 0));

    for (int i = 0; i < n; i++) {
        for (int k = i; k < n; k++) {
            double sum = 0.0;
            for (int j = 0; j < i; j++)
                sum += (L[i][j] * U[j][k]);
            U[i][k] = A[i][k] - sum;
        }

        for (int k = i; k < n; k++) {
            if (i == k)
                L[i][i] = 1;
            else {
                double sum = 0.0;
                for (int j = 0; j < i; j++)
                    sum += (L[k][j] * U[j][i]);
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }
}

vector<vector<double>> LU_Decomposition::forward_substitution(const vector<vector<double>> &L) {
    int n = B.size();
    vector<vector<double>> y(n, vector<double>(1, 0));
    for (int i = 0; i < n; i++) {
        y[i][0] = B[i][0];
        for (int j = 0; j < i; j++) {
            y[i][0] -= L[i][j] * y[j][0];
        }
        y[i][0] /= L[i][i];
    }
    return y;
}

vector<vector<double>>
LU_Decomposition::backward_substitution(const vector<vector<double>> &U, const vector<vector<double>> &y) {
    int n = U.size();
    vector<vector<double>> x(n, vector<double>(1, 0));
    for (int i = n - 1; i >= 0; i--) {
        x[i][0] = y[i][0];
        for (int j = i + 1; j < n; j++) {
            x[i][0] -= U[i][j] * x[j][0];
        }
        x[i][0] /= U[i][i];
    }
    return x;
}

vector<double> LU_Decomposition::solve() {
    vector<vector<double>> L, U;
    LU_decompose(L, U);

    auto y = forward_substitution(L);
    auto x = backward_substitution(U, y);

    vector<double> res;

    for (const auto &row: x) {
        for (const auto &elem: row)
            res.push_back(elem);
    }
    return res;
}