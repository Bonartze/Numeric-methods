#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double epsilon = 0.0000001;

class QR_decomposition {
private:
    vector<vector<double>> matrix = {{-9, 2, 2},
                                     {-2, 0, 7},
                                     {8,  2, 0}};
public:
    double scalar_mult(const vector<double> &v1, const vector<double> &v2) {
        double sum = 0.0;
        for (size_t i = 0; i < v1.size(); ++i) {
            sum += v1[i] * v2[i];
        }
        return sum;
    }

    vector<vector<double>> multiply_matrices(const vector<vector<double>> &A, const vector<vector<double>> &B) {
        size_t rowsA = A.size();
        size_t colsA = A[0].size();
        size_t colsB = B[0].size();

        vector<vector<double>> product(rowsA, vector<double>(colsB, 0));

        for (size_t i = 0; i < rowsA; ++i) {
            for (size_t j = 0; j < colsB; ++j) {
                for (size_t k = 0; k < colsA; ++k) {
                    product[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        return product;
    }

    vector<vector<double>> get_transposed_matrix(const vector<vector<double>> &m) {
        vector<vector<double>> trans(m[0].size(), vector<double>(m.size()));
        for (size_t i = 0; i < m.size(); ++i) {
            for (size_t j = 0; j < m[i].size(); ++j) {
                trans[j][i] = m[i][j];
            }
        }
        return trans;
    }

    vector<double> projection(const vector<double> &b, const vector<double> &a) {
        auto k = scalar_mult(a, b) / scalar_mult(b, b);
        vector<double> res(b.size());
        for (size_t i = 0; i < b.size(); i++)
            res[i] = b[i] * k;
        return res;
    }

    double get_norm(const vector<vector<double>> &m) {
        double norm = 0.0;
        for (size_t i = 0; i < m.size(); i++) {
            for (size_t j = 0; j < m[0].size(); j++) {
                if (i != j)
                    norm += m[i][j] * m[i][j];
            }
        }
        return sqrt(norm);
    }

    void get_normalized(vector<double> &norm) {
        double n = 0.0;
        for (auto &a: norm)
            n += a * a;
        n = sqrt(n);
        for (auto &a: norm)
            a /= n;
    }

    vector<vector<double>> get_Q(vector<vector<double>> m) {
        vector<vector<double>> result(m[0].size(), vector<double>(m.size()));
        auto t_m = get_transposed_matrix(m);
        for (size_t i = 0; i < m.size(); i++) {
            for (size_t j = 0; j < m.size(); j++) {
                result[i][j] += t_m[i][j];
                for (size_t k = 0; k < i; k++)
                    result[i][j] -= projection(result[k], t_m[i])[j];
            }
            get_normalized(result[i]);
        }
        return get_transposed_matrix(result);
    }

    vector<vector<double>> geg_QR() {

        bool flag = true;
        while (flag) {
            auto Q = get_Q(matrix);
            auto R = multiply_matrices(get_transposed_matrix(Q), matrix);
            matrix = multiply_matrices(R, Q);
            flag = false;
            for (size_t i = 0; i < matrix.size(); i++)
                for (size_t j = 0; j < matrix.size(); j++) {
                    if (i != j && matrix[i][j] > fabs(epsilon))
                        flag = true;
                }
        }
        return matrix;
    }
};

int main() {
    QR_decomposition q_r;
    auto a = q_r.geg_QR();
    for (size_t i = 0; i < a.size(); i++)
        cout << a[i][i] << ' ';
}