#include <vector>
#include <iostream>
#include <cmath>

using namespace std;
double epsilon = 0.0000000001;

class VibrationMethod {
private:
    vector<vector<double>> matrix = {{9,  -2, 3},
                                     {-2, 6,  8},
                                     {3,  8,  -6}};
public:
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


    vector<vector<double>> transposed_matrix(const vector<vector<double>> &m) {
        vector<vector<double>> trans(m[0].size(), vector<double>(m.size()));
        for (size_t i = 0; i < m.size(); ++i) {
            for (size_t j = 0; j < m[i].size(); ++j) {
                trans[j][i] = m[i][j];
            }
        }
        return trans;
    }

    pair<vector<vector<double>>, vector<vector<double>>> vibrate_method() {
        vector<vector<double>> q(matrix.size(), vector<double>(matrix.size(), 0.0));
        for (size_t i = 0; i < matrix.size(); ++i)
            q[i][i] = 1.0;
        bool flag = true;
        while (flag) {
            auto a = matrix[0][1];
            auto positions = make_pair(0, 1);
            for (size_t i = 0; i < matrix.size(); i++) {
                for (size_t j = i + 1; j < matrix.size(); j++) {
                    if (fabs(a) < fabs(matrix[i][j])) {
                        a = fabs(matrix[i][j]);
                        positions = make_pair(i, j);
                    }
                }
            }
            double teta;
            if (matrix[positions.second][positions.second] != matrix[positions.first][positions.first])
                teta = 0.5 * atan2(2 * matrix[positions.first][positions.second],
                                   matrix[positions.first][positions.first] -
                                   matrix[positions.second][positions.second]);
            else
                teta = M_PI / 4;

            double c = cos(teta);
            double s = sin(teta);
            vector<vector<double>> J(matrix.size(), vector<double>(matrix.size(), 0.0));
            for (size_t i = 0; i < matrix.size(); ++i)
                J[i][i] = 1.0;
            J[positions.second][positions.second] = c;
            J[positions.first][positions.first] = c;
            J[positions.first][positions.second] = -s;
            J[positions.second][positions.first] = s;
            matrix = multiply_matrices(multiply_matrices(transposed_matrix(J), matrix), J);
            q = multiply_matrices(q, J);
            flag = false;
            for (size_t i = 0; i < matrix.size(); i++)
                for (size_t j = 0; j < matrix.size(); j++) {
                    if (i != j && matrix[i][j] > fabs(epsilon))
                        flag = true;
                }
        }
        return make_pair(matrix, q);
    }
};

int main() {
    VibrationMethod vb;
    auto a = vb.vibrate_method();
    for (size_t i = 0; i < a.first.size(); i++)
        cout << a.first[i][i] << ' ';
    cout << endl << endl;
    for (size_t i = 0; i < a.second.size(); i++) {
        for (size_t j = 0; j < a.second.size(); j++)
            cout << a.second[i][j] << ' ';
        cout << endl;
    }
}