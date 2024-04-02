#include <iostream>
#include <vector>
#include "../Common/LU_decomposition.hpp"

using namespace std;

class MNS {
private:
    vector<double> x = {-0.7, -0.4, -0.1, 0.2, 0.5, 0.8};
    vector<double> f = {-0.7754, -0.41152, -0.10017, 0.20136, 0.5236, 0.9273};
public:
    double error_sum_sq_1() {
        double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
        for (size_t i = 0; i < x.size(); i++) {
            sum_x += x[i];
            sum_y += f[i];
            sum_x2 += x[i] * x[i];
            sum_xy += x[i] * f[i];
        }
        vector<vector<double>> matrix = {{double(x.size()), sum_x},
                                         {sum_x,            sum_x2}};
        vector<vector<double>> matrix_result = {{sum_y},
                                                {sum_xy}};
        LU_Decomposition lu(matrix, matrix_result);
        auto x_polynome = lu.solve();
        double phi = 0.0;
        for (size_t i = 0; i < x.size(); i++)
            phi += (x_polynome[0] + x_polynome[1] * x[i] - f[i]) * (x_polynome[0] + x_polynome[1] * x[i] - f[i]);
        return phi;
    }

    double error_sum_sq_2() {
        double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0, sum_x3 = 0, sum_x4 = 0, sum_y_x2 = 0;
        for (size_t i = 0; i < x.size(); i++) {
            sum_x += x[i];
            sum_y += f[i];
            sum_x2 += x[i] * x[i];
            sum_xy += x[i] * f[i];
            sum_x3 += x[i] * x[i] * x[i];
            sum_x4 += x[i] * x[i] * x[i] * x[i];
            sum_y_x2 += x[i] * x[i] * f[i];
        }
        vector<vector<double>> matrix = {{double(x.size()), sum_x,  sum_x2},
                                         {sum_x,            sum_x2, sum_x3},
                                         {sum_x2,           sum_x3, sum_x4}};
        vector<vector<double>> matrix_result = {{sum_y},
                                                {sum_xy},
                                                {sum_y_x2}};
        LU_Decomposition lu(matrix, matrix_result);
        auto x_polynome = lu.solve();
        double phi = 0.0;
        for (size_t i = 0; i < x.size(); i++)
            phi += (x_polynome[0] + x_polynome[1] * x[i] + x_polynome[2] * x[i] * x[i] - f[i]) *
                   (x_polynome[0] + x_polynome[1] * x[i] + x_polynome[2] * x[i] * x[i] - f[i]);
        return phi;
    }

};

int main() {
    MNS ms;
    cout << ms.error_sum_sq_1();
    cout << ms.error_sum_sq_2();
}