#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

class Diff {
private:
    vector<double> x = {-1.0, 0.0, 1.0, 2.0, 3.0};
    vector<double> f = {-0.7854, 0.0, 0.78540, 1.1071, 1.249};
public:
    double get_first_diff_second_accuracy(double val) {
        size_t i;
        for (i = 0; i < x.size() - 1; i++) {
            if (x[i] < val && x[i + 1] >= val)
                break;
        }
        return ((f[i + 1] - f[i]) / (x[i + 1] - x[i])) +
               ((((f[i + 2] - f[i + 1]) / (x[i + 2] - x[i + 1])) - (f[i + 1] - f[i]) / (x[i + 1] - x[i])) /
                (x[i + 2] - x[i])) * (2 * val - x[i] - x[i + 1]);
    }

    double get_second_diff(double val) {
        size_t i;
        for (i = 0; i < x.size() - 1; i++) {
            if (x[i] < val && x[i + 1] >= val)
                break;
        }
        return 2 * (((f[i + 2] - f[i + 1]) / (x[i + 2] - x[i + 1])) - ((f[i + 1] - f[i]) / (x[i + 1] - x[i]))) /
               (x[i + 2] - x[i]);
    }
};

int main() {
    Diff df;
    cout << df.get_first_diff_second_accuracy(1.0);
}