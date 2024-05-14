#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

class Diff {
private:
    vector<double> x = {0.0, 0.1, 0.2, 0.3, 0.4};
    vector<double> f = {1.0, 1.1052, 1.2214, 1.3499, 1.4918};
public:
    double get_first_diff_second_accuracy_left(double val) {
        size_t i;
        for (i = 0; i < x.size() - 1; i++) {
            if (x[i] < val && x[i + 1] >= val)
                break;
        }
        return ((f[i + 1] - f[i]) / (x[i + 1] - x[i])) +
               ((((f[i + 2] - f[i + 1]) / (x[i + 2] - x[i + 1])) - (f[i + 1] - f[i]) / (x[i + 1] - x[i])) /
                (x[i + 2] - x[i])) * (2 * val - x[i] - x[i + 1]);
    }

    double get_first_diff_second_accuracy_right(double val) {
        size_t i;
        for (i = 0; i < x.size() - 2; i++) {
            if (x[i] <= val && x[i + 1] > val)
                break;
        }
        return ((f[i + 2] - f[i + 1]) / (x[i + 2] - x[i + 1])) +
               ((((f[i + 2] - f[i + 1]) / (x[i + 2] - x[i + 1])) - (f[i + 1] - f[i]) / (x[i + 1] - x[i])) /
                (x[i + 2] - x[i])) * (2 * val - x[i + 1] - x[i + 2]);
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
    double val = 0.2;

    cout << "Left-side first derivative with second-order accuracy at " << val << " is: "
         << df.get_first_diff_second_accuracy_left(val) << endl;

    cout << "Right-side first derivative with second-order accuracy at " << val << " is: "
         << df.get_first_diff_second_accuracy_right(val) << endl;

    cout << "Second derivative at " << val << " is: "
         << df.get_second_diff(val) << endl;

    return 0;
}
