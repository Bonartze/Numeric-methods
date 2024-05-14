#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <initializer_list>
#include <iterator>

using namespace std;

class InterpolationPolinomsLagrange {
private:
    const vector<double> x_ = {-0.4, -0.1, 0.2, 0.5};

    double get_omega(double x) {
        double mult = 1;
        for (auto &x_el: x_)
            mult *= (x - x_el);
        return mult;
    }

    double get_omega_st(double x) {
        double mult = 1;
        for (auto &x_el: x_)
            if (x != x_el)
                mult *= (x - x_el);
        return mult;
    }

public:
    double f(double x) {
        return asin(x);
    }

    double get_L(double x) {
        double sum = 0;
        for (size_t i = 0; i < x_.size(); i++) {
            sum += f(x_[i]) * get_omega(x) / ((x - x_[i]) * get_omega_st(x_[i]));
        }
        return sum;
    }

    void print_L() {
        cout << "Lagrange Interpolation Polynomial:" << endl;
        for (size_t i = 0; i < x_.size(); i++) {
            cout << f(x_[i]) << " * ";
            for (size_t j = 0; j < x_.size(); j++) {
                if (i != j) {
                    cout << "(x - " << x_[j] << ") ";
                }
            }
            cout << "/ ";
            for (size_t j = 0; j < x_.size(); j++) {
                if (i != j) {
                    cout << "(" << x_[i] << " - " << x_[j] << ") ";
                }
            }
            if (i != x_.size() - 1) cout << "+ ";
        }
        cout << endl;
    }
};

class list_index : public list<double> {
public:
    list_index(std::initializer_list<double> init) : std::list<double>(init) {}

    template<typename Iter>
    list_index(Iter begin, Iter end) : std::list<double>(begin, end) {}

    double &operator[](size_t i) {
        if (i >= this->size()) {
            throw std::out_of_range("Index out of range");
        }

        auto it = this->begin();
        std::advance(it, i);

        return *it;
    }
};

class InterpolationPolinomsNewton {
private:

public:
    double f(list_index x) {
        if (x.size() == 1)
            return asin(x[0]);
        auto xl = x;
        auto xr = x;
        auto xi = xl.back();
        xl.pop_back();
        auto xk = xr.front();
        xr.pop_front();
        return (f(xl) - f(xr)) / (xk - xi);
    }

    double get_P(double x) {
        list_index x_ = {-0.4, -0.1, 0.2, 0.5};
        double p = f({x_[0]});
        for (size_t i = 1; i < x_.size(); i++) {
            double temp = 1;
            list_index push_x{};
            for (size_t j = 0; j < i; j++) {
                temp *= (x - x_[j]);
                push_x.push_back(x_[j]);
            }
            push_x.push_back(x_[i]);
            p += f(push_x) * temp;
        }
        return p;
    }

    void print_P() {
        list_index x_ = {-0.4, -0.1, 0.2, 0.5};
        cout << "Newton Interpolation Polynomial:" << endl;
        cout << f({x_[0]});
        for (size_t i = 1; i < x_.size(); i++) {
            cout << " + ";
            cout << f(list_index(x_.begin(), std::next(x_.begin(), i + 1)));
            for (size_t j = 0; j < i; j++) {
                cout << " * (x - " << x_[j] << ")";
            }
        }
        cout << endl;
    }
};

int main() {
    InterpolationPolinomsNewton ip;
    InterpolationPolinomsLagrange lg;
    cout << "P(0.1) using Newton: " << ip.get_P(0.1) << endl;
    cout << "L(0.1) using Lagrange: " << lg.get_L(0.1) << endl;

    ip.print_P();
    lg.print_L();

    return 0;
}
