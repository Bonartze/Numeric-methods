#include <iostream>
#include <vector>
#include <iomanip>
/*Построить кубический сплайн для функции, заданной в узлах
 * интерполяции, предполагая, что сплайн имеет нулевую кривизну
 * при  и вычислить значение функции в точке */
using namespace std;

/*Суть сплайн-интерполяции
заключается в определении интерполирующей функции по формулам одного типа для
различных непересекающихся промежутков и в стыковке значений функции и её
производных на их границах.*/
class RunMethod {
public:
    vector<vector<double>> matrix;
    vector<vector<double>> results;

    RunMethod(int size) : matrix(size, vector<double>(size, 0)), results(size, vector<double>(1, 0)) {}

    vector<double> run_through() {
        int n = matrix.size();
        vector<double> c(n);
        for (int i = 1; i < n; ++i) {
            double m = matrix[i][i - 1] / matrix[i - 1][i - 1];
            matrix[i][i] -= m * matrix[i - 1][i];
            results[i][0] -= m * results[i - 1][0];
        }
        c[n - 1] = results[n - 1][0] / matrix[n - 1][n - 1];
        for (int i = n - 2; i >= 0; --i) {
            c[i] = (results[i][0] - matrix[i][i + 1] * c[i + 1]) / matrix[i][i];
        }
        return c;
    }
};
/*между людбыми 2 точками разбиения строится многочлен n-ой степени  */


/*Его коэффициенты находятся из
условий равенства в узлах сетки значений сплайна и приближаемой функции, а также
равенства производных соответствующих многочленов.*/
class CubedSpline {
private:
    vector<double> f = {-0.41152, -0.10017, 0.20136, 0.52360, 0.92730};
    vector<double> x = {-0.4, -0.1, 0.2, 0.5, 0.8};
    vector<double> a, b, c, d;

    double h(int i) { return x[i] - x[i - 1]; }

public:
    CubedSpline() {
        int n = x.size();
        RunMethod mx(n);

        for (int i = 1; i < n - 1; ++i) {
            mx.matrix[i][i - 1] = h(i);
            mx.matrix[i][i] = 2 * (h(i) + h(i + 1));
            mx.matrix[i][i + 1] = h(i + 1);
            mx.results[i][0] = 3 * ((f[i + 1] - f[i]) / h(i + 1) - (f[i] - f[i - 1]) / h(i));
        }

        mx.matrix[0][0] = 1;
        mx.results[0][0] = 0;
        mx.matrix[n - 1][n - 1] = 1;
        mx.results[n - 1][0] = 0;

        c = mx.run_through();
        a = vector<double>(n);
        b = vector<double>(n);
        d = vector<double>(n);
        for (int i = 1; i < n; ++i) {
            a[i] = f[i - 1];
            d[i] = (c[i] - c[i - 1]) / (3 * h(i));
            b[i] = (f[i] - f[i - 1]) / h(i) - h(i) * (c[i] + 2 * c[i - 1]) / 3;
        }
    }

/* сплайн непрерывен вместе со своими n-1 производными*/
    double f_cubed(double val_x) {
        int n = x.size();
        if (val_x <= x[0]) return f[0];
        if (val_x >= x[n - 1]) return f[n - 1];

        int i = 1;
        for (; i < n; ++i) {
            if (x[i] > val_x) break;
        }
        double dx = val_x - x[i - 1];
        return a[i] + b[i] * dx + c[i - 1] * dx * dx + d[i] * dx * dx * dx;
    }

/*Для построения кубического сплайна необходимо построить n многочленов
третьей степени, т.е. определить 4n неизвестных . Эти коэффициенты ищутся из
условий в узлах сетки.*/
    void print_spline() {
        cout << "Spline coefficients:" << endl;
        cout << "i\t[xi-1, xi]\t\t\ta[i]\tb[i]\tc[i]\td[i]" << endl;
        for (size_t i = 1; i < x.size(); ++i) {
            cout << i << "\t" << "[" << x[i - 1] << "," << x[i] << "]\t" << a[i] << "\t" << b[i] << "\t" << c[i - 1]
                 << "\t" << d[i] << endl;
        }
    }
};
/*для построения кубического сплайна строится n многочленов 3 степени(4n неизвестных)*/
int main() {
    CubedSpline cb;
    cout.precision(5);
    cout << fixed << "f(0.1) = " << cb.f_cubed(0.1) << endl;
    cb.print_spline();
    return 0;
}
