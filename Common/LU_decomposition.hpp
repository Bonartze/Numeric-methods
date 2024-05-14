#include <iostream>
#include <vector>

using namespace std;

class LU_Decomposition {
private:
    vector<vector<double>> A;
    vector<vector<double>> B;

public:
    LU_Decomposition(const vector<vector<double>> &matrix, const vector<vector<double>> &result);

    void LU_decompose(vector<vector<double>> &L, vector<vector<double>> &U);

    vector<vector<double>> forward_substitution(const vector<vector<double>> &L);

    vector<vector<double>> backward_substitution(const vector<vector<double>> &U, const vector<vector<double>> &y);

    vector<double> solve();
};

