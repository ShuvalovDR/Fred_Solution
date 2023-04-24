#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <functional>
#include <fstream>
using namespace std;


double toch_solution(double x){
    return x - (8 * (M_PI - 1))/(pow(M_PI, 2) - 4) * sin(x) - 2 * (4 - 2 * M_PI + pow(M_PI, 2))/(pow(M_PI, 2) - 4) * cos(x);
}
vector<double> toch_solution_grid(vector<double> grid, const int n){
    vector<double> sol_grid(n);
    for (int i = 0; i < n; i++) {
        sol_grid.at(i) = toch_solution(grid.at(i));
    }
    return sol_grid;
}
vector<double> check_error(vector<double> a, vector<double> b, const int n){
    vector<double> error(n);
    for (int i = 0; i < n; i++){
        error.at(i) = abs(a.at(i) - b.at(i));
    }
    return error;
}
double max_error(vector<double> a, const int n){
    double max_value = a.at(0);
    for (int i = 1; i < n; i ++){
        if (a.at(i) > max_value){
            max_value = a.at(i);
        }
    }
    return max_value;
}

vector<double> linspace(double a, double b, const int n){
    vector<double> grid(n);
    for (int i = 0; i < n; i++) {
        grid.at(i) = (b - a) / (n - 1) * i;
    }
    return grid;
}
vector<double> fill_f(vector<double> grid, const int n, function<double(double, double)> my_f){
    vector<double> f(n);
    for(int i = 0; i < n; i++){
        f.at(i) = my_f(grid.at(i), 0);
    }
    return f;
}

vector<vector<double>> fill_matrix(double a, double b, vector<double> grid, double lambda, const int n, function<double(double, double)> my_kern){
    vector<vector<double>> A(n, vector<double>(n));
    double h = (b - a)/(n - 1);
    double w;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if ((j == 0) or (j == n - 1)){
                w = 1;
            }
            else if(j % 2 == 0){
                w = 2;
            }
            else{
                w = 4;
            }
            A[i][j] = -h/3 * lambda * w * my_kern(grid.at(i), grid.at(j));
        }
        A[i][i] += 1;
    }
    return A;
}
pair <vector<vector<double>>, vector<vector<double>>> LU(vector<vector<double>> A, const int n){
    vector<vector<double>> L(n, vector<double>(n));
    vector<vector<double>> U(n, vector<double>(n));
    pair <vector<vector<double>>, vector<vector<double>>> p;
    double total;
    for (int m = 0; m < n; m++){
        for (int j = m; j < n; j++){
            total = 0;
            for (int k = 0; k < m; k++){
                total += L[m][k] * U[k][j];
            }
            U[m][j] = A[m][j] - total;
        }
        for (int i = m; i < n; i++){
            total = 0;
            for(int k = 0; k < m; k++){
                total += L[i][k] * U[k][m];
            }
            L[i][m] = (A[i][m] - total)/U[m][m];
        }
    }
    p.first = L;
    p.second = U;
    return p;
}
vector<double> solve_utm(vector<vector<double>> A, vector<double> b, const int n){
    vector<double> x(n);
    x.at(n - 1) = b.at(n - 1)/A[n - 1][n - 1];
    double total;
    for (int i = n - 2; i >= 0; i--){
        total = 0;
        for (int k = i + 1; k < n; k++){
            total += A[i][k] * x.at(k);
        }
        x.at(i) = (b.at(i) - total)/A[i][i];
    }
    return x;
}
vector<double> solve_ltm(vector<vector<double>> A, vector<double> b, const int n){
    vector<double> x(n);
    x.at(0) = b.at(0)/A[0][0];
    double total;
    for (int i = 1; i < n; i++){
        total = 0;
        for (int k = 0; k < i; k++){
            total += A[i][k] * x.at(k);
        }
        x.at(i) = (b.at(i) - total)/A[i][i];
    }
    return x;
}

int main(){
    map <string, function<double(double, double)>> parameters{
            {"a", [](double x, double s) {return double(0);}},
            {"b", [](double x, double s) {return double(M_PI_2);}},
            {"lambda", [](double x, double s) {return double(1);}},
            {"n", [](double x, double s) {return double(50);}},
            {"f", [] (double x, double s) {return x;}},
            {"kern", [] (double x, double s) {return sin(x + s);}}
    };
    auto a = parameters["a"](0, 0);
    double b = parameters["b"](0, 0);
    double lambda = parameters["lambda"](0, 0);
    double n = parameters["n"](0, 0);; //кол-во точек (будет использована квадратурная формула Симпсона, поэтому узлов будет 2n - 1)
    n = round(n);
    n = 2 * n - 1;
    auto grid = linspace(a, b, n);
    auto my_f = parameters["f"];
    auto f = fill_f(grid, n, my_f);
    auto my_kern = parameters["kern"];
    auto A = fill_matrix(a, b, grid, lambda, n, my_kern);
    auto p = LU(A, n);
    auto L = p.first;
    auto U = p.second;
    auto y = solve_ltm(L, f, n);
    auto x = solve_utm(U, y, n);
    auto s = toch_solution_grid(grid, n);
    auto error = check_error(x, s, n);
    double maximum_error = max_error(error, n);
    ofstream fout;
    fout.open("Solution.txt");
    fout.clear();
    fout << "First column - grid, second column - solution:" << endl;
    for (int i = 0; i < n; i++){
        fout << grid.at(i) << "          " <<x.at(i) << endl;
    }
    fout << "Maximum error is : " << maximum_error << endl;
    fout.close();
    return 0;
}