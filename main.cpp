#include <iostream>
#include <cmath>
#include <vector>
using namespace std;
double my_function(double x){
    return x;
}
double kern(double x, double s){
    return sin(x + s);
}
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
void print_solution(vector<double> grid, vector<double> solution, const int n){
    cout << "Node:" << "        "<< "Value:" << endl;
    for (int i = 0; i < n; i++){
        cout << grid.at(i) << "      " <<solution.at(i) << endl;
    }
    cout << endl;
}
vector<double> linspace(double a, double b, const int n){
    vector<double> grid(n);
    for (int i = 0; i < n; i++) {
        grid.at(i) = (b - a) / (n - 1) * i;
    }
    return grid;
}
vector<double> fill_f(vector<double> grid, const int n){
    vector<double> f(n);
    for(int i = 0; i < n; i++){
        f.at(i) = my_function(grid.at(i));
    }
    return f;
}

vector<vector<double>> fill_matrix(double a, double b, vector<double> grid, double lambda, const int n){
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
            A[i][j] = -h/3 * lambda * w * kern(grid.at(i), grid.at(j));
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
    double a = 0;
    double b = M_PI_2;
    double lambda = 1;
    const int n = 11; //кол-во точек
    auto grid = linspace(a, b, n);
    auto f = fill_f(grid, n);
    auto A = fill_matrix(a, b, grid, lambda, n);
    auto p = LU(A, n);
    auto L = p.first;
    auto U = p.second;
    auto y = solve_ltm(L, f, n);
    auto x = solve_utm(U, y, n);
    auto s = toch_solution_grid(grid, n);
    auto error = check_error(x, s, n);
    print_solution(grid, s, n);
    cout << "Max error = "<< max_error(error, n) << endl;
    return 0;
}