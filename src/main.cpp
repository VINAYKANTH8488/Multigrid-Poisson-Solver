#include <iostream>
#include <vector>
#include <cmath>

#include "poisson.h"
#include "smoothers.h"

using namespace std;

//=====================================================
// Matrix Vector Product
//=====================================================

vector<double> matvec(
    const vector<vector<double>>& A,
    const vector<double>& x)
{
    int N = A.size();

    vector<double> y(N, 0.0);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            y[i] += A[i][j] * x[j];
        }
    }

    return y;
}

//=====================================================
// Residual
// r = b - Au
//=====================================================

vector<double> residual(
    const vector<vector<double>>& A,
    const vector<double>& u,
    const vector<double>& b)
{
    vector<double> Au = matvec(A, u);

    int N = b.size();

    vector<double> r(N);

    for (int i = 0; i < N; i++)
    {
        r[i] = b[i] - Au[i];
    }

    return r;
}

//=====================================================
// L2 Norm
//=====================================================

double norm2(const vector<double>& v)
{
    double sum = 0.0;

    for (double x : v)
    {
        sum += x * x;
    }

    return sqrt(sum);
}

//=====================================================
// Main
//=====================================================

int main()
{
    int nx = 4;
    int ny = 4;

    int N = nx * ny;

    double hx = 1.0 / (nx + 1);
    double hy = 1.0 / (ny + 1);

    double ax = 1.0 / (hx * hx);
    double ay = 1.0 / (hy * hy);

    //-------------------------------------------------
    // Build Poisson Matrix
    //-------------------------------------------------

    vector<vector<double>> A =
    buildPoissonMatrix(nx, ny);
    
    //-------------------------------------------------
    // Exact Solution and RHS
    //-------------------------------------------------

    vector<double> uExact(N);
    vector<double> b(N);

    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            double x = (i + 1) * hx;
            double y = (j + 1) * hy;

            int p = idx(i, j, nx);

            uExact[p] = exact(x, y);

            b[p] = rhs(x, y);
        }
    }
//-------------------------------------------------
// Solve using Solvers
//-------------------------------------------------

vector<double> uJ(N,0.0);
vector<double> uGS(N,0.0);
vector<double> uSOR(N,0.0);

Jacobi(A,uJ,b,100);

GaussSeidel(A,uGS,b,100);

SOR(A,uSOR,b,100,1.7);

//-------------------------------------------------
// Residual after solving
//-------------------------------------------------

vector<double> rJ =
    residual(A,uJ,b);

vector<double> rGS =
    residual(A,uGS,b);

vector<double> rSOR =
    residual(A,uSOR,b);

//-------------------------------------------------
// Compute L2 Error
//-------------------------------------------------

double errorJ = 0.0;

for(int i=0;i<N;i++)
{
    errorJ +=
        pow(uJ[i]-uExact[i],2);
}

errorJ = sqrt(errorJ);

//-------------------------------------------------
// Output
//-------------------------------------------------
cout << "\nSolver Comparison\n";
cout << "---------------------\n";

cout << "Jacobi Residual      = "
     << norm2(rJ) << endl;

cout << "Gauss-Seidel Residual= "
     << norm2(rGS) << endl;

cout << "SOR Residual         = "
     << norm2(rSOR) << endl;
return 0;
}