#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

//=====================================================
// Index Mapping
//=====================================================

int idx(int i, int j, int nx)
{
    return j * nx + i;
}

//=====================================================
// Exact Solution
// u = sin(pi*x)sin(pi*y)
//=====================================================

double exact(double x, double y)
{
    return sin(M_PI * x) * sin(M_PI * y);
}

//=====================================================
// RHS corresponding to
// -∇²u = f
//=====================================================

double rhs(double x, double y)
{
    return 2.0 * M_PI * M_PI *
           sin(M_PI * x) *
           sin(M_PI * y);
}

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

    vector<vector<double>> A(
        N,
        vector<double>(N, 0.0));

    for (int j = 0; j < ny; j++)
    {
        for (int i = 0; i < nx; i++)
        {
            int p = idx(i, j, nx);

            A[p][p] = 2.0 * (ax + ay);

            if (i > 0)
                A[p][idx(i - 1, j, nx)] = -ax;

            if (i < nx - 1)
                A[p][idx(i + 1, j, nx)] = -ax;

            if (j > 0)
                A[p][idx(i, j - 1, nx)] = -ay;

            if (j < ny - 1)
                A[p][idx(i, j + 1, nx)] = -ay;
        }
    }

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
    // Compute Residual
    //-------------------------------------------------

    vector<double> r =
        residual(A, uExact, b);

    //-------------------------------------------------
    // Output
    //-------------------------------------------------

    cout << "\n====================================\n";
    cout << " Manufactured Solution Test\n";
    cout << "====================================\n";

    cout << "Grid Size : "
         << nx << " x " << ny
         << endl;

    cout << "Residual Norm = "
         << norm2(r)
         << endl;

    cout << "====================================\n";

    return 0;
}