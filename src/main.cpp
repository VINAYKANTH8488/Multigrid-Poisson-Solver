#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

int idx(int i,int j,int nx)
{
    return j*nx+i;
}

int main()
{
    int nx=4;
    int ny=4;

    int N=nx*ny;

    vector<vector<double>> A(
        N,
        vector<double>(N,0.0)
    );

    double hx=1.0/(nx+1);
    double hy=1.0/(ny+1);

    double ax=1.0/(hx*hx);
    double ay=1.0/(hy*hy);

    for(int j=0;j<ny;j++)
    {
        for(int i=0;i<nx;i++)
        {
            int p=idx(i,j,nx);

            A[p][p]=2.0*(ax+ay);

            if(i>0)
                A[p][idx(i-1,j,nx)] = -ax;

            if(i<nx-1)
                A[p][idx(i+1,j,nx)] = -ax;

            if(j>0)
                A[p][idx(i,j-1,nx)] = -ay;

            if(j<ny-1)
                A[p][idx(i,j+1,nx)] = -ay;
        }
    }

    cout<<"Poisson Matrix (first 10 rows)\n";

    for(int i=0;i<min(N,10);i++)
    {
        for(int j=0;j<min(N,10);j++)
        {
            cout<<A[i][j]<<" ";
        }
        cout<<"\n";
    }

    return 0;
}