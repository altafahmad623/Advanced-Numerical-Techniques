#include <bits/stdc++.h>
using namespace std;
#define lli long long
#define PI 3.14159265358
// \frac{ \partial u}{\partial t} +  \frac{ \partial u}{\partial x} = \alpha \frac{ \partial^2 T}{\partial x^2} \; , 0 < x < 1;
// u = 0.1 , \alpha = 1 , u(-5,t) = 1 , T(5,t) = 0, T(x, 0) = if(x < 0) 1 ; else 0 
// using Crank-Nicolson scheme
double fin(double x, double t)
{
    return 0.0;
}
int main()
{
    double h = 1; // delta x
    double dt = 0.1 , dx = h;
    double r = 0.5, mu = 1;
    double u2 = 1, alpha = 1;
    //double dt = r * (h * h) / (mu);
    double x_0 = -5.0;
    double x_n = 5.0;
    int N = (x_n - x_0) / h + 0.5;
    cout << "N = " << N << endl;
    vector<vector<double>> u;
    vector<double> col;
    double *xx = new double[N + 1];
    double *tt = new double[2 * N + 1];
    double *aa = new double[N + 1];
    double *bb = new double[N + 1];
    double *cc = new double[N + 1];
    double *dd = new double[N + 1];
    double *c_dash = new double[N + 1];
    double *d_dash = new double[N + 1];
    tt[0] = 0.0;
    xx[0] = x_0;
    col.push_back(1);
    for (int i = 1; i <= N; i++)
    {
        xx[i] = xx[i - 1] + h;
        if(xx[i] < 0)
            col.push_back(1);
        else
        {
            col.push_back(0);
        }
        
    }
    for (int i = 1; i <= 2 * N; i++)
    {
        tt[i] = tt[i - 1] + dt;
    }

    u.push_back(col);
    for (int n = 1; n < N * 2; n++)
    {
        col[0] = 1;
        col[N] = 0; // these are the initial conditions that nee
        for (int i = 1; i < N; i++)
        {
            aa[i] = (- u2 /(4*dx)) - (alpha)/(2 * dx * dx);
            bb[i] = (1/dt) + alpha/(dx * dx);
            cc[i] = (u2 /(4*dx)) - (alpha)/(2 * dx * dx);
            dd[i] = (alpha * (u[n-1][i-1] - 2*u[n-1][i] +u[n-1][i+1])/(2 * dx * dx) ) + (u[n-1][i]/dt ) - (u2*(u[n-1][i+1] - u[n-1][i-1])/(4 * dx) ) ;
        }
        c_dash[1] = cc[1] / bb[1];
        d_dash[1] = (dd[1] - aa[1] * col[0]) / bb[1];
        for (int i = 2; i < (N - 1); i++)
        {
            c_dash[i] = cc[i] / (bb[i] - aa[i] * c_dash[i - 1]);
            d_dash[i] = (dd[i] - aa[i] * d_dash[i - 1]) / (bb[i] - aa[i] * c_dash[i - 1]);
        }
        col[N - 1] = (dd[N - 1] - cc[N - 1] * col[N] - aa[N - 1] * d_dash[N - 2]) / (bb[N - 1] - aa[N - 1] * c_dash[N - 2]);

        for (int i = N - 2; i >= 1; i--)
        {
            //cout << "c' " << i << " =" << c_dash[i] << " ;d' " << i << " =" << d_dash[i] << endl;
            col[i] = d_dash[i] - c_dash[i] * col[i + 1];
        }
        u.push_back(col);
    }
    cout << "n -> |";
    for (int j = 0; j < N ; j++)
    {
        //cout<<"\t"<<j;
        printf("\t\t %d", j);
    }
    printf("\n");
    cout << "t =  |";
    for (int j = 0; j < N  ; j++)
    {
        //cout<<"\t"<<tt[j];
        printf("\t\t %0.4lf", tt[j]);
    }
    printf("\n");
    for (int j = 0; j < 8 * N + 10; j++)
    {
        printf("--");
    }
    printf("\n");
    for (int j = 0; j <= N; j++)
    {
        cout << "x_" << j << "  |";
        for (int i = 0; i < N ; i++)
        {
            //cout<<"\t"<<u[i][j];
            printf("\t\t%0.4lf", u[i][j]);
        }
        cout << "\n";
    }

    return 0;
}