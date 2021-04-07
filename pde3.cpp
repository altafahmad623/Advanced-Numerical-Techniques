#include <bits/stdc++.h>
using namespace std;
#define lli long long
#define PI 3.14159265358
// u_t = u_{xx} , u(x,0) = sin (pi x), 0 < x < 1; u(0,t) , u(1,t) =0
// using Crank-Nicolson scheme
double fin(double x, double t)
{
    return 0.0;
}
int main()
{
    double h = 0.05; // delta x
    double r = 0.5, mu = 1;
    double dt = r * (h * h) / (mu);
    double x_0 = 0.0;
    double x_n = 1.0;
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
    col.push_back(sin(PI * x_0));
    for (int i = 1; i <= N; i++)
    {
        xx[i] = xx[i - 1] + h;
        col.push_back(sin(PI * xx[i]));
    }
    for (int i = 1; i <= 2 * N; i++)
    {
        tt[i] = tt[i - 1] + dt;
    }

    u.push_back(col);
    for (int n = 1; n < N * 2; n++)
    {
        col[0] = 0;
        col[N] = 0;
        for (int i = 1; i < N; i++)
        {
            aa[i] = r / 2.0;
            bb[i] = -1.0 - r;
            cc[i] = r / 2.0;
            dd[i] = -fin(xx[i], tt[n - 1]) - ( r* (u[n - 1][i - 1] + u[n - 1][i + 1]) / 2.0) + ((r-1)*u[n-1][i] );
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
        printf("\t %d", j);
    }
    printf("\n");
    cout << "t =  |";
    for (int j = 0; j < N ; j++)
    {
        printf("\t %0.4lf", tt[j]);
    }
    printf("\n");
    for (int j = 0; j < 4 * N + 10; j++)
    {
        printf("--");
    }
    printf("\n");
    for (int j = 0; j <= N; j++)
    {
        cout << "x_" << j << "  |";
        for (int i = 0; i < N ; i++)
        {
            printf("\t %0.4lf", u[i][j]);
        }
        cout << "\n";
    }

    return 0;
}