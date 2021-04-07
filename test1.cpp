#include <bits/stdc++.h>
using namespace std;
#define lli long long
//change the functions A,B,C, and D which are in the equation: A(x)y''
// + B(x)y' + C(x)y = D(x), they can be functions of x or constants.
double A(double x)
{
    return (1-x);
}
double B(double x)
{
    return x;
}
double C(double x)
{
    return -1.0;
}
double D(double x)
{
    return (1-x)*(1-x);
}
double a_i(double x, double h)
{
    double aii = (A(x) / (h * h)) - (B(x) / (2.0 * h));
    return aii;
}
double b_i(double x, double h)
{
    double aii = -(2.0 * A(x) / (h * h)) + C(x);
    return aii;
}
double c_i(double x, double h)
{
    double aii = (A(x) / (h * h)) + (B(x) / (2.0 * h));
    return aii;
}
int main()
{
    double h = 0.25; // provide the value of h
    double x_0 = 0.0; // the initial values of x 
    double x_n = 1.0;
    int n = (x_n - x_0) / h + 0.5;
    double *yy = new double[n + 1];
    yy[0] = 1.0; // and the boundary conditions
    yy[n] = 3.0;
    double *xx = new double[n + 1];
    double *aa = new double[n + 1];
    double *bb = new double[n + 1];
    double *cc = new double[n + 1];
    double *dd = new double[n + 1];
    xx[0] = x_0;
    for (int i = 1; i <= n; i++)
    {
        xx[i] = xx[i - 1] + h;
        //cout<<xx[i]<<endl;
        aa[i] = a_i(xx[i], h);
        bb[i] = b_i(xx[i], h);
        cc[i] = c_i(xx[i], h);
        dd[i] = D(xx[i]);
        //cout << aa[i] << " " << bb[i] << " " << cc[i] << endl;
    }
    for (int i = 1; i <= n; i++)
    {
        cout<<aa[i] <<" y_"<<i-1<<" + "<<bb[i] <<" y_"<<i<<" + "<<cc[i] <<" y_"<<i+1 <<" = "<<dd[i]<<endl; 
    }
    double *c_dash = new double[n + 1];
    double *d_dash = new double[n + 1];
    c_dash[1] = cc[1] / bb[1];
    d_dash[1] = (dd[1] - aa[1] * yy[0]) / bb[1];

    for (int i = 2; i < (n - 1); i++)
    {
        c_dash[i] = cc[i] / (bb[i] - aa[i] * c_dash[i - 1]);
        d_dash[i] = (dd[i] - aa[i] * d_dash[i - 1]) / (bb[i] - aa[i] * c_dash[i - 1]);
    }
    //cout << dd[n] << " " << cc[n - 1] << " " << yy[n] << " " <<aa[n]<< " " <<d_dash[n-1]<< " "
    // << bb[n - 1] << " " << aa[n - 1] << " " << c_dash[n - 2] << endl;
    yy[n - 1] = (dd[n-1] - cc[n - 1] * yy[n] - aa[n - 1] * d_dash[n - 2]) / (bb[n - 1] - aa[n - 1] * c_dash[n - 2]);

    for (int i = n - 2; i >= 1; i--)
    {
        //cout << "c' " << i << " =" << c_dash[i] << " ;d' " << i << " =" << d_dash[i] << endl;
        yy[i] = d_dash[i] - c_dash[i] * yy[i + 1];
    }
    for (int i = 0; i <= n; i++)
    {
        //cout << "y(" << xx[i] << ") = " << yy[i] << endl;
        cout<<xx[i] <<";"<<yy[i]<<endl;
    }

    return 0;
}