#include <bits/stdc++.h>
using namespace std;
#define lli long long
//change the functions A,B,C, and D which are in the equation: 
//A(x)y'' + B(x)y' + C(x)y = D(x), they can be functions of x or constants.
double A(double x)
{
    return 1.0;
}
double B(double x)
{
    return 0;
}
double C(double x)
{
    return -2.0;
}
double D(double x)
{
    return 0.0;
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
    double h = 0.2;
    double x_0 = 0.0;
    double x_n = 1.0;
    int n = (x_n - x_0) / h + 0.5;
    double *yy = new double[n + 1];
    yy[0] = 1.0;
    //yy[n] = 0.0566;
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
    
    
    // now we should apply the y_n-1 = y_n+1 in the last equation. so a_n(new ) = a_n + c_n
    aa[n] = aa[n] + cc[n];
    cc[n] = 0.0;
    for (int i = 1; i <= n; i++)
    {
        cout<<aa[i] <<" y_"<<i-1<<" + "<<bb[i] <<" y_"<<i<<" + "<<cc[i] <<" y_"<<i+1 <<" = "<<dd[i]<<endl; 
    }
    double *c_dash = new double[n + 1];
    double *d_dash = new double[n + 1];
    c_dash[1] = cc[1] / bb[1];
    d_dash[1] = (dd[1] - aa[1] * yy[0]) / bb[1];

    for (int i = 2; i <= (n ); i++)
    {
        c_dash[i] = cc[i] / (bb[i] - aa[i] * c_dash[i - 1]);
        d_dash[i] = (dd[i] - aa[i] * d_dash[i - 1]) / (bb[i] - aa[i] * c_dash[i - 1]);
    }
    //cout << dd[n] << " " << cc[n - 1] << " " << yy[n] << " " <<aa[n]<< " " <<d_dash[n-1]<< " " 
    //<< bb[n - 1] << " " << aa[n - 1] << " " << c_dash[n - 2] << endl;
    //yy[n ] = (dd[n] - cc[n] * yy[n] - aa[n] * d_dash[n - 1]) / (bb[n ] - aa[n ] * c_dash[n - 1]);
    yy[n] = d_dash[n];
    for (int i = n - 1; i >= 1; i--)
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