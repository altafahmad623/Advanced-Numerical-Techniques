#include <bits/stdc++.h>
using namespace std;
#define lli long long
#define PI 3.14159265358
//change the functions A,B,C, and D which are in the equation: A(x)y'' + B(x)y' + C(x)y = D(x),
//they can be functions of x or constants.
double f_i(double y_i_minus, double y_i, double y_i_plus, double h, double x_i)
{
    double f = (3 * y_i *((y_i_minus - 2 * y_i + y_i_plus) / (h * h)) ) +  (((y_i_plus - y_i_minus) / (2 * h))* ((y_i_plus - y_i_minus) / (2 * h))) ;
    return f;
}
double f_dash_y_minus(double y_i_minus, double y_i, double y_i_plus, double h, double x_i)
{
    double f = ((3 * y_i)/(h*h)) - ((y_i_plus - y_i_minus) / (2 * h * h));
    return f;
}
double f__dash_y(double y_i_minus, double y_i, double y_i_plus, double h, double x_i)
{
    double f = (3  *((y_i_minus - 4 * y_i + y_i_plus) / (h * h)) );
    return f;
}
double f_dash_y_plus(double y_i_minus, double y_i, double y_i_plus, double h, double x_i)
{
    double f = ((3 * y_i)/(h*h)) + ((y_i_plus - y_i_minus) / (2 * h * h));
    return f;
}
double absolute(double k)
{
    if(k < 0.0)
    {
        return -1* k;
    }
    else
    {
        return k;
    }
    
}
int main()
{
    double h = 0.3333333333333;
    double x_0 = 0.0;
    double x_n = 1;
    int n = (x_n - x_0) / h + 0.5;
    cout<<n<<endl;
    double *yy = new double[n + 1];
    yy[0] = 0.0;
    yy[n] = 1.0;
    double dely = (yy[n] - yy[0]) / n;
    double *xx = new double[n + 1];
    double *aa = new double[n + 1];
    double *bb = new double[n + 1];
    double *cc = new double[n + 1];
    double *dd = new double[n + 1];
    double *c_dash = new double[n + 1];
    double *d_dash = new double[n + 1];
    xx[0] = x_0;
    for (int i = 1; i <= n; i++)
    {
        xx[i] = xx[i - 1] + h;
        //cout<<xx[i]<<endl;
        yy[i] = yy[i - 1] + dely;
        //cout << aa[i] << " " << bb[i] << " " << cc[i] << endl;
    }
    cout << "Initial Guess is \n";
    for (int i = 0; i <= n; i++)
    {
        //cout << "y(" << xx[i] << ") = " << yy[i] << endl;
        cout << "y(" << xx[i] << ")^(0) = " << yy[i] <<" = y_"<<i<< endl;
    }
    double error = 0.00001;
    double diff = 1;
    double *Dell_y = new double[n + 1];
    Dell_y[0] = 0;
    Dell_y[n] = 0;
    int iter = 0;
    while (diff > error)
    {
        iter++;
        cout << "\n Iteration no: " << iter << " \n";
        for (int i = 1; i < n; i++)
        {
            aa[i] = f_dash_y_minus(yy[i - 1], yy[i], yy[i + 1], h, xx[i]);
            bb[i] = f__dash_y(yy[i - 1], yy[i], yy[i + 1], h, xx[i]);
            cc[i] = f_dash_y_plus(yy[i - 1], yy[i], yy[i + 1], h, xx[i]);
            dd[i] = -1 * f_i(yy[i - 1], yy[i], yy[i + 1], h, xx[i]);
        }
        for (int i = 1; i < n; i++)
        {
            cout << aa[i] << " del_y_" << i - 1 << " + " << bb[i] << "  del_y_" << i << " + " << cc[i] << " del_y_" << i + 1 << " = " << dd[i] << endl;
        }
        c_dash[1] = cc[1] / bb[1];
        d_dash[1] = (dd[1] - aa[1] * Dell_y[0]) / bb[1];
        for (int i = 2; i < (n - 1); i++)
        {
            c_dash[i] = cc[i] / (bb[i] - aa[i] * c_dash[i - 1]);
            d_dash[i] = (dd[i] - aa[i] * d_dash[i - 1]) / (bb[i] - aa[i] * c_dash[i - 1]);
        }
        Dell_y[n - 1] = (dd[n - 1] - cc[n - 1] * Dell_y[n] - aa[n - 1] * d_dash[n - 2]) / (bb[n - 1] - aa[n - 1] * c_dash[n - 2]);
        for (int i = n - 2; i >= 1; i--)
        {
            //cout << "c' " << i << " =" << c_dash[i] << " ;d' " << i << " =" << d_dash[i] << endl;
            Dell_y[i] = d_dash[i] - c_dash[i] * Dell_y[i + 1];
        }
        double ddc;
        for (int i = 1; i< n ; i ++ )
        {
            ddc = Dell_y[i-1] * aa[i] + Dell_y[i] * bb[i] + Dell_y[i+1] *cc[i] - dd[i];
            cout<<" The value of the Del equation " << i << " is " << ddc<<endl;
        }
        diff = 0.0000001;
        for (int i = 1; i < n; i++)
        {
            if ( absolute(Dell_y[i]) > diff)
            {
                diff = absolute(Dell_y[i]) ;
            }
            yy[i] = yy[i] + Dell_y[i];
        }
        cout << "\n Solution at iteration " << iter << " is :\n";
        for (int i = 0; i <= n; i++)
        {
            //cout << "y(" << xx[i] << ") = " << yy[i] << endl;
            cout << "y(" << xx[i] << ")^(0) = " << yy[i] <<" = y_"<<i<< endl;
        }
        cout << "\n And the values of the equations at " << iter << " is :\n";
        for (int i = 1; i< n ; i ++ )
        {
            ddc = f_i(yy[i-1],yy[i],yy[i+1],h,xx[i]);
            cout<<" The value of f _ " << i << " is " << ddc<<endl;
        }
        cout << "\n And maximum error is : " << diff << "\n";
        if (iter >= 20) // pseudo code for checking . remove this later
        {
            break;
        }
    }
    if (iter >= 20)
    {
        cout << "Did not converge :(\n";
    }
    else
    {
        cout << "The solution converged :)\n And the final answer is : \n";
        for (int i = 0; i <= n; i++)
        {
            //cout << "y(" << xx[i] << ") = " << yy[i] << endl;
            cout << "y(" << xx[i] << ")^(0) = " << yy[i] <<" = y_"<<i<< endl;
        }
        cout << "\n And the values of the equations at " << iter << " is :\n";
        double ddc;
        for (int i = 1; i< n ; i ++ )
        {
            cout<<" y_minus = "<<yy[i-1] <<" , y_i = " << yy[i] <<" y_plus = "<<yy[i+1]<<endl;
            /*ddc = (yy[i+1] - 2* yy[i] + yy[i-1] )/ (h*h);
            cout<<" a_1 = "<<ddc;
            ddc = ( yy[i] * (yy[i+1] - yy[i-1]) )/ (2*h);
            cout<<" a_2 = "<<ddc;
            ddc = - 4 - 4* xx[i] * xx[i];
            cout<<" a_1 = "<<ddc<<"\n";*/
            ddc = f_i(yy[i-1],yy[i],yy[i+1],h,xx[i]);
            cout<<" The value of f _ " << i << " is : " << ddc<<endl;
        }
    }

    return 0;
}