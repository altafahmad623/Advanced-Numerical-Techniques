#include <bits/stdc++.h>
using namespace std;
#define lli long long
#define PI 3.14159265358
// u_t + u * u_x = 0 ; u(x,0) = {sqrt(x) \; 0 <= x <= 1} else u(x,0) = 0
// mu = 1/2
double F(double x)
{
    return ((x*x)/2.0);
}
int main()
{
    double h = 0.2; // delta x
    double dt = 0.5*0.2, dx = h;
    double r = 0.5, mu = 0.5;
    //double dt = r * (h * h) / (mu);
    double x_0 = 0.0;
    double x_n = 1.0;
    int N = (x_n - x_0) / h + 0.5;
    cout << "N = " << N << endl;
    double** u =new double*[100];
    for (int i = 0; i < 99; i++)
    {
        u[i] = new double[100];
    }
    
    double *col = new double[100];
    double *xx = new double[100];
    double *tt = new double[200];
    tt[0] = 0.0;
    xx[0] = x_0;
    u[0][0] = sqrt(x_0);
    for (int i = 1; i <= N; i++)
    {
        xx[i] = xx[i - 1] + h;
        u[0][i] = sqrt(xx[i]);
    }
    for (int i = 1; i <= 2 * N; i++)
    {
        tt[i] = tt[i - 1] + dt;
    }
    for (int n = 1; n <= 2* N; n++)
    {
        for (int i = 0; i < N; i++)
        {
            col[i] = ((u[n-1][i] + u[n-1][i+1])/2.0 ) - 0.5 *0.5 * (F(u[n-1][i+1]) -F(u[n-1][i]) );
            cout<<"u_bar["<<n<<"] ["<<i<<"] = "<<col[i]<<" ";
        }
        cout<<"\n";
        u[n][0] = 0;
        for (int i = 1; i < N; i++)
        {
            u[n][i] = u[n-1][i] - mu*(F(col[i]) - F(col[i-1]));
        }
    }
    
    cout << "n -> |";
    for (int j = 0; j < N * 2; j++)
    {
        printf("\t %d", j);
    }
    printf("\n");
    cout << "t =  |";
    for (int j = 0; j < N * 2 ; j++)
    {
        printf("\t%0.5lf", tt[j]);
    }
    printf("\n");
    for (int j = 0; j < 8 * N + 10; j++)
    {
        printf("--");
    }
    printf("\n");
    for (int j = 0; j < N; j++)
    {
        cout << "x_" << j << "  |";
        for (int i = 0; i < N*2 ; i++)
        {
            printf("\t%0.5lf", u[i][j]);
        }
        cout << "\n";
    }
    return 0;
}