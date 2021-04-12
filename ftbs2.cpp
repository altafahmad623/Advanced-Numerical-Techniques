#include <bits/stdc++.h>
using namespace std;
#define lli long long
#define PI 3.14159265358
// \frac
// mu = 0.1
double fin(double x, double t)
{
    return 0.0;
}
int main()
{
    double h = 0.1; // delta x
    double u2 = 0.1, mu = 0.1;
    double dt = ((mu * h) / u2), dx = h;

    double x_0 = 0.0;
    double x_n = 1.0;
    int N = (x_n - x_0) / h + 0.5;
    cout << "N = " << N << endl;
    vector<vector<double>> u;
    vector<double> col;
    double *xx = new double[100];
    double *tt = new double[200];
    tt[0] = 0.0;
    xx[0] = x_0;
    col.push_back((200 * x_0));
    for (int i = 1; i <= N; i++)
    {
        xx[i] = xx[i - 1] + h;
        if (xx[i] <= 0.5)
            col.push_back((xx[i] * 200));
        else
            col.push_back(200 * (1 - xx[i]));
    }
    for (int i = 1; i <= 2 * N; i++)
    {
        tt[i] = tt[i - 1] + dt;
    }
    //cout<<"dfsfsf";

    u.push_back(col);
    for (int n = 1; n <= 2 * N; n++)
    {
        col[0] = 0;
        for (int i = 1; i <= N; i++)
        {
            col[i] = ((1 - mu) * u[n - 1][i]) + (mu * u[n - 1][i - 1]);
        }
        u.push_back(col);
    }

    cout << "n -> |";
    for (int j = 0; j < N; j++)
    {
        printf("\t %d", j);
    }
    printf("\n");
    cout << "t =  |";
    for (int j = 0; j < N; j++)
    {
        printf("\t%0.5lf", tt[j]);
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
        for (int i = 0; i < N; i++)
        {
            printf("\t%0.3lf", u[i][j]);
        }
        cout << "\n";
    }
    return 0;
}