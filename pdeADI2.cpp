#include <bits/stdc++.h>
using namespace std;
#define lli long long
#define PI 3.14159265358
// \frac{ \partial u}{\partial t}  = \nu [ \frac{ \partial^2 u}{\partial x^2} +  \frac{ \partial^2 u}{\partial y^2}] \; , 0 < x, y  < 1;
// \nu = 1,  u(x,y, 0 ) = sin (pi * x )* sin (pi * y ) \And u(t,  (1&0) , y) = 0 , u ( t, x,  (1&0) ) = 0
// using ALternating direction implicit scheme
double init(double x, double y)
{
    double k = (sin((x * PI))) * (sin((y * PI) ));
    return k;
}

int main()
{
    double h = 0.25; // delta x
    double dt = 0.1, dx = h, dy = 0.25;
    int timesteps = 4;
    double u2 = 0.1, alpha = 0.01;
    //double dt = r * (h * h) / (mu);
    double x_0 = 0.0, y_0 = 0;
    double x_n = 1.0, y_n = 1;
    int N = (x_n - x_0) / h + 0.5;
    int M = (y_n - y_0) / dy + 0.5;
    cout << "N = " << N << ", M = " << M << endl;
    vector<vector<vector<double>>> u;
    vector<vector<double>> instant;

    double *xx = new double[N + 1];
    double *yy = new double[M + 1];
    double *tt = new double[2 * N + 1];
    double *aa = new double[N + 1];
    double *bb = new double[N + 1];
    double *cc = new double[N + 1];
    double *dd = new double[N + 1];
    double *c_dash = new double[N + 1];
    double *d_dash = new double[N + 1];
    double *aya = new double[M + 1];
    double *byb = new double[M + 1];
    double *cyc = new double[M + 1];
    double *dyd = new double[M + 1];
    double *cy_dash = new double[M + 1];
    double *dy_dash = new double[M + 1];
    tt[0] = 0.0;
    xx[0] = x_0;
    //col.push_back(100 * x_0);
    for (int i = 1; i <= N; i++)
    {
        xx[i] = xx[i - 1] + h;
        //col.push_back((100 * xx[i]));
    }
    yy[0] = y_0;
    double val;
    for (int i = 1; i <= M; i++)
    {
        yy[i] = yy[i - 1] + dy;
    }
    for (int i = 0; i <= N; i++)
    {
        vector<double> line;
        for (int j = 0; j <= M; j++)
        {
            val = init(xx[i], yy[j]);
            if (val < 1e-7)
                val = 0;
            line.push_back(val);
        }
        instant.push_back(line);
    }
    u.push_back(instant); // u is stored as u[0][x][y] = u(0 , x, y) ; u[1][x][y] = u( dt/2 , x, y) ; u[2][x][y] = u( 2* dt/2 , x, y)  ...
    cout << "at t = 0, \n";
    for (int i = 0; i <= N; i++)
    {
        for (int j = 0; j <= M; j++)
        {
            //cout<<instant[i][j];  // <<"\t"<<"instant ["<<i<<"]["<<j<<"] = "
            printf("u(0 , %0.3lf, %0.3lf) = %0.5lf\t", xx[i], yy[j], u[0][i][j]);
        }
        cout << "\n";
    }
    cout<<"\n";
    for (int i = 1; i <= 2 * N; i++)
    {
        tt[i] = tt[i - 1] + dt;
    }

    for (int n = 1; n < timesteps; n ++)
    {
        //step 1
        // these are the initial conditions that need
        for (int i = 0; i <= N; i++)
        {
            instant[i][0] = 0;
            instant[i][M] = 0;
        }
        for (int j = 0; j <= M; j++)
        {
            instant[0][j] = 0;
            instant[N][j] = 0;
        }
        for (int j = 1; j < M; j++)
        {
            for (int i = 1; i < N; i++)
            {
                aa[i] = 1 / (dx * dx);
                bb[i] = (-2 / dt) - 2 / (dx * dx);
                cc[i] = 1 / (dx * dx);
                dd[i] = -((u[n - 1][i][j + 1] + u[n - 1][i][j - 1]) / (dy * dy)) + 2 * ((1 / dy * dy) - (1 / dt)) * u[n - 1][i][j];
            }
            c_dash[1] = cc[1] / bb[1];
            d_dash[1] = (dd[1] - aa[1] * instant[0][j]) / bb[1];
            for (int i = 2; i < (N - 1); i++)
            {
                c_dash[i] = cc[i] / (bb[i] - aa[i] * c_dash[i - 1]);
                d_dash[i] = (dd[i] - aa[i] * d_dash[i - 1]) / (bb[i] - aa[i] * c_dash[i - 1]);
            }
            instant[N - 1][j] = (dd[N - 1] - cc[N - 1] * instant[N][j] - aa[N - 1] * d_dash[N - 2]) / (bb[N - 1] - aa[N - 1] * c_dash[N - 2]);

            for (int i = N - 2; i >= 1; i--)
            {
                //cout << "c' " << i << " =" << c_dash[i] << " ;d' " << i << " =" << d_dash[i] << endl;
                instant[i][j] = d_dash[i] - c_dash[i] * instant[i + 1][j];
            }
            //u.push_back(col);
        }
        u.push_back(instant);
        cout << "at t = " << n << "*dt/2\n";
        for (int i = 0; i <= N; i++)
        {
            for (int j = 0; j <= M; j++)
            {
                //cout<<instant[i][j];  // <<"\t"<<"instant ["<<i<<"]["<<j<<"] = "
                printf("u(0 , %0.3lf, %0.3lf) = %0.5lf\t", xx[i], yy[j], u[n][i][j]);
            }
            cout << "\n";
        }
        cout << "\n";
        n++;
        //step 2
        // these are the initial conditions that need
        for (int i = 0; i <= N; i++)
        {
            instant[i][0] = 0;
            instant[i][M] = 0;
        }
        for (int j = 0; j <= M; j++)
        {
            instant[0][j] = 0;
            instant[N][j] = 0;
        }
        for (int i = 1; i < N; i++)
        {
            for (int j = 1; j < M; j++)
            {
                aya[j] = 1 / (dy * dy);
                byb[j] = (-2 / dt) - 2 / (dy * dy);
                cyc[j] = 1 / (dy * dy);
                dyd[j] = -((u[n - 1][i + 1][j] + u[n - 1][i-1][j ]) / (dx * dx)) + 2 * ((1 / dx * dx) - (1 / dt)) * u[n - 1][i][j];
            }
            cy_dash[1] = cyc[1] / byb[1];
            dy_dash[1] = (dyd[1] - aya[1] * instant[i][0]) / byb[1];
            for (int j = 2; j < (M - 1); j++)
            {
                cy_dash[j] = cyc[j] / (byb[j] - aya[j] * cy_dash[j - 1]);
                dy_dash[j] = (dyd[j] - aya[j] * dy_dash[j - 1]) / (byb[j] - aya[j] * cy_dash[j - 1]);
            }
            instant[i][M-1] = (dyd[M - 1] - cyc[M - 1] * instant[i][M] - aya[M - 1] * dy_dash[M - 2]) / (byb[M - 1] - aya[M - 1] * cy_dash[M - 2]);

            for (int j = M - 2; j >= 1; j--)
            {
                //cout << "c' " << i << " =" << cy_dash[j] << " ;d' " << i << " =" << dy_dash[j] << endl;
                instant[i][j] = dy_dash[j] - cy_dash[j] * instant[i][j+1];
            }
        }
        u.push_back(instant);
        cout << "at t = " << n << "*dt/2\n";
        for (int i = 0; i <= N; i++)
        {
            for (int j = 0; j <= M; j++)
            {
                //cout<<instant[i][j];  // <<"\t"<<"instant ["<<i<<"]["<<j<<"] = "
                printf("u(0 , %0.3lf, %0.3lf) = %0.5lf\t", xx[i], yy[j], u[n][i][j]);
            }
            cout << "\n";
        }
        cout << "\n";
    }
    /*
    cout << "n -> |";
    for (int j = 0; j < N * 2; j++)
    {
        cout<<"\t"<<j;
        //printf("\t\t %d", j);
    }
    printf("\n");
    cout << "t =  |";
    for (int j = 0; j < N * 2 ; j++)
    {
        cout<<"\t"<<tt[j];
        //printf("\t\t %0.4lf", tt[j]);
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
        for (int i = 0; i < N * 2 ; i++)
        {
            cout<<"\t"<<u[i][j];
            //printf("\t\t%0.4lf", u[i][j]);
        }
        cout << "\n";
    }*/

    return 0;
}