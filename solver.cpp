#include <stdio.h>
#include <math.h>

#define N 9 // row count in the overdetermined matrix
#define M 4 // x y z coordinates and time

int main()
{
    double a[N][M];
    double at[M][N];
    double b[N];
    double P[M][M + 1];

    // ground station coordinates
    double x[N + 1] = {0.0, 0.5, 0.7, 0.1, 0.5, 1.5, 2.0, 0.7, 0.8, 2.5};
    double y[N + 1] = {0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.5, 0.7, 0.2, 2.0};
    double z[N + 1] = {0.0, 0.1, 0.1, 0.5, 0.5, 1.5, 1.0, 1.3, 0.7, 2.0};
    
    // distances from the stations to signal origin
    double r[N + 1];
    
    // real signal coordinates
    double coord[3] = {2, 2, 5};
    
    // calculated signal coordinates
    double sol[M];

    // calculating distances from the stations
    for (int i = 0; i < N + 1; i++)
    {
    	r[i] = sqrt(pow(x[i] - coord[0], 2) + pow(y[i] - coord[1], 2) + pow(z[i] - coord[2], 2));
    }  
    
    // linear system of equations are obtained with method described in
    // http://www.diku.dk/~rfonseca/implementations/apollonius3d.pdf
    
    for (int i = 0; i < N; i++)
    {
        a[i][0] = 2 * (x[i + 1] - x[0]);
        a[i][1] = 2 * (y[i + 1] - y[0]);
        a[i][2] = 2 * (z[i + 1] - z[0]);
        a[i][3] = 2 * (r[i + 1] - r[0]);
        b[i] = (pow(x[i + 1], 2) + pow(y[i + 1], 2) + pow(z[i + 1], 2))
                  - (pow(x[0], 2) + pow(y[0], 2) + pow(z[0], 2))
                  + pow(r[0], 2) - pow(r[i + 1], 2);
    }

    // solving system a * sol = b with least-square method 

    // transpose matrix
    for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
        {
            at[i][j] = a[j][i];
        }


    // at * a * sol = at * b
    // P - extended matrix
    // P = ( at * a | at * b ) 
    
    for (int i = 0; i < M; i++)
        for (int j = 0; j < M; j++)
        {
            double sum = 0;

            for (int k = 0; k < N; k++)
            {
                sum += at[i][k] * a[k][j];
            }

            P[i][j] = sum;
        }

    for (int i = 0; i < M; i++)
    {
        double sum = 0;

        for (int k = 0; k < N; k++)
        {
            sum += at[i][k] * b[k];
        }

        P[i][M] = sum;
    }

    // Gaussian elimination

    for (int i = 0; i < M; i++)
    {
        double div = P[i][i];
        if (div != 0)
        {
            for (int j = i; j <= M; j++)
            {
                P[i][j] /= div;
            }
        }
        
        for (int j = i + 1; j < M; j++)
        {
            double mult = P[j][i];
            for (int k = i; k <= M; k++)
            {
                P[j][k] -= P[i][k] * mult;
            }
        }
    }
    
    for (int i = M - 1; i >= 0; i--)
    {
        double sum = 0;
        for (int j = i + 1; j < M; j++)
        {
            sum += sol[j] * P[i][j];
        }
        sol[i] = P[i][M] - sum;
    }

    printf("Solution:\n\rx = %f y = %f z = %f r = %f \n\r", sol[0], sol[1], sol[2], sol[3]);

}
