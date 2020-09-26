#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>


#define N 14000

// using Freivald’s Algorithm to check if matrix is correct, "1" means correct, "0" means wrong
int check(double matrixA[N][N], double matrixB[N][N], double matrixC[N][N]){
	double r[N][1];
	for (i = 0; i < N; i++)
    {
        r[i][0] = rand() % 2;
    }

    double br[N][1];
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < 1; j++)
            {
                for (k = 0; k < N; k++)
                    {
                        br[i][j] = br[i][j] + b[i][k] * r[k][j];
                    }
            }
    }
    double cr[N][1];
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < 1; j++)
            {
                for (k = 0; k < N; k++)
                    {
                        cr[i][j] = cr[i][j] + c[i][k] * r[k][j];
                    }
            }
    }
    double abr[N][1];
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < 1; j++)
            {
                for (k = 0; k < N; k++)
                    {
                        abr[i][j] = abr[i][j] + a[i][k] * br[k][j];
                    }
            }
    }
    for (i = 0; i < N; i++)
    {
        abr[i][0] -= cr[i][0];
    }
    int flag = 1;
    for (i = 0; i < N; i++)
    {
        if (abr[i][0] == 0)
            continue;
        else
            flag = 0;
    }
    return flag;
}