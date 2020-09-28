#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>


#define N 2

// using Freivald’s Algorithm to check if matrix is correct, "1" means correct, "0" means wrong
int check(double a[N][N], double b[N][N], double c[N][N]){
	//double r[N][1];
    double **r;
    r = (double **)malloc(N*sizeof(double*));
    for(int i = 0; i < N; i++){
        r[i] = (double*)malloc(sizeof(double));
    }
    int i,j,k;
	for (i = 0; i < N; i++)
    {
        r[i][0] = rand() % 2;
    }

    //ouble br[N][1];
    double **br;
    br = (double **)malloc(N*sizeof(double*));
    for(int i = 0; i < N; i++){
        br[i] = (double*)malloc(sizeof(double));
    }
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
    //double cr[N][1];
    double **cr;
    cr = (double **)malloc(N*sizeof(double*));
    for(int i = 0; i < N; i++){
        cr[i] = (double*)malloc(sizeof(double));
    }
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
    //double abr[N][1];
    double **abr;
    abr = (double **)malloc(N*sizeof(double*));
    for(int i = 0; i < N; i++){
        abr[i] = (double*)malloc(sizeof(double));
    }
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
        if (abr[i][0] <= 0.000005)
            continue;
        else
            flag = 0;
    }
    return flag;
}

int main(){
    double a[2][2] = {{1.234567,2.234567},{3.234567,4.234567}};
    double b[2][2] = {{1.234567,2.234567},{3.234567,4.234567}};
    double c[2][2] = {{8.75201,12.221146},{17.690280,25.159414}};
    int res = check(a,b,c);
    printf("%d\n",res);
}
