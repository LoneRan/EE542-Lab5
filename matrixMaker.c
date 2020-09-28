#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
double drand(double low, double high)
{   
    
    return ((double)rand()*(high-low))/(double)RAND_MAX+low;
}

int main(){
    int size;
    printf("Please enter matrix size\n");
    scanf("%d",&size);

    double **matrix_A;
    matrix_A = (double **)malloc(size*sizeof(double*));
    for(int i = 0; i < size; i++){
        matrix_A[i] = (double*)malloc(size*sizeof(double));
    }
    double **matrix_B;
    matrix_B = (double **)malloc(size*sizeof(double*));
    for(int i = 0; i < size; i++){
        matrix_B[i] = (double*)malloc(size*sizeof(double));
    }
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            matrix_A[i][j] = drand(1,10);
        }
    }
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            matrix_B[i][j] = drand(1,10);
        }
    }
    //matrix A
    FILE *fp;

    if((fp=fopen("matrixA.txt", "wb"))==NULL) {
        printf("Cannot open file.\n");
    }
    
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size-1; j++){
            fprintf(fp,"%f ",matrix_A[i][j]);
        }
        fprintf(fp,"%f",matrix_A[i][size-1]);
        if(i != size-1) fprintf(fp,"\n");
    }
    fclose(fp);

    //matrix B
    FILE *fp1; 

    if((fp1=fopen("matrixB.txt", "wb"))==NULL) {
        printf("Cannot open file.\n");
    }
    
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size-1; j++){
            fprintf(fp1,"%f ",matrix_B[i][j]);
        }
        fprintf(fp1,"%f",matrix_B[i][size-1]);
        if(i != size-1) fprintf(fp1,"\n");
    }
    fclose(fp1);

    return 0;

}