#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include"mpi.h"

#define N 4 /* dimension of the input matrix */

double matrixA[N][N];
double matrixB[N][N];
double matrixC[N][N];
double temp0[N][N];
double temp1[N][N];
double temp2[N][N];
int count = 0;
typedef struct {
	int p; /* number of processors */
	MPI_Comm comm; /* handle to global grid communicator */
	MPI_Comm row_comm; /* row communicator */
	MPI_Comm col_comm; /* column communicator */
	int q; /* dimension of the grid, = sqrt(p) */
	int my_row; /* row position of a processor in a grid */
	int my_col; /* column position of a procesor in a grid */
	int my_rank; /* rank within the grid */
}GridInfo;

double[N][N] addMatrix(double m1[N][N], m2[N][N]){
	for(int i = 0; i < N;  i++){
		for(int j = 0; j < N; j++){
			m1[i][j] += m2[i][j];	
		}
	}
	return m1;
}

int check(double matrixA[N][N], double matrixB[N][N], double matrixC[N][N]){
	double r[N][1];
	int i, j, k;
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
                        br[i][j] = br[i][j] + matrixB[i][k] * r[k][j];
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
                        cr[i][j] = cr[i][j] + matrixC[i][k] * r[k][j];
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
                        abr[i][j] = abr[i][j] + matrixA[i][k] * br[k][j];
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

void SetupGrid(GridInfo *grid)
{
	int old_rank;
	int dimensions[2];
	int wrap_around[2];
	int coordinates[2];
	int free_coords[2];
	
	/* get the overall information before overlaying cart_grid */

	MPI_Comm_size(MPI_COMM_WORLD,&(grid->p));
	MPI_Comm_rank(MPI_COMM_WORLD,&old_rank);
	
	/* Assumption: p is a perfect square */
	grid->q=(int)sqrt((double)grid->p);
	/* set the dimensions */
	dimensions[0]=dimensions[1]=grid->q;
	
	/* we want a torus on the second dimension, so set it appropriately */

	wrap_around[0]=0;
	wrap_around[1]=1;
	
	MPI_Cart_create(MPI_COMM_WORLD,2,dimensions,wrap_around,1,&(grid->comm));
	/* since we have set reorder to true, this might have changed the ranks */
	MPI_Comm_rank(grid->comm,&(grid->my_rank));
	/* get the cartesian coordinates for the current process */
	MPI_Cart_coords(grid->comm,grid->my_rank,2,coordinates);
	/* set the coordinate values for the current coordinate */
	grid->my_row=coordinates[0];
	grid->my_col=coordinates[1];

        /* create row communicators */
	free_coords[0]=0;
	free_coords[1]=1; /* row is gonna vary */
	MPI_Cart_sub(grid->comm,free_coords,&(grid->row_comm));
	
        /* create column communicators */
	free_coords[0]=1;
	free_coords[1]=0; /* row is gonna vary */
	MPI_Cart_sub(grid->comm,free_coords,&(grid->col_comm));
	
}

/* normal matrix multiplication stuff */

void matmul(double **a, double **b, double **c, int size)
{
	int i,j,k;
       
	double **temp = (double**) malloc(size*sizeof(double*));
	for(i=0;i<size;i++)
		*(temp+i)=(double*) malloc(size*sizeof(double));

	for(i=0;i<size;i++)
	{
			for(j=0;j<size;j++)
			{
				temp[i][j]=0;
				for(k=0;k<size;k++){
					temp[i][j]=temp[i][j]+ (a[i][k] * b[k][j]);
				}
			}
	}
	
	for(i=0;i<size;i++)
		for(j=0;j<size;j++)
			c[i][j]+=temp[i][j];
	
}
void transfer_data_from_buff(double *buff,double **a,int buffsize, int row, int col){
	
  	if(buffsize!=row*col)
	{
		printf("transfer_data_from_buf: buffer size does not match matrix size!\n");
		exit(1);
	}
	int count=0, i,j;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++){
			a[i][j]=buff[count];
			count++;
		}
	}
}

void transfer_data_to_buff(double *buff,double **a,int buffsize, int row, int col){
	
  	if(buffsize!=row*col)
	{
		printf("transfer_data_to_buf: buffer size does not match matrix size!");
		exit(1);
	}
	int count=0, i,j;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++){
			buff[count]=a[i][j];
			count++;
		}
	}
}

void Fox(int n,GridInfo *grid,double **a, double **b, double **c)
{
	double **tempa;
	double *buff; /* buffer for Bcast & send_recv */
	int stage;
	int root;
	int submat_dim; /* = n/q */
	int source;
	int dest;
	int i;
	MPI_Status status;
	
	submat_dim=n/grid->q;
	
	/* Initialize tempa */
	tempa=(double**) malloc(submat_dim*sizeof(double*));
	for(i=0;i<submat_dim;i++)
		*(tempa+i)=(double*) malloc(submat_dim*sizeof(double));
	/* initialize buffer */
	buff=(double*)malloc(submat_dim*submat_dim*sizeof(double));

        /* we are gonna shift the elements of matrix b upwards with the column fixed */
	source = (grid->my_row+1) % grid->q; /* pick the emmediately lower element */
	dest= (grid->my_row+grid->q-1) % grid->q; /* move current element to immediately upper row */
	
	
	for(stage=0;stage<grid->q;stage++)
	{
		root=(grid->my_col+stage)%grid->q;
		if(root==grid->my_col)
		{
			transfer_data_to_buff(buff,a,submat_dim*submat_dim, submat_dim,submat_dim);
			MPI_Bcast(buff,submat_dim*submat_dim,MPI_INT,root,grid->row_comm);
			transfer_data_from_buff(buff,a,submat_dim*submat_dim, submat_dim,submat_dim);
		
			matmul(a,b,c,submat_dim);
		}else
		{
			transfer_data_to_buff(buff,tempa,submat_dim*submat_dim, submat_dim,submat_dim);
			MPI_Bcast(buff,submat_dim*submat_dim,MPI_INT,root,grid->row_comm);
			transfer_data_from_buff(buff,tempa,submat_dim*submat_dim, submat_dim,submat_dim);
			
			matmul(tempa,b,c, submat_dim);
		}
		transfer_data_to_buff(buff,b,submat_dim*submat_dim, submat_dim,submat_dim);
		MPI_Sendrecv_replace(buff,submat_dim*submat_dim,MPI_INT,dest,0,source,0,grid->col_comm,&status);
		transfer_data_from_buff(buff,b,submat_dim*submat_dim, submat_dim,submat_dim);
	}

}

void initialiseAB()
{
	int i,j;
/* *****************************************************************************************
Initialize the input matrix 
Note: This initalization is deterministic & hence is done by every process in the same way 
      I wanted to design a fully distributed program, hence I took this strategy.
      A better strategy could have been to let master alone initialize the matrices & then
      send the slaves their local copy only 
*******************************************************************************************/
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			matrixA[i][j]=rand() % 10 + 1;
			matrixB[i][j]=rand() % 10 + 1;
		}
	}
	
			
	
}


int main(int argc, char *argv[])
{	
	
	int i,j,dim;
	double **localA;
	double **localB;
	double **localC;
	MPI_Init (&argc, &argv);
	
	GridInfo grid;
	/*initialize Grid */

	SetupGrid(&grid);
        /* Initialize matrix A & B */
	initialiseAB();
        /* calculate local matrix dimension */
	dim=N/grid.q;
	/* allocate space for the three matrices */		

	
	localA=(double**) malloc(dim*sizeof(double*));

	localB=(double**) malloc(dim*sizeof(double*));
	
	localC=(double**) malloc(dim*sizeof(double*));
	
	for(i=0;i<dim;i++)
	{
		*(localA+i)=(double*) malloc(dim*sizeof(double));
		*(localB+i)=(double*) malloc(dim*sizeof(double));
		*(localC+i)=(double*) malloc(dim*sizeof(double));
	}


/* Compute local matrices - Ideally the master should do this & pass it onto all the slaves */
/* At the same time initialize localC to all zeros */

	int base_row=grid.my_row*dim;
	int base_col=grid.my_col*dim;

	for(i=base_row;i<base_row+dim;i++)
	{
		for(j=base_col;j<base_col+dim;j++)
		{
		         localA[i-(base_row)][j-(base_col)]=matrixA[i][j];
			 localB[i-(base_row)][j-(base_col)]=matrixB[i][j];
			 localC[i-(base_row)][j-(base_col)]=0;
		}
	}



	Fox(N,&grid,localA, localB, localC);

/* print results */
	printf("rank=%d, row=%d col=%d\n",grid.my_rank,grid.my_row,grid.my_col);
	for(i=0;i<dim;i++)
	{
		for(j=0;j<dim;j++)
		{
			//printf("localC[%d][%d]=%d ", i,j,localC[i][j]);
			printf("%.2f ", localC[i][j]);
		}
		printf("\n");
	}
	
	
	for(i=base_row;i<base_row+dim;i++){
		for(j=base_col;j<base_col+dim;j++)
		{
			matrixC[i][j] = localC[i-(base_row)][j-(base_col)];
		}
	}
	if(grid.my_rank != 3)
	{
		MPI_Send(&matrixC([0][0]), N*N, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD);
	}
	if(grid.my_rank == 3)
	{
		MPI_Recv(&temp0([0][0]), N*N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(&temp1([0][0]), N*N, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv(&temp2([0][0]), N*N, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		matrixC = addMatrix(matrixC, temp0);
		matrixC = addMatrix(matrixC, temp1);
		matrixC = addMatrix(matrixC, temp2);
		
	}
	
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			//printf("localC[%d][%d]=%d ", i,j,localC[i][j]);
			printf("%.2f ", matrixC[i][j]);
		}
		printf("\n");
	}
	
	
	if(check(matrixA,matrixB,matrixC)==1) printf("result is correct!\n");
	else printf("result is wrong!\n");
	
		
	MPI_Finalize ();
	exit(0);

}
