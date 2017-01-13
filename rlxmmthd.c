#include <stdlib.h>	/* Standard Libary: malloc, calloc, free, ralloc functions */
#include <stdio.h> 	/* Standard I/O Library: printf */

struct Station
{
	int x;		/* x coordinate */
	int y;		/* y coordinate */
	double t;	/* temperature */
};

void fillB(int * matrix, int nx, int ny){
	int i = 0, row = 0, col = 0;
	for (i = 0; i < nx * ny; ++i){
		row = i / nx;
		col = i % nx;
		matrix[i] = ( row == 0 || col == 0  || row == ny - 1 || col == nx - 1) ? 0 : 1;
	}
		

}

void fillT(double * matrix, int dx, int dy, struct Station * s, int ne){
	int sum = 0;
	int i = 0;
	for (i = 0; i < ne; ++i) sum += s[i].t;
	double To = (1 / ne) * sum;
	
	for (i=0; i < dx * dy; ++i) matrix[i] = To;
	for (i = 0; i < ne ; ++i) matrix[ s[i].x + s[i].y ] = s[i].t;
}

void show(int * matrix, int dx, int dy){
	int i = 0;
	for (i = 1; i <= dx * dy; ++i)
	{
		printf("%d ", matrix[i-1]);
		if(i % dy == 0) printf("\n");
	}
}

int main(int argc, char const **argv)
{
	int Nx = 10;
	int Ny = 10;

	int * B = NULL;
	double * Ta = NULL;
	double * Tb = NULL;

	B = (int *) malloc( Nx * Ny * sizeof(int));
	Ta = (double *) malloc(Nx * Ny * sizeof(double));
	Tb = (double *) malloc(Nx * Ny * sizeof(double));

	fillB(B,Nx,Ny);
	show(B,Nx,Ny);



	free(B);
	free(Ta);
	free(Tb);


	return 0;
}
