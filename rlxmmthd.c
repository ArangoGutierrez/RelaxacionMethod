/*
	Relaxacion method for temperature interpolation

		Vivas,A;Arango,C;Arguelles,A; 
		Lascilab
		CIBioFi-QuanTIC
		  
		Last updated: January 13, 2017

	This program is in development. 

	Laboratory of Distributed systems and Networks
	https://http://Lascilab.univalle.edu.co/

*/
#include <stdlib.h>	/* Standard Libary: malloc, calloc, free, ralloc functions */
#include <stdio.h> 	/* Standard I/O Library: printf */

struct Station
{
	int x;		/* x coordinate */
	int y;		/* y coordinate */
	double t;	/* temperature */
};

void fillB(int * matrix, int nx, int ny, struct Station *s, int ne){
	int i = 0, row = 0, col = 0;
	for (i = 0; i < nx * ny; i++){
		row = i / nx, col = i % nx;
		matrix[i] = ( row == 0 || col == 0  || row == ny - 1 || col == nx - 1) ? 0 : 1;
	}

	int s_i = 0, s_j = 0;
	for (i = 0; i < ne; ++i){
		s_i = s[i].x * nx, s_j = s[i].y;
		matrix[ s_i + s_j ] = 2;
	}
}

void fillT(double * matrix, int nx, int ny,struct Station *s, int ne){
	double sum = 0;
	int i;
	for (i = 0; i < ne; i++) sum = sum + s[i].t;
	double To =  ( sum / ne );
	
	for (i=0; i < nx * ny; i++) matrix[i] = To;

	int s_i = 0, s_j = 0;
	for (i=0; i < ne ; i++) {
		s_i = s[i].x * nx;
		s_j = s[i].y;
		matrix[( s_i + s_j ) ] = s[i].t;
	}
}

void prntB(int * matrix, int nx, int ny){
	int i;
	for (i = 1; i <= nx * ny; i++)
	{
		printf("%d\t", matrix[i-1]);
		if(i % nx == 0) printf("\n");
	}
}

void prntT(double * matrix, int nx, int ny){
	int i;
	for (i = 1; i <= nx * ny; i++)
	{
		printf("%f\t", matrix[i-1]);
		if(i % nx == 0) printf("\n");
	}
}

double transitionFunc(int * B, double * T, int nx, int ny, int cell){
	int row = 0, col = 0;
	row = cell / nx, col = cell % nx;
	double up = ( row - 1 < 0 | B[ cell - nx ] == 0 ) ? 0 : T[cell - nx];
	double down = ( row + 1 >= ny | B[ cell + nx ] == 0 ) ? 0 : T[cell + nx];
	double left = ( col - 1 < 0 | B[ cell - 1 ] == 0 ) ? 0 : T[ cell - 1 ];
	double right = ( col + 1 >= nx | B[ cell + 1 ] == 0) ? 0 : T[ cell + 1];
	double sum = ((up > 0) ? 1:0) + ((down > 0) ? 1:0) + ((left > 0) ? 1:0) + ((right > 0) ? 1:0);
	double Tij = (up + down + left + right) / sum;
	printf("Row %d, Col %d, up %f,down %f, left %f, right %f, sum %f, Tij %f\n",
		row,col,up,down,left,right,sum,Tij);
	return Tij;
}

int main(int argc, char const **argv)
{
	int Nx = 10;
	int Ny = 10;
	int Ne = 4;
	
	struct Station s[4];
	//First Station
	s[0].x=2;
	s[0].y=2;
	s[0].t=26.5;

	//Second Station
	s[1].x=7;
	s[1].y=2;
	s[1].t=29.3;

	//Third Station
	s[2].x=2;
	s[2].y=7;
	s[2].t=28.7;

	//Fourth Station
	s[3].x=7;
	s[3].y=7;
	s[3].t=30.1;
	
	int * B = NULL;
	double * Ta = NULL;
	double * Tb = NULL;

	B = (int *) malloc( Nx * Ny * sizeof(int));
	Ta = (double *) malloc(Nx * Ny * sizeof(double));
	Tb = (double *) malloc(Nx * Ny * sizeof(double));

	fillB(B,Nx,Ny,s,Ne);
	fillT(Ta,Nx,Ny,s,Ne);
	fillT(Tb,Nx,Ny,s,Ne);
	printf("B\n");
	prntB(B,Nx,Ny);
	printf("Ta\n");
	prntT(Ta,Nx,Ny);
	printf("Tb\n");
	prntT(Tb,Nx,Ny);
	printf("transitionFunc\n");
	double result = transitionFunc(B,Ta,Nx,Ny,32);
	
	free(B);
	free(Ta);
	free(Tb);


	return 0;
}
