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

double transition(int * B, double * T, int nx, int ny, int cell){
	int row = 0, col = 0;
	row = cell / nx, col = cell % nx;

	double up = ( row - 1 < 0 || B[ cell - nx ] == 0 ) ? 0 : T[cell - nx];
	double down = ( row + 1 >= ny || B[ cell + nx ] == 0 ) ? 0 : T[cell + nx];
	double left = ( col - 1 < 0 || B[ cell - 1 ] == 0 ) ? 0 : T[ cell - 1 ];
	double right = ( col + 1 >= nx || B[ cell + 1 ] == 0) ? 0 : T[ cell + 1];
	double sum = ((up > 0) ? 1:0) + ((down > 0) ? 1:0) + ((left > 0) ? 1:0) + ((right > 0) ? 1:0);
	double Tij = (up + down + left + right) / sum;
	//printf("Row %d, Col %d, up %f,down %f, left %f, right %f, sum %f, Tij %f\n",row,col,up,down,left,right,sum,Tij);
	return (B[cell] == 0 || B[cell] == 2) ? T[cell] : Tij;
}

int test(double * Ta, double * Tb, int nx, int ny){
	int i = 0;
	double result = 0;
	for (i = 0; i < nx * ny; ++i){
		result = (1 / (nx*ny)) * (( (Tb[i] - Ta[i]) * (Tb[i] - Ta[i]) ) * 0.5);
		if( result > 0.1 ) return 0;
	}
	return 1;
}

void evolve(int * B, double * Tin, double * Tout, int nx, int ny){
	int i = 0;
	for (i = 0; i < nx * ny; ++i) Tout[i] = transition(B,Tin,nx,ny,i);
}

void printMatrixes(int * B, double * Ta, double * Tb, int nx, int ny){
	printf("B\n");
	prntB(B,nx,ny);
	printf("Ta\n");
	prntT(Ta,nx,ny);
	printf("Tb\n");
	prntT(Tb,nx,ny);
}

int main(int argc, char const **argv)
{
	int Nx = 10;
	int Ny = 15;
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
	printMatrixes(B,Ta,Tb,Nx,Ny);

	printf("Transition: %f\n",transition(B,Ta,Nx,Ny,32));

	printf("Evolve\n");
	
	printf("Generation 0\n");
	printf("Test: %d\n",test(Ta,Tb,Nx,Ny));
	printMatrixes(B,Ta,Tb,Nx,Ny);
	
	printf("Generation 1\n");
	evolve(B,Ta,Tb,Nx,Ny);
	double * temp = Ta;
	Ta = Tb;
	Tb = temp;
	printf("Test: %d\n",test(Ta,Tb,Nx,Ny));
	printMatrixes(B,Ta,Tb,Nx,Ny);

	printf("Generation 2\n");
	evolve(B,Ta,Tb,Nx,Ny);
	temp = Ta;
	Ta = Tb;
	Tb = temp;
	printf("Test: %d\n",test(Ta,Tb,Nx,Ny));
	printMatrixes(B,Ta,Tb,Nx,Ny);

	free(B);
	free(Ta);
	free(Tb);

	return 0;
}
