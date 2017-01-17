/*
	Relaxacion method for temperature interpolation

		Vivas,A;Arango,C;Arguelles,A; 
		Lascilab
		CIBioFi-QuanTIC
		  
		Last updated: January 17, 2017

	This program is in development. 

	Laboratory of Distributed systems and Networks
	http://Lascilab.univalle.edu.co/

*/
#include <stdlib.h>	/* Standard Libary: malloc, calloc, free, ralloc functions */
#include <stdio.h> 	/* Standard I/O Library: printf */
#include <math.h>	/* Standard Math Library */

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
	for (i = 0; i < ne; ++i) matrix[ s[i].y * nx + s[i].y ] = 2;
}

void fillT(double * matrix, int nx, int ny,struct Station *s, int ne){
	double sum = 0;
	int i;
	for (i = 0; i < ne; i++) sum = sum + s[i].t;
	double To =  ( sum / ne );
	
	for (i=0; i < nx * ny; i++) matrix[i] = To;
	for (i=0; i < ne ; i++) matrix[( s[i].y * nx + s[i].x ) ] = s[i].t;
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

int test(int * B, double * Ta, double * Tb, int nx, int ny){
	int i = 0;
	double result = 0;
	for (i = 0; i < nx * ny; ++i) result += (B[i] == 0 || B[i] == 2) ? 0 : pow(Tb[i] - Ta[i],2);
	printf("Test Result delta: %.10f\n",result);
	//result = (1 / (nx*ny)) * sqrt(result);
	result =  sqrt(result);
	//printf("Test Result: %.10f\n",result);
	return ( result <= 0.0000001 ) ? 1 : 0;
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
	int Nx = 100;
	int Ny = 90;
	int Ne = 10;
	
	struct Station s[Ne];
	//First Station
	s[0].x=65;
	s[0].y=10;
	s[0].t=22.6;

	//Second Station
	s[1].x=55;
	s[1].y=20;
	s[1].t=23.2;

	//Third Station
	s[2].x=50;
	s[2].y=25;
	s[2].t=22.7;

	//Fourth Station
	s[3].x=50;
	s[3].y=35;
	s[3].t=22.9;

	//Fith Station
	s[4].x=40;
	s[4].y=45;
	s[4].t=22.0;

	//Sixth Station
	s[5].x=50;
	s[5].y=60;
	s[5].t=23.7;

	//Seventh Station
	s[6].x=35;
	s[6].y=65;
	s[6].t=23.1;

	//Eighth Station
	s[7].x=50;
	s[7].y=70;
	s[7].t=23.1;

	//Ninth Station
	s[8].x=35;
	s[8].y=75;
	s[8].t=22.9;

	//Tenth Station
	s[9].x=45;
	s[9].y=80;
	s[9].t=23.5;
	
	int * B = NULL;
	double * Ta = NULL;
	double * Tb = NULL;

	B = (int *) malloc( Nx * Ny * sizeof(int));
	Ta = (double *) malloc(Nx * Ny * sizeof(double));
	Tb = (double *) malloc(Nx * Ny * sizeof(double));

	fillB(B,Nx,Ny,s,Ne);
	fillT(Ta,Nx,Ny,s,Ne);
	fillT(Tb,Nx,Ny,s,Ne);
	//printMatrixes(B,Ta,Tb,Nx,Ny);

	for (int i = 0 ; i < 10; i++) {

	evolve(B,Ta,Tb,Nx,Ny);
	
	FILE *dskw1;
	char FileName[50];
	int file;
	file=sprintf(FileName,"./OutputData/RlxMthd_v1.0_%d.dat",i);
	file++;
  	dskw1=fopen(FileName,"w+");

	for (int c = 0; c < Nx * Ny; c++)
	{		
	int x=( c % Nx ) + 1;
	int y=( c / Nx ) + 1;
	fprintf(dskw1,"%d\t%d\t%f\n",x,y,Tb[c]);
	
	}

	double * temp = Ta;
	Ta = Tb;
	Tb = temp;

	}
	
	
	/*int i = 0;
	for (i = 0; i < 20000; ++i)
	{
		printf("Generation %d\n",i+1);
		evolve(B,Ta,Tb,Nx,Ny);
		double * temp = Ta;
		Ta = Tb;
		Tb = temp;
		printf("Test: %d\n",test(B,Ta,Tb,Nx,Ny));
	}*/

	//printMatrixes(B,Ta,Tb,Nx,Ny);

	/*
	printf("Transition: %f\n",transition(B,Ta,Nx,Ny,32));

	printf("Evolve\n");
	
	printf("Generation 0\n");
	printf("Test: %d\n",test(B,Ta,Tb,Nx,Ny));
	printMatrixes(B,Ta,Tb,Nx,Ny);
	
	printf("Generation 1\n");
	evolve(B,Ta,Tb,Nx,Ny);
	double * temp = Ta;
	Ta = Tb;
	Tb = temp;
	printf("Test: %d\n",test(B,Ta,Tb,Nx,Ny));
	printMatrixes(B,Ta,Tb,Nx,Ny);

	printf("Generation 2\n");
	evolve(B,Ta,Tb,Nx,Ny);
	temp = Ta;
	Ta = Tb;
	Tb = temp;
	printf("Test: %d\n",test(B,Ta,Tb,Nx,Ny));
	printMatrixes(B,Ta,Tb,Nx,Ny);
	*/

	free(B);
	free(Ta);
	free(Tb);

	return 0;
}
