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

void fillB(int * matrix, int nx, int ny){
	int i = 0, row = 0, col = 0;
	for (i = 0; i < nx * ny; i++){
		row = i / nx;
		col = i % nx;
		matrix[i] = ( row == 0 || col == 0  || row == ny - 1 || col == nx - 1) ? 0 : 1;
	}
		

}

void fillT(double * matrix, int dx, int dy, struct Station *s, int ne){
	int sum = 0;
	int i;
	for (i = 0; i < ne; i++) {sum = sum + s[i].t;};
	double To =  ( sum / ne );
	
	for (i=0; i < dx * dy; i++) {matrix[i] = To;};
	for (i=0; i < ne ; i++) {matrix[( s[i].x * s[i].y) ] = s[i].t;};
	printf("To=%f\nS[0]=%f\nS[1]=%f\nS[2]=%f\nS[3]=%f\n",To,s[0],s[1],s[2],s[3]);
}

void prnt(double * matrix, int dx, int dy){
	int i;
	for (i = 1; i <= dx * dy; i++)
	{
		printf("%f\t", matrix[i-1]);
		if(i % dx == 0) printf("\n");
	}
}

int main(int argc, char const **argv)
{
	int Nx = 10;
	int Ny = 10;
	
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

	fillB(B,Nx,Ny);
	fillT(Ta,Nx,Ny,s,4);
	prnt(Ta,Nx,Ny);



	free(B);
	free(Ta);
	free(Tb);


	return 0;
}
