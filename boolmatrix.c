/*
    Relaxacion method for temperature interpolation
		Boolean matrix generator
        Vivas,A;Arango,C;Arguelles,A; 
        Lascilab
        CIBioFi-QuanTIC
          
        Last updated: January 23, 2017

    This program is in development. 

    Laboratory of Distributed systems and Networks
    http://Lascilab.univalle.edu.co/

*/
#include <stdlib.h> /* Standard Libary: malloc, calloc, free, ralloc functions */
#include <stdio.h>  /* Standard I/O Library: printf */
#include <math.h>   /* Standard Math Library */

void prnt(int * matrix, int nx, int ny){
    	int i;
    	for (i = 1; i <= nx * ny; i++)
    		{
        int b_i = (i-1)/nx, b_j = (i-1)%nx; 
        printf("%d\t",matrix[i-1]);
        if(i % nx == 0) printf("\n");
    		}
	}

void fileload(int* matrix,int Nx,int Ny){
   	int file;
   	char FileName[50];
   	for (int i = 0; i < Nx * Ny; ++i) matrix[i] = 0;
   	for(int f=2; f <=3; f++)	
    		{
    	file=sprintf(FileName,"Maps/Contour_VALLE_960_Border_%d.dat",f);
    	file++;
    	FILE* file = fopen(FileName, "r");
    	int x = 0, y = 0;
    	for(; fscanf(file, "%d", &x) && fscanf(file, "%d", &y) && !feof(file);) matrix [ y * Nx + x ] = 1;
 	fclose(file);
		}
}

int main(int argc, char const **argv){

    	int Nx = 500;
    	int Ny = 1397;
       	int * matrix = (int*) malloc( Nx * Ny * sizeof(int));
    
    	fileload(matrix,Nx,Ny);
    	prnt(matrix,Nx,Ny);
    	
    	free(matrix);
    
    	return 0;
}
