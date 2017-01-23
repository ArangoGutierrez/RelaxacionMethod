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

int main(int argc, char const **argv){

    int Nx = 500;
    int Ny = 1397;

    int * matrix = (int*) malloc( Nx * Ny * sizeof(int));
    int i = 0;
    for (i = 0; i < Nx * Ny; ++i) matrix[i] = 0;


    FILE* file = fopen("Contour_VALLE_960_Border_2.dat", "r");

    int x = 0, y = 0;
    // https://www.quora.com/How-do-you-read-integers-from-a-file-in-C
    for(; fscanf(file, "%d", &x) && fscanf(file, "%d", &y) && !feof(file);) matrix [ y * Nx + x ] = 1;

    prnt(matrix,Nx,Ny);

    fclose(file);
    free(matrix);
    
    return 0;
}