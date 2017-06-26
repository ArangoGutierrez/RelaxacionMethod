/*
    Relaxacion method for temperature interpolation
        Boolean matrix generator
        Vivas,A;Arango,C;Arguelles,A; 
        Lascilab
        CIBioFi-QuanTIC
          
        Last updated: June 26, 2017

    This program is in development. 

    Laboratory of Distributed systems and Networks
    http://Lascilab.univalle.edu.co/

*/
#include <stdlib.h> /* Standard Libary: malloc, calloc, free, ralloc functions */
#include <stdio.h>  /* Standard I/O Library: printf */
#include <math.h>   /* Standard Math Library */
#include <time.h>   /* to take execution time */

void mapload(int* matrix,int Nx,int Ny){
    int file;
    char FileName[50];
    for(int f=1; f <=3; f++){
        file=sprintf(FileName,"../Maps/Contour_VALLE_960_Border_%d.dat",f);
        file++;
        FILE* file = fopen(FileName, "r");
        int x = 0, y = 0;
        for(; fscanf(file, "%d\t%d", &x, &y) && !feof(file);) matrix [ y * Nx + x ] = 2;
        fclose(file);
    }
}

void estload(int* matrix,int Nx,int Ny){
    // FILE* file = fopen("../Maps/estaciones.dat", "r");
    int x = 0, y = 0;
    double t = 0;
    matrix[(1400 - 481) * Nx + 297] = 1;
    matrix[175 * Nx + 54] = 1;
    // for(; fscanf(file, "%d\t%d\t%f", &x, &y, &t) && !feof(file);) matrix [ y * Nx + x ] = 1;
    // fclose(file);       
}

int propagation(int* M,int Nx,int Ny, int i){

    int row = i / Nx;
    int col = i % Nx;
    
    double up =     ( row - 1 < 0   || M[i - Nx]  == 0  || M[i - Nx] == 2 ) ? 0 : 1;    
    double down =   ( row + 1 >= Ny || M[i + Nx]  == 0  || M[i + Nx] == 2 ) ? 0 : 1;    
    double left =   ( col - 1 < 0   || M[i - 1]  == 0   || M[i - 1] == 2 )  ? 0 : 1;    
    double right =  ( col + 1 >= Nx || M[i + 1]  == 0   || M[i + 1] == 2 )  ? 0 : 1;    
    int sum = (up + down + left + right);
    int B = ( sum < 1 ) ? 0 : 1;
    
    return ((M[i]==2) ? 2 : B); 
    }

void saveAsMatrix(int * matrix, int nx, int ny){
    FILE * pf = fopen("../Maps/boolmap_matrix.dat","w");
    int i;
    for (i = 1; i <= nx * ny; i++){
        fprintf(pf,"%d\t",matrix[i-1]);
        if(i % nx == 0) fprintf(pf,"\n");
    }
    fclose(pf);
}

void saveAsTable(int * matrix, int nx, int ny){
    FILE * pf = fopen("../Maps/boolmap_table.dat","w");
    int x = 0, y = 0;
    int i;
    for (i = 1; i <= nx * ny; i++){
        // int b_i = (i-1)/nx, b_j = (i-1)%nx; /* Creeo que hay un error */
        // fprintf(pf,"%d\t%d\t%d\n",b_i,b_j,(matrix[i-1]==2) ? 0 :matrix[i-1] );
        x = ( (i-1) % nx );
        y = ( (i-1) / nx );
        fprintf(pf,"%d\t%d\t%d\n",x,y,(matrix[i-1]==2) ? 0 :matrix[i-1] );
    }
    fclose(pf);
}

void savetime(double time){
    FILE * pf = fopen("time_bool.txt","a+");
    fprintf(pf, "%f\n",time);
    fclose(pf); 
}

void printAsMatrix(int * matrix, int nx, int ny){
    int i;
    for (i = 1; i <= nx * ny; i++){
        int b_i = (i-1)/nx, b_j = (i-1)%nx; 
        printf("%d\t",matrix[i-1]);
        if(i % nx == 0) printf("\n");
    }
}

void printAsTable(int * matrix, int nx, int ny){
    int i;
    for (i = 1; i <= nx * ny; i++){
        int b_i = (i-1)/nx, b_j = (i-1)%nx; 
        printf("%d\t%d\t%d\n",b_i,b_j,(matrix[i-1]==2) ? 0 :matrix[i-1] );
    }
}

int main(int argc, char const **argv){

    #ifndef ITERATIONS
    printf("Number of iterations is needed please compile with -D ITERATIONS=<int>\n");
    return 0;
    #endif

    int Nx = 500;
    int Ny = 1397;
    int * Ma = (int*) malloc( Nx * Ny * sizeof(int));
    for (int i = 0; i < Nx * Ny; ++i) Ma[i] = 0;
    int * Mb = (int*) malloc( Nx * Ny * sizeof(int));
    for (int i = 0; i < Nx * Ny; ++i) Mb[i] = 0;

    mapload(Ma,Nx,Ny);/*Load de Map contour from the .dat file */
    estload(Ma,Nx,Ny);/*Load de Coordinates of the estations from the .dat file*/
    
    clock_t start = 0.0, end = 0.0;
    double sum = 0.0;
    start = clock();
    for(int iters=0;iters<ITERATIONS;iters++){
        for(int i=0;i< Nx*Ny; i++) 
            Mb[i] = ( Ma[i]==1 || Ma[i]==2 ) ? Ma[i] : propagation(Ma,Nx,Ny,i);

        int * temp = Ma;
        Ma = Mb;
        Mb = temp;
    }
    end = clock();
    sum = (end -start) / (double) CLOCKS_PER_SEC;

    #ifdef TIME 
    savetime(sum);
    #endif

    #if SAVE == 1
    saveAsTable(Ma,Nx,Ny);
    #elif SAVE == 2
    saveAsMatrix(Ma,Nx,Ny);
    #else
    saveAsTable(Ma,Nx,Ny);
    saveAsMatrix(Ma,Nx,Ny);
    #endif
        
    free(Ma);
    free(Mb);
    
    return 0;
}
