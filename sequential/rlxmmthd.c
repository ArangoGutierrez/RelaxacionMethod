/*
    Relaxacion method for temperature interpolation

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
#include <time.h>   /* To performance evaluation */

struct Station
{
    int x;      /* x coordinate */
    int y;      /* y coordinate */
    double t;   /* temperature */
};


void boolmapload(int* B,int Nx,int Ny){
    FILE *file = fopen("../OutputData/boolmap_table.dat","r");
    int x = 0, y = 0, b = 0;
    for(; fscanf(file, "%d\t%d\t%d", &y, &x, &b) && !feof(file);) B[ y * Nx + x ] = b;
    fclose(file);   
}

void loadstations(struct Station * s, int ne){
    // FILE * file = fopen("../Maps/estaciones.dat", "r");
    // int x = 0, y = 0, i = 0;
    // for(; fscanf(file, "%d %d %lf",&x,&y,&tmp) && !feof(file);){
    //     s[i].x = x;
    //     s[i].y = y;
    //     s[i].t = tmp;
    //     i++;
    // }
    // double tmp = 0.0;
    // fclose(file);
    s[0].x = 297;
    s[0].y = (1400 - 481);
    s[0].t = 38.2490;

    s[1].x = 54;
    s[1].y = 175;
    s[1].t = 20.7490;
}

void putstations(int* B, int Nx,int Ny, struct Station * s, int ne){
    int i = 0;
    for(int i=0; i<ne; ++i) B[ s[i].y * Nx + s[i].x ] = 2;
}

void puttemperatures(int * B, double * matrix, int nx, int ny, struct Station * s, int ne){
    double sum = 0;
    int i;
    for (i = 0; i < ne; i++) sum = sum + s[i].t;
    double Tprom =  ( sum / ne );
    for (i=0; i < nx * ny; ++i) matrix[i] = ( B[i] == 0 )? 0 : Tprom;
    for (i=0; i < ne ; ++i) matrix[( s[i].y * nx + s[i].x ) ] = s[i].t;
}

void prntB(int * matrix, int nx, int ny){
    int i;
    for (i = 1; i <= nx * ny; i++)
    {
        printf("%d\t",matrix[i-1]);
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

void savetemperatures(const char * filename, double * matrix, int nx, int ny){
    FILE * pf = fopen(filename,"w");
    int i = 0, x = 0, y = 0;
    for (i = 0; i < nx * ny; ++i){
        x = ( i % nx ) + 1;
        y = ( i / nx ) + 1;
        fprintf(pf,"%d\t%d\t%f\n",x,y,matrix[i]);
    }
    fclose(pf);
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
    // printf("Test Result delta: %.10f\n",result);
    // result = (1 / (nx*ny)) * sqrt(result);
    result =  sqrt(result);
    //printf("Test Result: %.10f\n",result);
    return ( result <= 0.0000001 ) ? 1 : 0;
}

void nextstate(int * B, double * Tin, double * Tout, int nx, int ny){
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

void savetime(double time){
    FILE * pf = fopen("time_relaxation.txt","a+");
    fprintf(pf, "%f\n",time);
    fclose(pf); 
}

int main(int argc, char const **argv)
{

    #ifndef GENERATIONS
    printf("Number of generations is needed please compile with -D GENERATIONS=<int>\n");
    return 0;
    #endif
    
    int Nx = 500;
    int Ny = 1397;
    int Ne = 2;
    char filename[50];
    
    struct Station s[Ne];

    int * B = NULL;
    double * Ta = NULL;
    double * Tb = NULL;

    B = (int *) malloc( Nx * Ny * sizeof(int));
    Ta = (double *) malloc(Nx * Ny * sizeof(double));
    Tb = (double *) malloc(Nx * Ny * sizeof(double));
    
    boolmapload(B,Nx,Ny);
    loadstations(s,Ne);
    putstations(B,Nx,Ny,s,Ne);
    puttemperatures(B,Ta,Nx,Ny,s,Ne);

    #ifdef SAVEINIT
    sprintf(filename,"../OutputData/RlxMthd_v1.0_%d.dat",0);
    savetemperatures(filename,Ta,Nx,Ny);
    #endif    

    clock_t start = 0.0, end = 0.0;
    double sum = 0.0;

    printf("Checking %d\n",4);
    #if GENERATIONS == 0
    /* The system evolves until it reaches a steady state (convergence) */
    start = clock();
    while( test(B,Ta,Tb,Nx,Ny) != 1 ){

        nextstate(B,Ta,Tb,Nx,Ny);
        
        double * temp = Ta;
        Ta = Tb;
        Tb = temp;

        #ifdef SAVEALL 
        sprintf(filename,"../OutputData/RlxMthd_v1.0_%d.dat",i);
        savetemperatures(filename,Ta,Nx,Ny);
        #endif
    }
    end = clock();
    sum = (end -start) / (double) CLOCKS_PER_SEC;
    #else
    /* The system evolves until it reaches a given number of generations */
    start = clock();
    for (int i = 0 ; i < GENERATIONS; i++) {
        nextstate(B,Ta,Tb,Nx,Ny);
        double * temp = Ta;
        Ta = Tb;
        Tb = temp;
        #ifdef SAVEALL 
        sprintf(filename,"../OutputData/RlxMthd_v1.0_%d.dat",i);
        savetemperatures(filename,Ta,Nx,Ny);
        #endif
    }
    end = clock();
    sum = (end -start) / (double) CLOCKS_PER_SEC;
    #endif

    #ifdef TIME 
    savetime(sum);
    #endif
    
    #ifdef SAVELAST
    sprintf(filename,"../OutputData/RlxMthd_v1.0_%d.dat",GENERATIONS);
    savetemperatures(filename,Ta,Nx,Ny);
    #endif
    
    free(B);
    free(Ta);
    free(Tb);

    return 0;
}
