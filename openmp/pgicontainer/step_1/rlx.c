void nextstate(int * B, double * Tin, double * Tout, int nx, int ny){
    int i = 0;
    #pragma omp parallel num_threads(NUM_THREADS)
    {
    #pragma omp for
    for (i = 0; i < nx * ny; ++i) Tout[i] = transition(B,Tin,nx,ny,i);
    }
}