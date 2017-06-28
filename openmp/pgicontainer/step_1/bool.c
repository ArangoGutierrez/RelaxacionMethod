    for(int iters=0;iters<ITERATIONS;iters++){
        #pragma omp parallel num_threads(NUM_THREADS)
        {
        #pragma omp for
        for(int i=0;i< Nx*Ny; i++) 
            Mb[i] = ( Ma[i]==1 || Ma[i]==2 ) ? Ma[i] : propagation(Ma,Nx,Ny,i);
        }
        int * temp = Ma;
        Ma = Mb;
        Mb = temp;
    }