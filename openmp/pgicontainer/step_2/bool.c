    int iters = 0;
    #pragma omp parallel private(iters) shared(Ma,Mb,Nx,Ny) num_threads(NUM_THREADS)
    {
        for(iters=0;iters<ITERATIONS;iters++){
            #pragma omp for
            for(int i=0;i< Nx*Ny; i++) 
                Mb[i] = ( Ma[i]==1 || Ma[i]==2 ) ? Ma[i] : propagation(Ma,Nx,Ny,i);
            
            #pragma omp single
            {
                int * temp = Ma;
                Ma = Mb;
                Mb = temp;
            }
        }
    }