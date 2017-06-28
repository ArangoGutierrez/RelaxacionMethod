    #pragma acc data copy(Nx,Ny),copy(Ma[0:(500*1397)],Mb[0:(500*1397)]), create(Mtemp)
    {
        for(int iters = 0; iters < ITERATIONS; ++iters){
            #pragma acc parallel loop
            for(int i = 0; i < (500*1397); ++i){ 
                Mb[i] = ( Ma[i] == 1 || Ma[i] == 2 ) ? Ma[i] : propagation(Ma,Nx,Ny,i);;
            }

            #pragma acc parallel loop
            for(int i = 0; i < (500*1397); ++i){
                Ma[i] = Mb[i];
            }
        }
    }