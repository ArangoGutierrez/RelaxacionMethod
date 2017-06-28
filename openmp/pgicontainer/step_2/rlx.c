    int gen = 0;
    #pragma omp parallel num_threads(NUM_THREADS) private(gen) shared(B,Ta,Tb,Nx,Ny)
    {
        for (int gen = 1 ; gen <= GENERATIONS; gen++) {

            #pragma omp for
            for (int i = 0; i < Nx * Ny; ++i) Tb[i] = transition(B,Ta,Nx,Ny,i);

            #pragma omp single
            {    
                double * temp = Ta;
                Ta = Tb;
                Tb = temp;
            }
            #ifdef SAVEALL 
            sprintf(filename,"../OutputData/RlxMthd_v1.0_%d.dat",gen);
            savetemperatures(filename,Ta,Nx,Ny);
            #endif
        }
    }