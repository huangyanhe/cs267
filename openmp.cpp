#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "decomposition.h"
#include "omp.h"

int InitialCapacity = 5;
double RegionSize = 0.01;
int main( int argc, char **argv )
{
    int navg,nabsavg=0,numthreads;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );

#pragma omp parallel
    {
    numthreads = omp_get_num_threads();
    }

    decomp decomposition(n, particles);
    decomposition.sub_decomp(numthreads);

    omp_lock_t* region_lock = (omp_lock_t*)malloc(decomposition.Num_region*sizeof(omp_lock_t));
    for(int i = 0; i < decomposition.Num_region; i++){
        omp_init_lock(&region_lock[i]);
    }
    //
    //  simulate a number of time steps
    //

    double simulation_time = read_timer( );
#pragma omp parallel private(dmin)
    {
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        davg = 0.0;
        dmin = 1.0;
        //
        //  compute forces
        //
#pragma omp for schedule(static) reduction(+:navg) reduction(+:davg)
        for(int i = 0; i < n; i++){
            particles[i].ax = particles[i].ay = 0;
            int m = particles[i].x/RegionSize, n = particles[i].y/RegionSize;
            for(int mInd = max(m-1,0); mInd <= min(m+1,decomposition.M-1); mInd++){
                for(int nInd = max(n-1,0); nInd <= min(n+1,decomposition.M-1); nInd++){
                    region& temp = decomposition(mInd, nInd);
                    for(int k = 0; k < temp.Num; k++){
                        apply_force(particles[i], particles[temp.ind[k]], &dmin, &davg, &navg);
                    }
                }
            }
        }
        //
        // move particles
        //
#pragma omp for schedule(dynamic) collapse(2)
    for(int i = 0; i < decomposition.num_sub_M; i++){
        for(int j = 0; j < decomposition.num_sub_N; j++){

            for(int m = decomposition.grid_M[i]; m < decomposition.grid_M[i+1]; m++){
                for(int n = decomposition.grid_N[j]; n < decomposition.grid_N[j+1]; n++){
                    for(int s = 0, k = 0; s < decomposition.region_length[m+n*decomposition.M]; s++){
                        int index = decomposition(m,n).ind[k];
                        move(particles[index]);
                        int m_new = particles[index].x/RegionSize, n_new = particles[index].y/RegionSize;
                        if(m_new != m || n_new != n){
                            if(m_new >= decomposition.grid_M[i] && m_new < decomposition.grid_M[i+1] &&
                                    n_new >= decomposition.grid_N[j] && n_new < decomposition.grid_N[j+1]){
                                decomposition.delete_particle(index, m, n);
                                decomposition.add_particle(index, m_new, n_new);
                            }
                            else{
                                omp_set_lock(&region_lock[m+n*decomposition.M]);
                                decomposition.delete_particle(index, m, n);
                                omp_unset_lock(&region_lock[m+n*decomposition.M]);

                                omp_set_lock(&region_lock[m_new+n_new*decomposition.M]);
                                decomposition.add_particle(index, m_new, n_new);
                                omp_unset_lock(&region_lock[m_new+n_new*decomposition.M]);
                            }
                        }
                        else{
                            k++;
                        }
                    }
                }
            }
        }
    }
#pragma omp for schedule(static)
        for(int i = 0; i < decomposition.Num_region; i++){
            decomposition.region_length[i] = decomposition.region_list[i].Num;
        }


//#pragma omp for schedule(static)
        //for( int i = 0; i < n; i++ ){
            //int m_old = particles[i].x/RegionSize, n_old = particles[i].y/RegionSize;
            //move(particles[i]);
            //int m_new = particles[i].x/RegionSize, n_new = particles[i].y/RegionSize;
            //if( m_new != m_old || n_new != n_old ){
                //omp_set_lock(&region_lock[m_old+n_old*decomposition.M]);
                //decomposition.delete_particle(i, m_old, n_old);
                //omp_unset_lock(&region_lock[m_old+n_old*decomposition.M]);

                //omp_set_lock(&region_lock[m_new+n_new*decomposition.M]);
                //decomposition.add_particle(i, m_new, n_new);
                //omp_unset_lock(&region_lock[m_new+n_new*decomposition.M]);
            //}
        //}

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
#pragma omp master
        if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
#pragma omp critical
          if (dmin < absmin) absmin = dmin;

          //
          //  save if necessary
          //
#pragma omp master
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    }
    simulation_time = read_timer( ) - simulation_time;

    for(int i = 0; i < decomposition.Num_region; i++){
        omp_destroy_lock(&region_lock[i]);
    }
    free(region_lock);

    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    //
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");

    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );
    free( particles );
    if( fsave )
        fclose( fsave );

    return 0;
}
