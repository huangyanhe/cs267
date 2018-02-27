#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include "decomposition_mpi.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{
    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg;

    //
    //  process command line parameters
    //
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

    //
    //  set up MPI
    //
    int num_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &num_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;


    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );

    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );

    //
    //  set up the data partitioning across processors
    //
    int particle_per_proc = (n + num_proc - 1) / num_proc;
    int *partition_offsets = (int*) malloc( (num_proc+1) * sizeof(int) );
    for( int i = 0; i < num_proc+1; i++ )
        partition_offsets[i] = min( i * particle_per_proc, n );

    int *partition_sizes = (int*) malloc( num_proc * sizeof(int) );
    for( int i = 0; i < num_proc; i++ )
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];

    //
    //  allocate storage for local partition
    //
    int nlocal = partition_sizes[rank];
    particle_t *local = (particle_t*) malloc( nlocal * sizeof(particle_t) );

    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    set_size( n );
    if( rank == 0 )
        init_particles( n, particles );
    //MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
    MPI_Bcast(particles, n, PARTICLE, 0, MPI_COMM_WORLD);
    //printf("rank %d complete Bcast\n", rank); fflush(stdout);
    decomp_proc local_proc(num_proc, rank, n, particles);

    //printf("rank %d complete initialization\n", rank); fflush(stdout);

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;
        //
        //  collect all global data locally (not good idea to do)
        //
        //MPI_Allgatherv( local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );

        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if( find_option( argc, argv, "-no" ) == -1 )
          if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );

        //
        //  compute all forces
        //
        for(int i = 0; i < local_proc.local_ind.size(); i++){
            particle_t& temp = particles[local_proc.local_ind[i]];
            temp.ax = 0; temp.ay = 0;
            int m_ind = temp.x/local_proc.region_size_m, n_ind = temp.y/local_proc.region_size_n;
            for(int m = max(m_ind-1, local_proc.closure_l_m); m <= min(m_ind+1, local_proc.closure_u_m-1); m++){
                for(int n = max(n_ind-1, local_proc.closure_l_n); n <= min(n_ind+1, local_proc.closure_u_n-1); n++){
                    for(int k = 0; k < local_proc(m,n).size(); k++){
                        apply_force(temp, particles[local_proc(m,n)[k]], &dmin, &davg, &navg);
                    }
                }
            }
        }


        if( find_option( argc, argv, "-no" ) == -1 )
        {

          MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);


          if (rank == 0){
            //
            // Computing statistical data
            //
            if (rnavg) {
              absavg +=  rdavg/rnavg;
              nabsavg++;
            }
            if (rdmin < absmin) absmin = rdmin;
          }
        }

        // clean the boundary
        local_proc.clean_boundary();

        //
        //  move particles
        //
        int num_local_particles = local_proc.local_ind.size();
        for(int i = 0, k = 0; i < num_local_particles; i++, k++){
            particle_t& temp = particles[local_proc.local_ind[k]];
            int m_old = temp.x/local_proc.region_size_m, n_old = temp.y/local_proc.region_size_n;
            move(temp);
            int m_new = temp.x/local_proc.region_size_m, n_new = temp.y/local_proc.region_size_n;
            local_proc.pre_sync(local_proc.local_ind[k], m_new, n_new);
            if(m_new != m_old || n_new != n_old){
                local_proc.delete_particle(m_old, n_old, local_proc.local_ind[k]);
                if(m_new >= local_proc.interior_l_m && m_new < local_proc.interior_u_m
                        && n_new >= local_proc.interior_l_n && n_new < local_proc.interior_u_n){
                    local_proc(m_new, n_new).push_back(local_proc.local_ind[k]);
                }
                else{
                    local_proc.local_delete_particle(k);
                    k--;
                }
            }
        }
        local_proc.synchronization(particles, step);
    }
    simulation_time = read_timer( ) - simulation_time;

    if (rank == 0) {
      printf( "n = %d, simulation time = %g seconds", n, simulation_time);

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
        fprintf(fsum,"%d %d %g\n",n,num_proc,simulation_time);
    }

    //
    //  release resources
    //
    if ( fsum )
        fclose( fsum );
    free( partition_offsets );
    free( partition_sizes );
    free( local );
    free( particles );
    if( fsave )
        fclose( fsave );

    MPI_Finalize( );

    return 0;
}
