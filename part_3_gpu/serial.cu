#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

#define InitialCapacity 10
#define RegionSize 0.01
extern double size;
int indexing(int region_ind_m, int region_ind_n, int lda, int local_ind){
    return (region_ind_m+lda*region_ind_n)*InitialCapacity+local_ind;
}

void init_region_num(int* region_num, int total_region){
    for(int tid = 0; tid < total_region; tid++)
    region_num[tid] = 0;
}

void init_region_list(int* region_list, int* region_num, particle_t* particles, int n, int lda){
    for(int tid = 0; tid < n; tid++){
        particle_t& temp = particles[tid];
        int m_ind = temp.x/RegionSize, n_ind = temp.y/RegionSize;
        int local_ind = region_num[m_ind+n_ind*lda];
        region_num[m_ind+n_ind*lda]++;
        region_list[indexing(m_ind, n_ind, lda, local_ind)] = tid;
    }
}

void compute_forces(int* region_list, int*region_num, particle_t* particles, int n, int lda){
    for(int tid = 0; tid < n; tid++){
        particles[tid].ax = particles[tid].ay = 0;
        int m_ind = particles[tid].x/RegionSize, n_ind = particles[tid].y/RegionSize;
        for(int m = max(m_ind-1, 0); m <= min(m_ind+1, lda-1); m++){
            for(int n = max(n_ind-1, 0); n <= min(n_ind+1, lda-1); n++){
                for(int i = 0; i < region_num[m+n*lda]; i++){
                    apply_force(particles[tid], particles[region_list[indexing(m, n, lda, i)]]);
                }
            }
        }
    }
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );

    int lda_region = ceil(size/RegionSize);
    int total_region = lda_region*lda_region;

    init_particles( n, particles );

    int* region_list = (int*) malloc(total_region*InitialCapacity*sizeof(int)); 
    int* region_num = (int*) malloc(total_region*sizeof(int));
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        init_region_num(region_num, total_region);
        init_region_list(region_list, region_num, particles, n, lda_region);
        //
        //  compute forces
        //
        compute_forces(region_list, region_num, particles, n, lda_region);
        
        //
        //  move particles
        //
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );
        
        //
        //  save if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    free(region_list); free(region_num);
    if( fsave )
        fclose( fsave );
    
    return 0;
}
