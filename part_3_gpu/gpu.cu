#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>
#include "common.h"

#define NUM_THREADS 128
#define InitialCapacity 10
#define binSize 0.01

extern double size;
//
//  benchmarking program
//

__device__ void apply_force_gpu(particle_t &particle, particle_t &neighbor)
{
  double dx = neighbor.x - particle.x;
  double dy = neighbor.y - particle.y;
  double r2 = dx * dx + dy * dy;
  if( r2 > cutoff*cutoff )
      return;
  //r2 = fmax( r2, min_r*min_r );
  r2 = (r2 > min_r*min_r) ? r2 : min_r*min_r;
  double r = sqrt( r2 );

  //
  //  very simple short-range repulsive force
  //
  double coef = ( 1 - cutoff / r ) / r2 / mass;
  particle.ax += coef * dx;
  particle.ay += coef * dy;

}

__device__ int indexing(int bin_ind_m, int bin_ind_n, int lda, int local_ind){
    return (bin_ind_m+lda*bin_ind_n)*InitialCapacity+local_ind;
}

__global__ void move_gpu (particle_t * particles, int n, double size, int* bin_list, int* bin_num, int lda)
{

  // Get thread (particle) ID
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid >= n) return;

    particle_t * p = &particles[tid];
    int m_old = (p->x)/binSize, n_old = (p->y)/binSize;
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p->vx += p->ax * dt;
    p->vy += p->ay * dt;
    p->x  += p->vx * dt;
    p->y  += p->vy * dt;

    //
    //  bounce from walls
    //
    while( p->x < 0 || p->x > size )
    {
        p->x  = p->x < 0 ? -(p->x) : 2*size-p->x;
        p->vx = -(p->vx);
    }
    while( p->y < 0 || p->y > size )
    {
        p->y  = p->y < 0 ? -(p->y) : 2*size-p->y;
        p->vy = -(p->vy);
    }

    int m_new = (p->x)/binSize, n_new = (p->y)/binSize;
    if(m_new != m_old || n_new != n_old){
        for(int i = 0; i < bin_num[m_old + n_old * lda]; i++){
            if(bin_list[indexing(m_old, n_old, lda, i)] == tid){
                bin_list[indexing(m_old, n_old, lda, i)] = -1;
                break;
            }
        }

        int index = atomicAdd(&bin_num[m_new + n_new * lda], 1);
        bin_list[indexing(m_new, n_new, lda, index)] = tid;
    }
}

__global__ void clean_bins(int* bin_list, int* bin_num, int lda){
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid >= lda * lda) return;

    int ind = 0;
    for(int i = 0; i < bin_num[tid]; i++){
        if(bin_list[tid * InitialCapacity + i] != -1){
            bin_list[tid * InitialCapacity + ind] = bin_list[tid * InitialCapacity + i];
            ind++;
        }
    }
    bin_num[tid] = ind;
}

__global__ void init_bin_num(int* bin_num, int total_bin){
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid >= total_bin) return;
    bin_num[tid] = 0;
}

__global__ void init_bin_list(int* bin_list, int* bin_num, particle_t* particles,
        int n, int lda){
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid >= n) return;
    particle_t& temp = particles[tid];
    int m_ind = temp.x/binSize, n_ind = temp.y/binSize;
    int local_ind = atomicAdd(&bin_num[m_ind+n_ind*lda], 1);
    bin_list[indexing(m_ind, n_ind, lda, local_ind)] = tid;
}

__global__ void compute_forces_gpu(int* bin_list, int*bin_num, particle_t* particles,
        int n, int lda){
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid >= n) return;
    particles[tid].ax = particles[tid].ay = 0;
    int m_ind = particles[tid].x/binSize, n_ind = particles[tid].y/binSize;
    for(int m = max(m_ind-1, 0); m <= min(m_ind+1, lda-1); m++){
        for(int n = max(n_ind-1, 0); n <= min(n_ind+1, lda-1); n++){
            for(int i = 0; i < bin_num[m+n*lda]; i++){
                apply_force_gpu(particles[tid], particles[bin_list[indexing(m, n, lda, i)]]);
            }
        }
    }
}

int main( int argc, char **argv )
{
    // This takes a few seconds to initialize the runtime
    cudaThreadSynchronize();

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

    // GPU particle data structure
    particle_t * d_particles;
    cudaMalloc((void **) &d_particles, n * sizeof(particle_t));

    set_size( n );

    int lda_bin = ceil(size/binSize);
    int total_bin = lda_bin*lda_bin;

    init_particles( n, particles );

    int* bin_list; cudaMalloc((void **) &bin_list, total_bin*InitialCapacity*sizeof(int));
    int* bin_num; cudaMalloc((void **) &bin_num, total_bin*sizeof(int));
    int particle_blks = (n + NUM_THREADS - 1) / NUM_THREADS;
    int bin_blks = (total_bin + NUM_THREADS - 1) / NUM_THREADS;

    cudaThreadSynchronize();
    double copy_time = read_timer( );

    // Copy the particles to the GPU
    cudaMemcpy(d_particles, particles, n * sizeof(particle_t), cudaMemcpyHostToDevice);

    cudaThreadSynchronize();
    copy_time = read_timer( ) - copy_time;

    //
    //  simulate a number of time steps
    //
    cudaThreadSynchronize();
    init_bin_num <<< bin_blks, NUM_THREADS >>> (bin_num, total_bin);
    init_bin_list <<< particle_blks, NUM_THREADS >>> (bin_list, bin_num, d_particles, n, lda_bin);
    double simulation_time = read_timer( );

    for( int step = 0; step < NSTEPS; step++ )
    {
        //
        //  compute forces
        //
	compute_forces_gpu <<< particle_blks, NUM_THREADS >>> (bin_list, bin_num, d_particles, n, lda_bin);

        //
        //  move particles
        //
	move_gpu <<< particle_blks, NUM_THREADS >>> (d_particles, n, size, bin_list, bin_num, lda_bin);
    clean_bins <<< bin_blks, NUM_THREADS >>> (bin_list, bin_num, lda_bin);


    /*init_bin_num <<< bin_blks, NUM_THREADS >>> (bin_num, bin_add_num, total_bin);*/
    /*init_bin_list <<< particle_blks, NUM_THREADS >>> (bin_list, bin_num, d_particles, n, lda_bin);*/
        //
        //  save if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 ) {
	    // Copy the particles back to the CPU
            cudaMemcpy(particles, d_particles, n * sizeof(particle_t), cudaMemcpyDeviceToHost);
            save( fsave, n, particles);
	}
    }
    cudaThreadSynchronize();
    simulation_time = read_timer( ) - simulation_time;

    printf( "CPU-GPU copy time = %g seconds\n", copy_time);
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );

    free( particles );
    cudaFree(d_particles);
    cudaFree(bin_list); cudaFree(bin_num);
    if( fsave )
        fclose( fsave );

    return 0;
}
