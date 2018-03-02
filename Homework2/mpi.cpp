#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

// void particleComputeMove(particle_t particles, int n)
// {
//   for( int i = 0; i < n; i++ )
//         {
//             particles[i].ax = particles[i].ay = 0;
//             for (int j = 0; j < n; j++ )
// 				apply_force( particles[i], particles[j],&dmin,&davg,&navg);
//         }
 
//         //
//         //  move particles
//         //
//         for( int i = 0; i < n; i++ ) 
//             move( particles[i] );
// };

//
//  benchmarking program
//
void findNeighbors(std::vector<int >& a_neighbors, std::vector< std::vector< std::vector< int> > >& bins, int i, int j, int numCells )
{
  a_neighbors.insert(a_neighbors.end(), bins[i][j].begin(), bins[i][j].end());	  
 if (i == 0)
    {
      if (j == 0)
	{
	  if (numCells == 1)
	    {
	      return;
	    }
	  else
	    {
	      a_neighbors.insert(a_neighbors.end(), bins[i+1][j].begin(), bins[i+1][j].end());
	      a_neighbors.insert(a_neighbors.end(), bins[i+1][j+1].begin(),bins[i+1][j+1].end());
	      a_neighbors.insert(a_neighbors.end(), bins[i][j+1].begin(), bins[i][j+1].end());
	    }
	}
      else if (j == numCells-1)
	{
	  a_neighbors.insert(a_neighbors.end(), bins[i+1][j].begin(), bins[i+1][j].end());
	  a_neighbors.insert(a_neighbors.end(), bins[i+1][j-1].begin(),bins[i+1][j-1].end());
	  a_neighbors.insert(a_neighbors.end(), bins[i][j-1].begin(), bins[i][j-1].end());

	}
      else
	{
	  a_neighbors.insert(a_neighbors.end(), bins[i+1][j].begin(), bins[i+1][j].end());
	  a_neighbors.insert(a_neighbors.end(), bins[i+1][j+1].begin(),bins[i+1][j+1].end());
	  a_neighbors.insert(a_neighbors.end(), bins[i+1][j-1].begin(),bins[i+1][j-1].end());
	  a_neighbors.insert(a_neighbors.end(), bins[i][j+1].begin(), bins[i][j+1].end());
	  a_neighbors.insert(a_neighbors.end(), bins[i][j-1].begin(), bins[i][j-1].end());
	}
    }
  else if ( i == numCells-1)
    {
      if (j == 0)
	{
	  a_neighbors.insert(a_neighbors.end(), bins[i-1][j].begin(), bins[i-1][j].end());
	  a_neighbors.insert(a_neighbors.end(), bins[i-1][j+1].begin(),bins[i-1][j+1].end());
	  a_neighbors.insert(a_neighbors.end(), bins[i][j+1].begin(), bins[i][j+1].end());
	}
      else if (j == numCells-1)
	{
	  a_neighbors.insert(a_neighbors.end(), bins[i-1][j].begin(), bins[i-1][j].end());
	  a_neighbors.insert(a_neighbors.end(), bins[i-1][j-1].begin(),bins[i-1][j-1].end());
	  a_neighbors.insert(a_neighbors.end(), bins[i][j-1].begin(), bins[i][j-1].end());

	}
      else
	{
	  a_neighbors.insert(a_neighbors.end(), bins[i-1][j].begin(), bins[i-1][j].end());
	  a_neighbors.insert(a_neighbors.end(), bins[i-1][j+1].begin(),bins[i-1][j+1].end());
	  a_neighbors.insert(a_neighbors.end(), bins[i-1][j-1].begin(),bins[i-1][j-1].end());
	  a_neighbors.insert(a_neighbors.end(), bins[i][j+1].begin(), bins[i][j+1].end());
	  a_neighbors.insert(a_neighbors.end(), bins[i][j-1].begin(), bins[i][j-1].end());
	}
    }
  else if (j == 0)
    {
      a_neighbors.insert(a_neighbors.end(), bins[i-1][j].begin(), bins[i-1][j].end());
      a_neighbors.insert(a_neighbors.end(), bins[i-1][j+1].begin(),bins[i-1][j+1].end());
      a_neighbors.insert(a_neighbors.end(), bins[i][j+1].begin(), bins[i][j+1].end());
      a_neighbors.insert(a_neighbors.end(), bins[i+1][j].begin(), bins[i+1][j].end());
      a_neighbors.insert(a_neighbors.end(), bins[i+1][j+1].begin(),bins[i+1][j+1].end());

    }
  else if (j == numCells-1)
    {
      a_neighbors.insert(a_neighbors.end(), bins[i-1][j].begin(), bins[i-1][j].end());
      a_neighbors.insert(a_neighbors.end(), bins[i-1][j-1].begin(),bins[i-1][j-1].end());
      a_neighbors.insert(a_neighbors.end(), bins[i][j-1].begin(), bins[i][j-1].end());
      a_neighbors.insert(a_neighbors.end(), bins[i+1][j].begin(), bins[i+1][j].end());
      a_neighbors.insert(a_neighbors.end(), bins[i+1][j-1].begin(),bins[i+1][j-1].end());
    }
  else
    {
      a_neighbors.insert(a_neighbors.end(), bins[i-1][j+1].begin(),bins[i-1][j+1].end());
      a_neighbors.insert(a_neighbors.end(), bins[i-1][j].begin(), bins[i-1][j].end());
      a_neighbors.insert(a_neighbors.end(), bins[i-1][j-1].begin(),bins[i-1][j-1].end());
      a_neighbors.insert(a_neighbors.end(), bins[i][j+1].begin(), bins[i][j+1].end());
      a_neighbors.insert(a_neighbors.end(), bins[i][j-1].begin(), bins[i][j-1].end());
      a_neighbors.insert(a_neighbors.end(), bins[i+1][j].begin(), bins[i+1][j].end());
      a_neighbors.insert(a_neighbors.end(), bins[i+1][j-1].begin(),bins[i+1][j-1].end());
      a_neighbors.insert(a_neighbors.end(), bins[i+1][j+1].begin(),bins[i+1][j+1].end()); 
    }
  
  
};
void bin(std::vector< std::vector< std::vector< int> > > &bins, particle_t *particles, int n, double h)
{
  int iposx, iposy;
  for(int p=0; p< n; p++)
    {
      iposx = floor(particles[p].x/h);
      iposy = floor(particles[p].y/h);  
      bins[iposx][iposy].push_back(p);
      // set acceleration to zero
      particles[p].ax = 0;
      particles[p].ay = 0;
    }
  return;
};
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
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
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
    int particle_per_proc = (n + n_proc - 1) / n_proc;
    int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
    for( int i = 0; i < n_proc+1; i++ )
        partition_offsets[i] = min( i * particle_per_proc, n );
    
    int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
    for( int i = 0; i < n_proc; i++ )
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
    MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );

    double density = 0.0005;	
    double cutoff = 0.01;
    //double cutoff = sqrt(n*density);
    double h = cutoff;
    double maxXY = sqrt(n*density);
    int numCells = ceil(maxXY/h);
    
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
        MPI_Allgatherv( local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );

	std::vector< std::vector< std::vector< int> > > bins(numCells, std::vector< std::vector< int> >(numCells));

	
    bin(bins, particles, n, h);
    
	
        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if( find_option( argc, argv, "-no" ) == -1 )
          if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
        
        //
        //  compute all forces
        //
	for( int i = 0; i < nlocal; i++ )
        {
	  local[i].ax = local[i].ay = 0;
	  int iposx = floor(local[i].x/h);
	  int iposy = floor(local[i].y/h);
	  std::vector<int> Neighbors;
	  findNeighbors(Neighbors, bins, iposx, iposy, numCells);
	  for(int q=0; q<Neighbors.size(); q++)
	    {			
	      apply_force( local[i], particles[Neighbors[q]],&dmin,&davg,&navg);
	    }
	}
	
	
        // for( int i = 0; i < nlocal; i++ )
        // {
        //     local[i].ax = local[i].ay = 0;
        //     for (int j = 0; j < n; j++ )
        //         apply_force( local[i], particles[j], &dmin, &davg, &navg );
        // }
     
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

        //
        //  move particles
        //
        for( int i = 0; i < nlocal; i++ )
            move( local[i] );
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
        fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
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
