#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include "common.h"
#include "omp.h"

void findNeighbors(std::vector<int > &a_neighbors, std::vector< std::vector< std::vector< int> > > &bins, int i, int j, int numCells)
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
  for(int i = 0 ; i < bins.size(); i++} {
    for(int j = 0; j < bins[i].size(); i++)
      bins[i][j].clear();
  }
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
    int navg,nabsavg=0,numthreads; 
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
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

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    double density = 0.0005;
    double cutoff = 0.01;
    double h = cutoff;
    double maxXY = sqrt(n*density);
    int numCells = ceil(maxXY/h);
    int BlockSizeM = floor(numCells/4.0);
    int BlockSizeN = floor(numCells/8.0);

    double bin_time = 0;
    double force_time = 0;
    double move_time = 0;
    double before_bin, before_move, before_force;
    std::vector< std::vector< std::vector< int> > > bins(numCells, std::vector< std::vector< int> >(numCells));    
    #pragma omp parallel private(dmin) 
    {
    numthreads = omp_get_num_threads();
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        davg = 0.0;
	dmin = 1.0;
#pragma omp master
	{
	  before_bin = omp_get_wtime();
	}

#pragma omp single
	{
	  bin(bins, particles, n, h);
	}
        //
        //  compute all forces
        //
#pragma omp master
	{
	  before_force = omp_get_wtime();
	}
#pragma omp for reduction (+:navg) reduction(+:davg) collapse(2)
	for( int i = 0; i< numCells; i++)
	  {
	    for( int j = 0; j< numCells; j++)
	      {
		std::vector<int> Neighbors;
		findNeighbors(Neighbors, bins, i, j, numCells);
		for(int k = 0; k<bins[i][j].size(); k++ )
		  {
		    for(int q=0; q<Neighbors.size(); q++)
		      {			
			apply_force( particles[bins[i][j][k]], particles[Neighbors[q]],&dmin,&davg,&navg);
		      }
		  }
	      }
	  }
#pragma omp master
	{
	  double after_force = omp_get_wtime();
	  force_time += (after_force - before_force)/NSTEPS;
	}      
		
        //
        //  move particles
        //
#pragma omp master
	{
	  before_move = omp_get_wtime();
	}
#pragma omp for schedule(static)
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );
#pragma omp master
	{
	  double after_move = omp_get_wtime();
	  move_time += (after_move - before_move)/NSTEPS;
	}

        if( find_option( argc, argv, "-no" ) == -1 ) 
        {
          //
          //  compute statistical data
          //
#pragma omp master
          if (navg) { 
            absavg += davg/navg;
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
    printf("bin_time =  %lf, force_time = %lf, move_time = %lf", bin_time, force_time, move_time);

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
