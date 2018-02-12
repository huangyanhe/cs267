#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <iostream>
#include <algorithm>
#include <vector>


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
void findNeighbors(std::vector<int > &a_neighbors, std::vector< std::vector< std::vector< int> > > &bins, int i, int j, int numCells)
{
  a_neighbors.insert(a_neighbors.end(), bins[i][j].begin(), bins[i][j].end());	  
  if (i == 0)
    {
      if (j == 0)
	{
	  a_neighbors.insert(a_neighbors.end(), bins[i+1][j].begin(), bins[i+1][j].end());
	  a_neighbors.insert(a_neighbors.end(), bins[i+1][j+1].begin(),bins[i+1][j+1].end());
	  a_neighbors.insert(a_neighbors.end(), bins[i][j+1].begin(), bins[i][j+1].end());
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
    }
  return;
};
int main( int argc, char **argv )
{    
    int navg,nabsavg=0;
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
    
    int n = read_int( argc, argv, "-n", 10 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );

    
    // for(int p=0; p< n; p++)
    //   {
    // 	maxX = std::max(particles[p].x, maxX);
    // 	maxY = std::max(particles[p].y, maxY);
    // 	minX = std::min( particles[p].x, minX);
   // 	minY = std::min( particles[p].y, minY);
    //   }

    // std::cout<<"Max = ("<< maxX<<","<<maxY<< ")"<<std::endl;
    // std::cout<<"Min = ("<< minX<<","<<minY<< ")"<<std::endl;

    //double xedgeL = minX - 0.5;
    //double xedgeR = maxX + 0.5;
    //double yedgeB = minY - 0.5;
    //double yedgeT = maxY + 0.5;
    
    //int numYcells = std::ceil(std::abs(yedgeT - yedgeB)/ );

    double density = 0.0005;
    double cutoff = 0.01;
    double h = cutoff;
    double maxXY = sqrt(n*density)+ h;
    int numCells = ceil(maxXY/h);
    //int numCells =2;
    //std::vector< int> temp1;
    //std::vector< std::vector< int> > temp2(numCells,  temp1);
    //std::vector< std::vector< std::vector< int> > > bins(numCells, std::vector< std::vector< int> >(numCells));

    
    //bin(bins, particles, n, h);
    
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
	
    for( int step = 0; step < NSTEPS; step++ )
    {
	navg = 0;
        davg = 0.0;
	dmin = 1.0;

	std::vector< std::vector< std::vector< int> > > bins(numCells, std::vector< std::vector< int> >(numCells));
    bin(bins, particles, n, h);
    
	//
        // compute forces
        //
	
        // for( int i = 0; i < n; i++ )
        // {
        //     particles[i].ax = particles[i].ay = 0;
        //     for (int j = 0; j < n; j++ )
	// 			apply_force( particles[i], particles[j],&dmin,&davg,&navg);
        // }
    for (int p = 0; p<n; p++)
      {
    	std::cout<<p<<","<<particles[p].x<<std::endl;
      }
    
	for( int i = 0; i< numCells; i++)
	  {
	    for( int j = 0; j< numCells; j++)
	      {
		std::vector<int> Neighbors;
		findNeighbors(Neighbors, bins, i, j, numCells);
		std::cout<<"Neigbors ij ("<<i<<","<<j<<") = ";
		for(auto pr = Neighbors.begin(); pr!= Neighbors.end(); ++pr )
		  {
		    std::cout<<*pr<<",";
		  }
		std::cout<<std::endl;
		for(int k = 0; k<bins[i][j].size(); k++ )
		  {
		    for(int q=0; q<Neighbors.size(); q++)
		      {			
			apply_force( particles[bins[i][j][k]], particles[Neighbors[q]],&dmin,&davg,&navg);
		      }
		  }
	      }
	  }
	
        //
        //  move particles
        //
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );
	

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;
		
          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;
    
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
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
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
