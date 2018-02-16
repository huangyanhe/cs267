#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <iostream>
#include <algorithm>
#include <vector>




//
//  benchmarking program
//

void findNeighbors(std::vector<particle_t > &a_neighbors, std::vector< std::vector<particle_t> >& particleBins, int i, int j, int numCells)
{
  a_neighbors.insert(a_neighbors.end(), particleBins[i+numCells*j].begin(), particleBins[i + numCells*j].end());	  
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
	      a_neighbors.insert(a_neighbors.end(), particleBins[i+1 + numCells*j].begin(), particleBins[i+1 + numCells*j].end());
	      a_neighbors.insert(a_neighbors.end(), particleBins[i+1 + (j+1)*numCells].begin(),particleBins[i+1 + (j+1)*numCells].end());
	      a_neighbors.insert(a_neighbors.end(), particleBins[i + (j+1)*numCells].begin(), particleBins[i + (j+1)*numCells].end());
	    }
	}
      else if (j == numCells-1)
	{
	  a_neighbors.insert(a_neighbors.end(), particleBins[i+1 + numCells*j].begin(), particleBins[i+1 + numCells*j].end());
	  a_neighbors.insert(a_neighbors.end(), particleBins[i+1 + (j-1)*numCells].begin(),particleBins[i+1 + (j-1)*numCells].end());
	  a_neighbors.insert(a_neighbors.end(), particleBins[i + (j-1)*numCells].begin(), particleBins[i+(j-1)*numCells].end());

	}
      else
	{
	  a_neighbors.insert(a_neighbors.end(), particleBins[i+1 + numCells*j].begin(), particleBins[i+1 + numCells*j].end());
	  a_neighbors.insert(a_neighbors.end(), particleBins[i+1 + (j+1)*numCells].begin(),particleBins[i+1 + (j+1)*numCells].end());
	  a_neighbors.insert(a_neighbors.end(), particleBins[i+1 + (j-1)*numCells].begin(),particleBins[i+1 +(j-1)*numCells].end());
	  a_neighbors.insert(a_neighbors.end(), particleBins[i + (j+1)*numCells].begin(), particleBins[i+ (j+1)*numCells].end());
	  a_neighbors.insert(a_neighbors.end(), particleBins[i + (j-1)*numCells].begin(), particleBins[i + (j-1)*numCells].end());
	}
    }
  else if ( i == numCells-1)
    {
      if (j == 0)
	{
	  a_neighbors.insert(a_neighbors.end(), particleBins[i-1 + numCells*j].begin(), particleBins[i-1 + numCells*j].end());
	  a_neighbors.insert(a_neighbors.end(), particleBins[i-1 + (j+1)*numCells].begin(),particleBins[i-1 + (j+1)*numCells].end());
	  a_neighbors.insert(a_neighbors.end(), particleBins[i + (j+1)*numCells].begin(), particleBins[i + (j+1)*numCells].end());
	}
      else if (j == numCells-1)
	{
	  a_neighbors.insert(a_neighbors.end(), particleBins[i-1 + numCells*j].begin(), particleBins[i-1 + numCells*j].end());
	  a_neighbors.insert(a_neighbors.end(), particleBins[i-1+ (j-1)*numCells].begin(),particleBins[i-1 + (j-1)*numCells].end());
	  a_neighbors.insert(a_neighbors.end(), particleBins[i + (j-1)*numCells].begin(), particleBins[i + (j-1)*numCells].end());

	}
      else
	{
	  a_neighbors.insert(a_neighbors.end(), particleBins[i-1 + numCells*j].begin(), particleBins[i-1 + numCells*j].end());
	  a_neighbors.insert(a_neighbors.end(), particleBins[i-1 + (j+1)*numCells].begin(),particleBins[i-1 + (j+1)*numCells].end());
	  a_neighbors.insert(a_neighbors.end(), particleBins[i-1+(j-1)*numCells].begin(),particleBins[i-1+(j-1)*numCells].end());
	  a_neighbors.insert(a_neighbors.end(), particleBins[i+(j+1)*numCells].begin(), particleBins[i + (j+1)*numCells].end());
	  a_neighbors.insert(a_neighbors.end(), particleBins[i + (j-1)*numCells].begin(), particleBins[i + (j-1)*numCells].end());
	}
    }
  else if (j == 0)
    {
      a_neighbors.insert(a_neighbors.end(), particleBins[i-1 + numCells*j].begin(), particleBins[i-1 + numCells*j].end());
      a_neighbors.insert(a_neighbors.end(), particleBins[i-1 + (j+1)*numCells].begin(),particleBins[i-1 + (j+1)*numCells].end());
      a_neighbors.insert(a_neighbors.end(), particleBins[i + (j+1)*numCells].begin(), particleBins[i + (j+1)*numCells].end());
      a_neighbors.insert(a_neighbors.end(), particleBins[i+1 + numCells*j].begin(), particleBins[i+1 + numCells*j].end());
      a_neighbors.insert(a_neighbors.end(), particleBins[i+1 + (j+1)*numCells].begin(),particleBins[i+1 + (j+1)*numCells].end());

    }
  else if (j == numCells-1)
    {
      a_neighbors.insert(a_neighbors.end(), particleBins[i-1 + numCells*j].begin(), particleBins[i-1 + numCells*j].end());
      a_neighbors.insert(a_neighbors.end(), particleBins[i-1 + (j-1)*numCells].begin(),particleBins[i-1 + (j-1)*numCells].end());
      a_neighbors.insert(a_neighbors.end(), particleBins[i + (j-1)*numCells].begin(), particleBins[i + (j-1)*numCells].end());
      a_neighbors.insert(a_neighbors.end(), particleBins[i+1 + numCells*j].begin(), particleBins[i+1 + numCells*j].end());
      a_neighbors.insert(a_neighbors.end(), particleBins[i+1 + (j-1)*numCells].begin(),particleBins[i+1 + (j-1)*numCells].end());
    }
  else
    {
      a_neighbors.insert(a_neighbors.end(), particleBins[i-1+(j+1)*numCells].begin(),particleBins[i-1 + (j+1)*numCells].end());
      a_neighbors.insert(a_neighbors.end(), particleBins[i-1 + numCells*j].begin(), particleBins[i-1 + numCells*j].end());
      a_neighbors.insert(a_neighbors.end(), particleBins[i-1 + (j-1)*numCells].begin(),particleBins[i-1 + (j-1)*numCells].end());
      a_neighbors.insert(a_neighbors.end(), particleBins[i + (j+1)*numCells].begin(), particleBins[i + (j+1)*numCells].end());
      a_neighbors.insert(a_neighbors.end(), particleBins[i + (j-1)*numCells].begin(), particleBins[i + (j-1)*numCells].end());
      a_neighbors.insert(a_neighbors.end(), particleBins[i+1 + numCells*j].begin(), particleBins[i+1 + numCells*j].end());
      a_neighbors.insert(a_neighbors.end(), particleBins[i+1 + (j-1)*numCells].begin(),particleBins[i+1 + (j-1)*numCells].end());
      a_neighbors.insert(a_neighbors.end(), particleBins[i+1 + (j+1)*numCells].begin(),particleBins[i+1 + (j+1)*numCells].end()); 
    }
  
  
};
void bin( std::vector< std::vector<particle_t> >& particleBins, particle_t *particles, int n, int numCells,double h)
{
  //Need to think about whether this needs to be here now
  // for(int i = 0 ; i < particleBins.size(); i++} {
  //   for(int j = 0; j < bins[i].size(); j++)
  //     bins[i][j].clear();
  // }
  
  int iposx, iposy;
  for(int p=0; p< n; p++)
    {
      iposx = floor(particles[p].x/h);
      iposy = floor(particles[p].y/h);  
      particleBins[iposx + iposy*numCells].push_back(particles[p]);
    }
  return;
};
void move(std::vector< std::vector<particle_t> >& particleBins, std::vector<std::vector<int> >& leftCell, int numCells, double h)
{
  int iposx, iposy;
  for( int i = 0; i< numCells; i++)
    {
      for( int j = 0; j< numCells; j++)
	{
	  for(int k = 0; k<particleBins[i + numCells*j].size(); k++ )
	    {
	      if (particleBins[i + numCells*j][k].exists == true)
		{
		  move( particleBins[i + numCells*j][k]);
		  // set acceleration to zero
		  particleBins[i + numCells*j][k].ax = 0;
		  particleBins[i + numCells*j][k].ay = 0;
		  iposx = floor(particleBins[i + numCells*j][k].x/h);
		  iposy = floor(particleBins[i + numCells*j][k].y/h);  
		  if (iposx != i || iposy != j)
		    {
		      std::vector<int> Coordinates(5, 0);
		      Coordinates[0] = i;
		      Coordinates[1] = j;
		      Coordinates[2] = k;
		      Coordinates[3] = iposx;
		      Coordinates[4] = iposy;
		      leftCell.push_back(Coordinates);
		    }
		}
	    }
	}
    }
}
/* Need to think about what the correct thing to do in this case is. One option is to delete the particle but maybe better is to leave it alone and just flag it as not used. I could move it outside the domain and anytime I came to the particle for force calculations, I could just ignore it. 
 */ 
void particleArithmetic(std::vector< std::vector<particle_t> >& particleBins, std::vector<std::vector<int> >& leftCell, int numCells, int numDeletedParticles)
{
  for(int j = 0; j<leftCell.size(); j++)
    {
      numDeletedParticles++;
      particleBins[leftCell[j][3] + leftCell[j][3]*numCells].push_back(particleBins[leftCell[j][0] + leftCell[j][1]*numCells][leftCell[j][2]]);
      
      particleBins[leftCell[j][0] + leftCell[j][1]*numCells][leftCell[j][2]].exists = false;
    }
}
void computeForces(std::vector< std::vector<particle_t> >& particleBins, int numCells, double& dmin, double& davg, int& navg)
{
  	for( int i = 0; i< numCells; i++)
	  {
	    for( int j = 0; j< numCells; j++)
	      {
		std::vector<particle_t> Neighbors;
		findNeighbors(Neighbors, particleBins, i, j, numCells);
		for(int k = 0; k<particleBins[i + numCells*j].size(); k++ )
		  {
		    for(int q=0; q<Neighbors.size(); q++)
		      {
			if (particleBins[i+numCells*j][k].exists == true && Neighbors[q].exists == true)
			apply_force( particleBins[i+numCells*j][k], Neighbors[q],&dmin,&davg,&navg);
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
    
    int n = read_int( argc, argv, "-n", 100 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );


    double density = 0.0005;
    double cutoff = 0.01;
    double h = cutoff;
    double maxXY = sqrt(n*density);
    int numCells = ceil(maxXY/h);
    int numDeletedParticles = 0;
    
    int avgParticlesPerCell = ceil(n/numCells);
    int overStorageFactor = 2;
    int maxNumParticlesinBin = avgParticlesPerCell*overStorageFactor;
    //std::vector<particle_t> initVector;
    //initVector.reserve(maxNumParticlesinBin);
    std::vector< std::vector<particle_t> > particleBins(numCells*numCells, std::vector<particle_t>(maxNumParticlesinBin));
    //    particleBins.reserve(numCells*numCells);
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
	
    for( int step = 0; step < NSTEPS; step++ )
    {
	navg = 0;
        davg = 0.0;
	dmin = 1.0;
	//
        //  compute forces
        //
	computeForces( particleBins, numCells,  dmin, davg, navg);
        //
        //  move particles
        //
	std::vector<std::vector<int> > leftCell;
	move( particleBins,  leftCell, numCells, h);
	//
        //  transfer particles between cells
        //	
	particleArithmetic( particleBins,  leftCell, numCells, numDeletedParticles);
	//
        //  copy and delete false particles if too many false particles
        //
	

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
