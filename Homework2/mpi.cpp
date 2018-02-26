#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg; 


    double density = 0.0005;
    double cutoff = 0.01;
    double h = cutoff;
    double maxXY = sqrt(n*density);
    int numCells = ceil(maxXY/h);
    int numDeletedParticles = 0;
    int copyFactor = 10; 
    int avgParticlesPerCell = ceil(n/numCells);
    int overStorageFactor = 2;
    int maxNumParticlesinBin = avgParticlesPerCell*overStorageFactor;
    //std::vector<particle_t> initVector;
    //initVector.reserve(maxNumParticlesinBin);

    
    
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
    MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
    
    std::vector< std::vector<particle_t> > particleBins;
    bin( particleBins, particles, n, numCells, h);
    std::vector<std::vector<int>> neighbors;
    neighbors.reserve(numCells*numCells);
    buildNeighbors(neighbors, numCells);
    std::vector<vector<int>> cellIndexMapping;
    cellIndexMapping.reserve(numCells*numCells);
    oneDtwoDMapping( numCells, cellIndexMapping);
    std::vector<std::vector<int>> assignedCells;
    assignedCells.reserve(n_proc);
    std::vector<int> cellProcessorMapping;
    cellProcessorMapping.reserve(numCells*numCells);
    assignCells2P(numCells, n_proc, assignedCells, cellProcessorMapping);

    int BlockSize = ceil(numCells*numCells/n_proc);
    
    if( rank == 0 )
      {
        init_particles( n, particles );
	for (int j = 1; j<n_proc; j++)
	  {
	    MPI_Send(&assignedCells[j].front(), &assignedCells[j].size(), MPI_INT, j, MPI_COMM_WORLD);
	  }
      }
    else
      {
	//Buffer the size to be half a size larger just in case.
	int myCells[BlockSize + BlockSize/2];
	MPI_IRecv(&myCells, BlockSize + BlockSize/2, MPI_INT, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    
    std::set<int> neighborProcessors;
    findNeighborProcessors( myCells, neighbors, cellProcessorMapping, neighborProcessors);

    
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
        //  MPI_Allgatherv( local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );
        
        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if( find_option( argc, argv, "-no" ) == -1 )
          if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );

	//
        //  compute forces
        //
	computeForces( particleBins, neighbors, myCells, numCells,  dmin, davg, navg);
        //
        //  move particles
        //
	std::vector<std::vector<int> > leftCell;
	moveandcheck( particleBins, myCells, cellIndexMapping, leftCell, numCells, h);
	//
        //  transfer particles between cells
        //
	numDeletedParticles += leftCell.size();
	vector< transfered_particle> particles2Transfer;
	particleArithmetic( particleBins,  leftCell, myCells, cellProcessorMapping, neighborProcessors, particles2Transfer,numCells, numDeletedParticles);
	//
        //  Communication of Particles that changed cells and processors
        //
	for (auto it = neighborProcessors.begin(); it != neighborProcessors.end(); ++it)
	  {
	    vector< particle_t> particles4Processorj;
	    for (int k=0; k < particles2Transfer.size(); k++ )
	      {
		if (particles2Transfer[k].processor == *it)
		  {
		    particles4Processorj.push_back();
		    particles4Processorj.end().x = particles2Transfer[k].x;
		    particles4Processorj.end().y = particles2Transfer[k].y;
		    particles4Processorj.end().vx = particles2Transfer[k].vx;
		    particles4Processorj.end().vy = particles2Transfer[k].vy;
		    particles4Processorj.end().ax = 0.0;
		    particles4Processorj.end().ay = 0.0;
		  }
	      }
	    if (avgParticlesPerCell < particles4Processorj.size())
	      {
		printf ("\n Particles sent was longer than particles received.");
		fflush(stdout);
	      }
	    MPI_Isend(&particles4Processorj.begin(), &particles4Processorj.size(), PARTICLE, *it, step, MPI_COMM_WORLD);
	  }
	vector< particle_t> particlesFromProcessorj;
	for (auto it = neighborProcessors.begin(); it != neighborProcessors.begin(); ++it)
	  {
	    vector< particle_t> particlesFromProcessorjTemp;
	    // This could be dangerous but can just empirically check that it is okay.
	    particlesFromProcessorjTemp.resize(avgParticlesPerCell);
	    MPI_IRecv(&particlesFromProcessorjTemp, particlesFromProcessorjTemp.end(), PARTICLE, *it, step,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    particlesFromProcessorj.insert( particlesFromProcessorj.end(), particlesFromProcessorjTemp.begin(), particlesFromProcessorjTemp.end());
	  }

	unwrapReceivedParticles( particleBins, particlesFromProcessorj, h);

	//
        //  Communication of Particles on boundary needed for next apply force.
        //
	//Does the send of all particles that processor touches to all neighboring processors (not finished)

	vector< particle_t> particles4ProcessorNeighbors;
	for (auto it = neighborCellstoProcessors.begin(); it != neighborCellstoProcessors.end(); ++it)
	  {
	    for( int k = 0; k<particleBin[*it].end(); k++)
	      {
		particles4ProcessorNeighbors.push_back();
		particles4ProcessorNeighbors.end().x =particleBin[*it][k].particle.x;
		particles4ProcessorNeighbors.end().y = particleBin[*it][k].particle.y;
		particles4ProcessorNeighbors.end().vx = particleBin[*it][k].particle.vx;
		particles4ProcessorNeighbors.end().vy = particleBin[*it][k].particle.vy;
		particles4ProcessorNeighbors.end().ax = 0.0;
		particles4ProcessorNeighbors.end().ay = 0.0;
	      }
	  }
	for (auto it = neighborProcessors.begin(); it != neighborProcessors.end(); ++it)
	  {
	    MPI_Isend(&particles4ProcessorNeighbors.begin(), &particles4ProcessorNeighbors.size(), PARTICLE, *it, step, MPI_COMM_WORLD);
	  }
	//Does Receive (probably working)
	vector< particle_t> particlesFromProcessorjNeighbors;
	for (auto it = neighborCellsfromProcessors.begin(); it != neighborCellsfromProcessors.end(); ++it)
	  {
	    vector< particle_t> particlesFromProcessorjTemp;
	    // This could be dangerous but can just empirically check that it is okay.
	    particlesFromProcessorjTemp.resize(3*avgParticlesPerCell);
	    MPI_IRecv(&particlesFromProcessorjTemp, particlesFromProcessorjTemp.end(), PARTICLE, cellProcessorMapping[*it], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    particlesFromProcessorjNeighbors.insert( particlesFromProcessorjNeighbors.end(), particlesFromProcessorjTemp.begin(), particlesFromProcessorjTemp.end());
	  }
	

	//Need a Waitall statement somewhere around here. 
	  
	
	
         copy and delete false particles if too many false particles
        
	if (step%10 == 0)
	  {
	    for (int j =0; j<numCells*numCells; j++)
	      {
		particleBins[j].erase(std::remove_if(particleBins[j].begin(), particleBins[j].end(), removeParticle), end(particleBins[j])); 
	      }
	    //numDeletedParticles = 0;
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
    free( myCells );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
