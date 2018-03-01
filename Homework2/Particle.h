#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include "common.h"
#include <algorithm>
#include <vector>
#include <set>

//
// particle data structure
//
typedef struct 
{
  particle_t particle;
  bool exists;		 
} particle_wrapper;

typedef struct 
{
  std::vector< double>  x;
  std::vector< double>  y;
  std::vector< double>  vx;
  std::vector< double>  vy;
  int bin; 
  int processor;		 
} transfered_particle;

/* typedef struct  */
/* { */
/*   int processor; */
/*   vector<int> myCells */
/*   std::vector<std::vector<int> > assignedCells; */
/*   std::set<int> neighborProcessors; */
/*   std::vector<int> neighborCellsfromProcessors; */
/*   std::vector<int> neighborCellstoProcessors; */
  		 
/* } Processor; */

auto removeParticle = [](particle_t particle) -> bool
{
  bool removeIfFalse =  !particle.exists;
  return removeIfFalse;
};

void assignCells2P(int numCells, int numProcessors, std::vector<std::vector<int>> assignedCells, std::vector<int> cellProcessorMapping);
void findNeighborProcessorsandCells( int rank, int *myCells, std::vector<std::vector<int>> &a_neighbors,std::vector<int> cellProcessorMapping, std::set<int>& neighborProcessors, std::set<int>& neighborCellsfromProcessors, std::set<int> neighborCellstoProcessors);


void oneDtwoDMapping(int numCells, std::vector<vector<int>> cellIndexMapping);


void buildNeighbors(std::vector<std::vector<int>> &a_neighbors, int numCells);
	

void bin( std::vector< std::vector<particle_wrapper> >& particleBins, particle_t *particles, int n, int numCells,double h);

void moveandcheck(std::vector< std::vector<particle_wrapper> >& particleBins,  int *myCells, std::vector<std::vector<int> >& cellIndexMapping, std::vector<std::vector<int> >& leftCell, int numCells, double h)
/* Need to think about what the correct thing to do in this case is. One option is to delete the particle but maybe better is to leave it alone and just flag it as not used. I could move it outside the domain and anytime I came to the particle for force calculations, I could just ignore it. 
 */ 
void particleArithmetic(std::vector< std::vector<particle_wrapper> >& particleBins, std::vector<std::vector<int> >& leftCell, int *myCells, std::vector<int>& cellProcessorMapping, std::set<int>& neighborProcessors, vector< transfered_particle>& particles2Transfer,int numCells, int numDeletedParticles);

void computeForces(std::vector< std::vector<particle_t> >& particleBins, std::vector<std::vector<int>> &a_neighbors, int *myCells ,int numCells, double& dmin, double& davg, int& navg)

void unwrapReceivedParticles(std::vector< std::vector<particle_wrapper> >& particleBins, vector< particle_t> particlesFromProcessorj, int numCells,double h);
void unwrapReceivedParticlesNeighbors(std::vector< std::vector<particle_wrapper> >& particleBins, int *myCells,vector< particle_t> particlesFromProcessorj, int numCells, double h);



#endif

