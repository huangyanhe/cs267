#include "common.h"

void assignCells2P(int numCells, int numProcessors, std::vector<std::vector<int>> assignedCells, std::vector<int> cellProcessorMapping)
{
  int BlockSize = ceil(numCells*numCells/numProcessors);
  int numProcessorPerRow = ceil(numCells/BlockSize);
  int iproc = 0, jproc = 0;
  for (int i = 0; i<numCells; i += BlockSize, iproc++)
    {
      M = min(BlockSize, numCells - i)
	for (int j = 0; j<numCells; j += BlockSize, jproc++)
	{
	  N = min(BlockSize, numCells - j)
	  for (int isub = 0; isub<M; isub++)
	    {
	      for (int jsub = 0; jsub<N; jsub++)
		{
		  assignedCells[iproc + jproc*numProcessorPerRow].push_back((i + isub) + (j + jsub)*numCells);
		  cellProcessorMapping[(i + isub) + (j + jsub)*numCells] = iproc + jproc*numProcessorPerRow;
		}
	    }
	  std::sort(assignedCells[iproc + jproc*numProcessorPerRow].begin(), assignedCells[iproc + jproc*numProcessorPerRow].end()); 
	}
    }
}
//This seems wrong
void findNeighborProcessorsandCells( int rank, int *myCells, std::vector<std::vector<int>> &a_neighbors,std::vector<int> cellProcessorMapping, std::set<int>& neighborProcessors, std::set<int>& neighborCellsfromProcessors, std::set<int> neighborCellstoProcessors)
{
  for (int j = 0; j< myCells.size(); j++)
    {
      for (int k = 0; k < a_neighbors[myCells[i]].size(); k++)
	{
	  if (cellProcessorMapping[a_neighbors[myCells[i]][k]] != rank)
	    neighborProcessors.insert(cellProcessorMapping[a_neighbors[myCells[i]][k]]);
       }
    }
  //std::set<int>::iterator iter; 
  //iter = std::set_difference(neighborProcessorsTemp.begin(), neighborProcessorsTemp.end(), myCells.begin(), myCells.end(), neighborProcessors.begin());
  //neighborProcessors.resize( iter - neighborProcessors.begin());
  //Computes who I get cells from
  //std::vector<int> neighborCellsfromProcessors; 
  for (int j = 0; j< myCells.size(); j++)
    {
      std::vector<int>::iterator iter;
      iter = std::set_difference(a_neighbors[myCells[j]].begin(), a_neighbors[myCells[j]].end(), myCells.begin(),  myCells.end(), neighborCellsfromProcessors.end());
      //      neighborCellsfromProcessors.resize(iter - neighborCellsfromProcessors.begin());
    }
  //Computes who I send cells to
  for(auto it neighborProcessors.begin(); it != neighborProcessors.end(); ++it)
    {
      std::vector<int> neighborCellstoProcessors; 
      std::vector<int>::iterator iter;
      for (int r =0; r< a_neighbors[assignedCells[*it]].end(), r++ )
	{
	  iter = std::set_intersection(a_neighbors[assignedCells[*it][r]].begin(), a_neighbors[assignedCells[*it][r]].end(), myCells.begin(),  myCells.end(), neighborCellstoProcessors.end());
	  //	  neighborCellstoProcessors.resize(iter - neighborCellstoProcessors.begin());
	}
      
    }
  
}
void oneDtwoDMapping(int numCells, std::vector<vector<int>> cellIndexMapping)
{
  for (int i =0; i< numCells; i++ )
    {
      for (int j =0; j< numCells; j++ )
	{
	  cellIndexMapping[i + numCells*j][0] = i;
	  cellIndexMapping[i + numCells*j][1] = j;
	}
    }
}
void buildNeighbors(std::vector<std::vector<int>> &a_neighbors, int numCells)
{
  if (numCells == 1)
    {
      std::vector<int> a_temp(numCells*numCells);
      for (int i = 0; i< numCells; i++)
	{
	  for (int j =0; j<numCells; j++ )
	    {
	      a_temp[i + numCells*j] = i + numCells*j;
	    }
	}
      a_neighbors.assign(numCells*numCells, a_temp);
      return;
    }
  
  for (int i = 0; i< numCells; i++)
    {
      for (int j =0; j<numCells; j++ )
	{
	  int index = i+ numCells*j;
	  a_neighbors[i + numCells*j].push_back(i + numCells*j);
	  if (i == 0)
	    {
	      if (j == 0)
		{
		  a_neighbors[index].push_back(i + numCells*(j+1));
		  a_neighbors[index].push_back(i+1 + numCells*(j+1));
		  a_neighbors[index].push_back(i+1 + numCells*j);
		}
	      else if (j == numCells-1)
		{
		  a_neighbors[index].push_back(i + numCells*(j-1));
		  a_neighbors[index].push_back(i+1 + numCells*(j-1));
		  a_neighbors[index].push_back(i+1 + numCells*j);
		}
	      else
		{
		  a_neighbors[index].push_back(i + numCells*(j+1));
		  a_neighbors[index].push_back(i + numCells*(j-1));
		  a_neighbors[index].push_back(i+1 + numCells*(j+1));
		  a_neighbors[index].push_back(i+1 + numCells*j);
		  a_neighbors[index].push_back(i+1 + numCells*(j-1));
		}
	    }
	  else if ( i == numCells-1)
	    {
	      if (j == 0)
		{
		  a_neighbors[index].push_back(i + numCells*(j+1));
		  a_neighbors[index].push_back(i-1 + numCells*(j+1));
		  a_neighbors[index].push_back(i-1 + numCells*j);
		}
	      else if (j == numCells-1)
		{
		  a_neighbors[index].push_back(i + numCells*(j-1));
		  a_neighbors[index].push_back(i-1 + numCells*(j-1));
		  a_neighbors[index].push_back(i-1 + numCells*j);
		}
	      else
		{
		  a_neighbors[index].push_back(i + numCells*(j+1));
		  a_neighbors[index].push_back(i + numCells*(j-1));
		  a_neighbors[index].push_back(i-1 + numCells*(j+1));
		  a_neighbors[index].push_back(i-1 + numCells*j);
		  a_neighbors[index].push_back(i-1 + numCells*(j-1));
		}
	    }
	  else if (j == 0)
	    {
	      a_neighbors[index].push_back(i + numCells*(j+1));
	      a_neighbors[index].push_back(i+1 + numCells*(j+1));
	      a_neighbors[index].push_back(i-1 + numCells*(j+1));
	      a_neighbors[index].push_back(i+1 + numCells*j);
	      a_neighbors[index].push_back(i-1 + numCells*j);
	    }
	  else if (j == numCells-1)
	    {
	      a_neighbors[index].push_back(i + numCells*(j-1));
	      a_neighbors[index].push_back(i+1 + numCells*(j-1));
	      a_neighbors[index].push_back(i-1 + numCells*(j-1));
	      a_neighbors[index].push_back(i+1 + numCells*j);
	      a_neighbors[index].push_back(i-1 + numCells*j);
	    }
	  else
	    {
	      a_neighbors[index].push_back(i + numCells*(j-1));
	      a_neighbors[index].push_back(i+1 + numCells*(j-1));
	      a_neighbors[index].push_back(i-1 + numCells*(j-1));
	      a_neighbors[index].push_back(i+1 + numCells*j);
	      a_neighbors[index].push_back(i-1 + numCells*j);
	      a_neighbors[index].push_back(i + numCells*(j+1));
	      a_neighbors[index].push_back(i+1 + numCells*(j+1));
	      a_neighbors[index].push_back(i-1 + numCells*(j+1));
	    }
	}
    }
}
void bin( std::vector< std::vector<particle_wrapper> >& particleBins, particle_t *particles, int n, int numCells,double h)
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
      //bool doesExist = particles[p].exists;
      particleBins[iposx + iposy*numCells].push_back();
      particleBins[iposx + iposy*numCells].end().particle = particles[p];
      particleBins[iposx + iposy*numCells].end().exists = true; 
    }
};

void moveandcheck(std::vector< std::vector<particle_t> >& particleBins,  int *myCells, std::vector<std::vector<int> >& cellIndexMapping, std::vector<std::vector<int> >& leftCell, int numCells, double h)
{
  int iposx, iposy, ix, jy;
  for( int i = 0; i< myCells.size(); i++)
    {
	  for(int k = 0; k<particleBins[myCells[i]].size(); k++ )
	    {
	      bool thisParticleExists = particleBins[myCells[i]][k].exists;
	      if (particleBins[myCells[i]][k].exists == true)
		{
		  move( particleBins[myCells[i]][k]);
		  // set acceleration to zero
		  particleBins[myCells[i]][k].ax = 0;
		  particleBins[myCells[i]][k].ay = 0;
		  iposx = floor(particleBins[myCells[i]][k].x/h);
		  iposy = floor(particleBins[myCells[i]][k].y/h);
		  ix = cellIndexMapping[myCells[i]][0];
		  jy = cellIndexMapping[myCells[i]][1];
		  if (iposx != i || iposy != j)
		    {
		      std::vector<int> Coordinates(5, 0);
		      Coordinates[0] = ix;
		      Coordinates[1] = iy;		      
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
void particleArithmetic(std::vector< std::vector<particle_wrapper> >& particleBins, std::vector<std::vector<int> >& leftCell, int *myCells, std::vector<int>& cellProcessorMapping, std::set<int>& neighborProcessors, vector< transfered_particle>& particles2Transfer,int numCells, int numDeletedParticles)
{
  for(int j = 0; j<leftCell.size(); j++)
    {
      //Particle stays in same processor.
      if ( std::any_of( myCells.begin(), myCells.end(), [](int a) {return a ==leftCell[j][3] + leftCell[j][4]*numCells}) )
	{
      //numDeletedParticles += 1;
	  particleBins[leftCell[j][3] + leftCell[j][4]*numCells].push_back(particleBins[leftCell[j][0] + leftCell[j][1]*numCells][leftCell[j][2]]);
	  particleBins[leftCell[j][0] + leftCell[j][1]*numCells][leftCell[j][2]].exists = false;
	}
      else
	{
	  int newCellIndex, newProcessor;
	  newCellIndex = leftCell[j][3] + leftCell[j][4]*numCells;
	  newProcessor = cellProcessorMapping[newCellIndex];
	  if ( std::any_of( neighborProcessors.begin(), neighborProcessors.end(), [](int a) {return a == newProcessor}) )
	    {
	      particles2Transfer.push_back();
	      particles2Transfer.end().x = particleBins[leftCell[j][0] + leftCell[j][1]*numCells][leftCell[j][2]].particle.x;
	      particles2Transfer.end().y = particleBins[leftCell[j][0] + leftCell[j][1]*numCells][leftCell[j][2]].particle.y;
	      particles2Transfer.end().vx = particleBins[leftCell[j][0] + leftCell[j][1]*numCells][leftCell[j][2]].particle.vx;
	      particles2Transfer.end().vy = particleBins[leftCell[j][0] + leftCell[j][1]*numCells][leftCell[j][2]].particle.vy;
	      particles2Transfer.end().bin = newCellIndex;
	      particles2Transfer.end().processor = newProcessor;
	    }
	  else
	    {
	      printf ("\n A Particle traveled more than a processor.");
	      fflush(stdout);
	    }
	}
    }
}
void unwrapReceivedParticles(std::vector< std::vector<particle_wrapper> >& particleBins, vector< particle_t> particlesFromProcessorj, int numCells, double h)
{
  for (int j = 0; j < particlesFromProcessorj.size(); j++ )
    {
      iposx = floor(particlesFromProcessorj[j].x/h);
      iposy = floor(particlesFromProcessorj[j].y/h);
      //bool doesExist = particles[p].exists;
      particleBins[iposx + iposy*numCells].push_back();
      particleBins[iposx + iposy*numCells].end().particle = particlesFromProcessorj[j];
      particleBins[iposx + iposy*numCells].end().exists = true; 
    }
}
void unwrapReceivedParticlesNeighbors(std::vector< std::vector<particle_wrapper> >& particleBins, int *myCells,vector< particle_t> particlesFromProcessorj, int numCells, double h)
{
  for (auto it = neighborCellsfromProcessors.begin(); it != neighborCellsfromProcessors.end(); ++it)
    {
      for (int k = 0; k< particleBins[*it].size(); k++)
	particleBins[*it][k].exists = false;
    }
  for (int j = 0; j < particlesFromProcessorj.size(); j++ )
    {
      iposx = floor(particlesFromProcessorj[j].x/h);
      iposy = floor(particlesFromProcessorj[j].y/h);
      //bool doesExist = particles[p].exists;
      particleBins[iposx + iposy*numCells].push_back();
      particleBins[iposx + iposy*numCells].end().particle = particlesFromProcessorj[j];
      particleBins[iposx + iposy*numCells].end().exists = true;
      if (iposx + numCells*iposy)
    }
}

void computeForces(std::vector< std::vector<particle_t> >& particleBins, std::vector<std::vector<int>> &a_neighbors, int *myCells ,int numCells, double& dmin, double& davg, int& navg)
{
  for( int i = 0; i< myCells.size(); i++)
    {
      for(int q=0; q<a_neighbors[myCells[i]].size(); q++)
	{
	  for(int k = 0; k<particleBins[myCells[i]].size(); k++ )
	    {
	      for(int r = 0; r<particleBins[a_neighbors[myCells[i]][q]].size(); r++ )
		{
		  if (particleBins[myCells[i]][k].exists == true && particleBins[a_neighbors[myCells[i]][q]][r].exists == true)
		    apply_force( particleBins[myCells[i]][k], particleBins[a_neighbors[myCells[i]][q]][r],&dmin,&davg,&navg);

		}
	    }
	}
    }
}
