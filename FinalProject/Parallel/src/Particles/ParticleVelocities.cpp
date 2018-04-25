#include "ParticleVelocities.H"

ParticleVelocities::ParticleVelocities()
{
}
ParticleVelocities::ParticleVelocities(ParticleSet& a_state)
{
  // m_box = a_state.m_box;
  // m_dx = a_state.m_dx;
  // m_supportSize = a_state.m_W.supportSize();
  // cout<<"m_box print: ";
  // m_box.print();
  // for (int j = 0; j<DIM; j++)
  //   {
  //     //m_domainSize[j] = m_box.getHighCorner()[0] + 1;
  //     m_domainSize[j] = m_box.getHighCorner()[0];
  //     cout<<"domain size = " <<m_domainSize[j]<<endl; 
  //   }
  
  // for (int k = 0; k < DIM; k++)
  //   {
  //     m_bdry[2*k] = m_box.shift(getUnitv(k)*(-1*m_supportSize))
  //       &m_box.shift(getUnitv(k)*(-1*m_domainSize[k]));
  //     m_bdry[2*k+1] = m_box.shift(getUnitv(k)*(m_supportSize))
  //       &m_box.shift(getUnitv(k)*m_domainSize[k]);
  //   }
  // int m = log2(m_box.getHighCorner()[0]);
  // PS.define(m_dx, m, m_box);
  define(a_state);
}
void ParticleVelocities::define(ParticleSet& a_state)
{
  m_box = a_state.m_box;
  m_dx = a_state.m_dx;
  m_L = a_state.m_L;
  m_supportSize = max(a_state.m_W.supportSize(), 2);
  //cout<<"m_box print: ";
  //m_box.print();
  for (int j = 0; j<DIM; j++)
    {
      //correct for ghost stuff to work.
      m_domainSize[j] = m_box.getHighCorner()[0] + 1;
      //m_domainSize[j] = m_box.getHighCorner()[0];
      //cout<<"domain size = " <<m_domainSize[j]<<endl; 
    }
  
  for (int k = 0; k < DIM; k++)
    {
      m_bdry[2*k] = m_box.shift(getUnitv(k)*(-1*m_supportSize))
        &m_box.shift(getUnitv(k)*(-1*m_domainSize[k]));
      m_bdry[2*k+1] = m_box.shift(getUnitv(k)*(m_supportSize))
        &m_box.shift(getUnitv(k)*m_domainSize[k]);
    }
  //Need the plus 1 to be correct.
  int m = log2(m_box.getHighCorner()[0] + 1);
  //cout<<"In PV  m_box.getHighCorner = " <<m_box.getHighCorner()[0]<<endl;
  //cout<<"In PV define m = "<< m<<endl;
  PS.define(m_dx, m, m_L, m_box);
}
void ParticleVelocities::getGhostDeposition(RectMDArray<double>& enlargedGrid)
{
  // cout<<"m_domainSize =";
  // for (int dir = 0; dir < DIM; dir++)
  //   cout<<m_domainSize[dir];
  // cout<<endl;

  //cout<<"At getGhostDeposition"<<endl;
  for (int k = 0; k < 2*DIM; k++)
    {
      DBox bx = m_bdry[k];
      //bx.print();
      for (Point pt=bx.getLowCorner(); bx.notDone(pt);bx.increment(pt))
        {
          int image[DIM];
          for (int dir = 0; dir < DIM; dir++)
            {
              image[dir] = (pt[dir] + m_domainSize[dir])%m_domainSize[dir];
            }
          Point ptimage(image);
	  // cout<<"pt =";
	  // pt.print();
	  // cout<<"ptimage =";
	  // ptimage.print();
	  // cout<<"value ="<< enlargedGrid[pt]<<endl;
          enlargedGrid[ptimage] += enlargedGrid[pt];
        }
    }
}
// void ParticleVelocities::setGhost(RectMDArray<double>& enlargedGrid)
// {
//   for (int k = 0; k < 2*DIM; k++)
//     {
//       DBox bx = m_bdry[k];
//       for (Point pt=bx.getLowCorner(); bx.notDone(pt);bx.increment(pt))
//         {
//           int image[DIM];
//           for (int dir = 0; dir < DIM; dir++)
//             {
//               image[dir] = (pt[dir] + m_domainSize[dir])%m_domainSize[dir];
//             }
//           Point ptimage(image);
// 	  enlargedGrid[ptimage] = enlargedGrid[pt];
//         }
//     }
// }
void ParticleVelocities::setGhost(RectMDArray<double>& enlargedGrid)
{
  for (int k = 0; k < 2*DIM; k++)
    {
      DBox bx = m_bdry[k];
      for (Point pt=bx.getLowCorner(); bx.notDone(pt);bx.increment(pt))
        {
          int image[DIM];
          for (int dir = 0; dir < DIM; dir++)
            {
              image[dir] = (pt[dir] + m_domainSize[dir])%m_domainSize[dir];
            }
          Point ptimage(image);
	  for (int j=0; j<DIM; j++)
	    {
	      enlargedGrid(pt, j) = enlargedGrid(ptimage, j);
	    }
        }
    }
}
void ParticleVelocities::setGhostMD(RectMDArray<double, DIM>& enlargedGrid)
{
  for (int k = 0; k < 2*DIM; k++)
    {
      DBox bx = m_bdry[k];
      for (Point pt=bx.getLowCorner(); bx.notDone(pt);bx.increment(pt))
        {
          int image[DIM];
          for (int dir = 0; dir < DIM; dir++)
            {
              image[dir] = (pt[dir] + m_domainSize[dir])%m_domainSize[dir];
            }
          Point ptimage(image);
	  for (int j=0; j<DIM; j++)
	    {
	      enlargedGrid(pt, j) = enlargedGrid(ptimage, j);
	    }
        }
    }
}
void ParticleVelocities::operator()(ParticleShift& a_k, 
                     const double& a_time, const double& dt, 
                     ParticleSet& a_state)
{
  //First need to compute a_state.m_x + dt*a_state.m_v + a_k
  //dt is used to control how much of the velocity is added since it isn't used anywhere else
  vector<Particle> t_particles = a_state.m_particles;
  for (int j = 0; j<a_k.m_particles.size(); j++)
    {
      t_particles[j].addVelocity(dt);
      t_particles[j].increment(a_k.m_particles[j]);
    }
  a_state.wrapParticles(t_particles);
  //Then use interpolation made up of the interpolating function given
  RectMDArray<double> density(a_state.m_box.grow(m_supportSize));
  density.setVal(0.0);
  DBox phiBox = a_state.m_box;
  RectMDArray<double> phi(phiBox);
  RectMDArray<double, DIM> EField(a_state.m_box.grow(m_supportSize));
  //phi.setVal(0.0);
  //cout<<"At Deposit"<<endl;
  a_state.deposit(density,  t_particles);

  
  // for (Point p=density.getDBox().getLowCorner(); density.getDBox().notDone(p); density.getDBox().increment(p))
  //   {
  //     p.print();
  //     cout<<density[p]<<endl;
  //   }
  // Deals with Ghost Cells 
  getGhostDeposition(density);
  //setGhost(density);
  // Solve Poisson's Equation with Periodic Boundary Conditions 
  //write density into phu
  for (Point p=phiBox.getLowCorner(); phiBox.notDone(p); phiBox.increment(p))
    {
      phi[p] = density[p];
      // cout<<"phi[ ";
      // p.print();
      // cout<< "]"<<phi[p]<<endl;
    }

   // MPI Communication
  cout<<"Made it to Allreduce"<<endl;
  double temp_array1[phi.dataSize()];
  double temp_array2[phi.dataSize()];
  for (Point p=phiBox.getLowCorner(); phiBox.notDone(p); phiBox.increment(p))
    {
  //      cout<<"Linear index"<<phi.getDBox().getIndex(p)<<endl;
        temp_array2[phi.getDBox().getIndex(p)] = phi[p]; 
    }
  cout<<"Phi Data Size = "<<phi.dataSize()<<endl;
  cout<< "dbox high corner = "<<endl;
  phi.getDBox().getHighCorner().print();
  MPI_Allreduce(&temp_array2, &temp_array1, phi.dataSize(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  //Could potentially replace this with a fucntion redefining the pointer for RectMDArray
  cout<<"Made it past Allreduce"<<endl;
  for (Point p=phiBox.getLowCorner(); phiBox.notDone(p); phiBox.increment(p))
    {
      phi[p] = temp_array1[phi.getDBox().getIndex(p)]; 
    }
  
  
  //Poisson solve
  PS.Solve( phi);
  // cout<<"Made it out of solve"<<endl;
  //     for (Point p=phi.getDBox().getLowCorner(); phi.getDBox().notDone(p); phi.getDBox().increment(p))
  //   {
  //     p.print();
  //     cout<<phi[p]<<endl;
  //   }
  //Write Phi values out back into dbox with ghost cells
  for (Point p=phiBox.getLowCorner(); phiBox.notDone(p); phiBox.increment(p))
  {
   density[p] = phi[p];
   // cout<<"phi[ ";
   // p.print();
   // cout<< "]"<<phi[p]<<endl;
  }
  //Set ghost cells for computing the gradient
  //cout<<"Set Ghost"<<endl;
  setGhost(density);
  // for (Point p=density.getDBox().getLowCorner(); density.getDBox().notDone(p); density.getDBox().increment(p))
  // {
  //   cout<<"density[ ";
  //   p.print();
  //   cout<< "]"<<density[p]<<endl;
  // }
  // Finite Difference 4th order first derivative
  //cout<<"Made it to FD step"<<endl;
  for (Point p=m_box.getLowCorner(); m_box.notDone(p); m_box.increment(p))
    {
      //p.print();
      for (int j = 0; j<DIM; j++)
	{
	  Point ej = getUnitv(j);
	  Point ej2 = ej;
	  ej2 *= 2;
	  //Double check signs
	  EField(p, j)= -(-density[p+ej2] + 8*density[p+ej] - 8*density[p-ej] + density[p- ej2])/(12*m_dx*m_L);
	  // cout<<"EField[ ";
	  // p.print();
	  // cout<< "] ="<<EField(p, j)<<endl;
	      
	}
    }
  setGhostMD(EField);
  //Computes Electric Field Energy for plotting in LLD case.
  int myRank; 
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  if (myRank == 0)
  {
  if (abs(dt) <=pow(10.0, -16))
    {
      double EField_Amplitude = 0.0;
      for (Point p=m_box.getLowCorner(); m_box.notDone(p); m_box.increment(p))
	{
	  for (int j = 0; j<DIM; j++)
	    {
	      EField_Amplitude += pow(EField(p, j), 2);    
	    }
	}
      EField_Amplitude*=m_dx;
      EField_Amplitude = sqrt(EField_Amplitude);
      cout<< EField_Amplitude<<endl;
    }
  }
  //Interpolate back and return particle fields in a_k
  //cout<<"Made it out of FD step"<<endl;
  a_k.zeroEField();
  a_state.InterpolateForce(EField, a_k.m_particles);

  // for (auto iter = a_k.m_particles.begin(); iter!= a_k.m_particles.end(); ++iter)
  //   {
  //     iter->print();
  //   }
  
}
