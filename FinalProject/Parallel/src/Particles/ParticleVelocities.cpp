#include "ParticleVelocities.H"
#include "omp.h"
#include <stdlib.h>
#include <stdio.h>
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

// void ParticleVelocities::getGhostDeposition(RectMDArray<double>& enlargedGrid)
// {
//   // cout<<"m_domainSize =";
//   // for (int dir = 0; dir < DIM; dir++)
//   //   cout<<m_domainSize[dir];
//   // cout<<endl;
 
//   //cout<<"At getGhostDeposition"<<endl;
//   //cout << "m_box = " << m_box.sizeOf() << endl;

  
//   int l = enlargedGrid.dataSize();
//   omp_lock_t* box_lock = (omp_lock_t*)malloc(l*sizeof(omp_lock_t));
//   for(int i = 0; i < l; i++)
//     {
//       omp_init_lock(&box_lock[i]);
//     }
// #pragma omp parallel shared(enlargedGrid, box_lock) 
//   {
//     for (int k = 0; k < 2*DIM; k++)
//       {
// 	DBox bx = m_bdry[k];
// 	//bx.print();
// 	//      for (Point pt=bx.getLowCorner(); bx.notDone(pt);bx.increment(pt))
// #pragma omp for   
// 	for (int i = 0; i < bx.sizeOf(); i++)
// 	  {
// 	    Point pt = bx.getPoint(i);
// 	    int image[DIM];
// 	    for (int dir = 0; dir < DIM; dir++)
// 	      {
// 		image[dir] = (pt[dir] + m_domainSize[dir])%m_domainSize[dir];
// 	      }
// 	    Point ptimage(image);
// 	    // cout<<"pt =";
// 	    // pt.print();
// 	    // cout<<"ptimage =";
// 	    // ptimage.print();
// 	    // cout<<"value ="<< enlargedGrid[pt]<<endl;
// 	    int j = m_box.getIndex(ptimage);
// 	    //cout << "ptimage= " << j << "m_box= " << m_box.sizeOf() << endl;
	    
// 	    //int num = omp_get_num_threads();
// 	    //cout << num << endl;
// 	    //int ID = omp_get_thread_num();
	    
// 	    //	  cout << "test locks =" << omp_test_lock(&box_lock[j]) << endl;
// 	    omp_set_lock(&box_lock[j]);
// 	    //cout <<" lock using= " << j << endl;
// 	    enlargedGrid[ptimage] += enlargedGrid[pt];
// 	    omp_unset_lock(&box_lock[j]);
// 	  }
//       }
//   }
//   for (int i = 0; i < l; i++)
//     {
//       omp_destroy_lock(&box_lock[i]);
//     }
//   free(box_lock);
// }
// Void
// Particlevelocities
//::setGhost(RectMDArray<double>& enlargedGrid)
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
#pragma omp parallel for 
      for (int i = 0; i < bx.sizeOf(); i++)
        {
          Point pt = bx.getPoint(i);
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
#pragma omp parallel for 
      for (int i = 0; i < bx.sizeOf(); i++)
        {
          Point pt = bx.getPoint(i);
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
#pragma omp parallel for schedule(dynamic)
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
  //cout << "stopped here" << endl;
  getGhostDeposition(density);
  // cout << "finished once" << endl;
  //setGhost(density);
  // Solve Poisson's Equation with Periodic Boundary Conditions 
  //write density into phu
  int l = phiBox.sizeOf();
#pragma omp parallel for
  for (int i = 0; i < l; i++)
    {
      Point p = phiBox.getPoint(i);
      phi[p] = density[p];
      // cout<<"phi[ ";
      //p.print();
      //cout<< "]"<<phi[p]<<endl;
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
#pragma omp parallel for
  for (int i = 0; i < l; i++)
  {
    Point p = phiBox.getPoint(i);
    density[p] = phi[p];
   // cout<<"phi[ ";
   // p.print();
   // cout<< "]"<<phi[p]<<endl;
  }
  //Set ghost cells for computing the gradient
  //cout<<"Set Ghost"<<endl;
  //cout << "setGhost start" << endl;
  setGhost(density);
  //cout << "setGhost finished" << endl;
  // for (Point p=density.getDBox().getLowCorner(); density.getDBox().notDone(p); density.getDBox().increment(p))
  // {
  //   cout<<"density[ ";
  //   p.print();
  //   cout<< "]"<<density[p]<<endl;
  // }
  // Finite Difference 4th order first derivative
  //cout<<"Made it to FD step"<<endl;
  int r = m_box.sizeOf();
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < r; i++)
    {
      Point p = m_box.getPoint(i);
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

  //cout << "setGhostMD start" << endl;
  setGhostMD(EField);
  //cout << "setGhostMD finish" << endl;
  //Computes Electric Field Energy for plotting in LLD case.
  if (abs(dt) <=pow(10.0, -16))
    {
      double EField_Amplitude = 0.0;
#pragma omp parallel for collapse(2) reduction(+:EField_Amplitude)
      for (int i = 0; i < r; i++)
	{
	  for (int j = 0; j<DIM; j++)
	    {
	      Point p = m_box.getPoint(i);
	      EField_Amplitude += pow(EField(p, j), 2);    
	    }
	}
      EField_Amplitude*=m_dx;
      EField_Amplitude = sqrt(EField_Amplitude);
      //      cout<< " EField_Amplitude = " << EField_Amplitude <<endl;
      cout<<  EField_Amplitude <<endl;
    }
  //Interpolate back and return particle fields in a_k
  //cout<<"Made it out of FD step"<<endl;
  a_k.zeroEField();
  //cout << "InterpolateForce start" << endl;
  a_state.InterpolateForce(EField, a_k.m_particles);
  // cout << "InterpolateForce finish" << endl;
  //  for (auto iter = a_k.m_particles.begin(); iter!= a_k.m_particles.end(); ++iter)
  //   {
  //    iter->print();
  //   }  
}
