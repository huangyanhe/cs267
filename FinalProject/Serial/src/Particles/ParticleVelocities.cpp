#include "ParticleVelocities.H"

ParticleVelocities::ParticleVelocities():
  m_box{},
  m_dx{}
{
  m_supportSize = 0.0;
}
ParticleVelocities::ParticleVelocities(ParticleSet& a_state)
{
  m_box = a_state.m_box;
  m_dx = a_state.m_dx;
  m_supportSize = a_state.m_W.supportSize();
  for (int j = 0; j<DIM; j++)
    {
      m_domainSize[j] = m_box.getHighCorner()[0] + 1;
    }
  for (int k = 0; k < DIM; k++)
    {
      m_bdry[2*k] = m_box.shift(getUnitv(k)*(-1*m_supportSize))
        &m_box.shift(getUnitv(k)*(-1*m_domainSize[k]));
      m_bdry[2*k+1] = m_box.shift(getUnitv(k)*(m_supportSize))
        &m_box.shift(getUnitv(k)*m_domainSize[k]);
    }
}
void ParticleVelocities::getGhostDeposition(RectMDArray<double>& enlargedGrid)
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
          enlargedGrid[pt] += enlargedGrid[ptimage];
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
void ParticleVelocities::setGhost(RectMDArray<double, DIM>& enlargedGrid)
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
	      enlargedGrid(ptimage, j) = enlargedGrid(pt, j);
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
  //Then use interpolation made up of the interpolating function given
  RectMDArray<double> density(a_state.m_box.grow(a_state.m_W.supportSize()));
  density.setVal(0.0);
  DBox phiBox = a_state.m_box;
  RectMDArray<double> phi(a_state.m_box.grow(a_state.m_W.supportSize()));
  RectMDArray<double, DIM> EField(a_state.m_box.grow(a_state.m_W.supportSize()));
  phi.setVal(0.0);
  a_state.deposit(density);
  // Deals with Ghost Cells 
  getGhostDeposition(density);
  setGhost(density);
  // Solve Poisson's Equation with Periodic Boundary Conditions 

  
  
  // Finite Difference 4th order first derivative 
  for (Point p=phiBox.getLowCorner(); phiBox.notDone(p); phiBox.increment(p))
        {
	  for (int j = 0; j<DIM; j++)
	    {
	      Point ej = getUnitv(j);
	      Point ej2 = ej;
	      ej2 *= 2;
	      //Double check signs
	      EField(p, j)= -(-phi[p+ej2] + 8*phi[p+ej] - 8*phi[p-ej] + phi[p- ej2])/(12*m_dx);
	    }
	}
  //Interpolate back and return particle fields in a_k
  a_k.setToZero();
  a_state.InterpolateForce(EField, a_k.m_particles);
}
