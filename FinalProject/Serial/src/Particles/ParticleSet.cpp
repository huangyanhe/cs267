#include "ParticleSet.H"

//Constructors
ParticleSet::ParticleSet(
              shared_ptr<ConvKernel>& a_kerptr,
              DBox& a_box,
              double& a_dx, 
              array<double, DIM>& a_lowCorner,
              int a_M, int order, int smoothness):
  m_particles{},
  m_dx{a_dx},
  m_box{a_box},
  m_lowCorner{a_lowCorner},
  m_W{order, smoothness}
{
  m_hockney.define(a_kerptr, a_dx, a_M);
  for (int k = 0; k < DIM; k++)
    {
      m_bdry[2*k] = m_box.shift(getUnitv(k)*(-1))
        &m_box.shift(getUnitv(k)*(-m_domainSize));
      m_bdry[2*k+1] = m_box.shift(getUnitv(k))
        &m_box.shift(getUnitv(k)*m_domainSize);
    }
}
ParticleSet::ParticleSet(
	      DBox& a_box,
              double& a_dx, 
              array<double, DIM>& a_lowCorner,
              int a_M, int order, int smoothness):
  m_particles{},
  m_dx{a_dx},
  m_box{a_box},
  m_lowCorner{a_lowCorner},
  m_W{order, smoothness}
{
  for (int k = 0; k < DIM; k++)
    {
      m_bdry[2*k] = m_box.shift(getUnitv(k)*(-1))
        &m_box.shift(getUnitv(k)*(-m_domainSize));
      m_bdry[2*k+1] = m_box.shift(getUnitv(k))
        &m_box.shift(getUnitv(k)*m_domainSize);
    }
}
// Functions
void ParticleSet::increment(const ParticleShift& a_shift)
{
  for (int j = 0; j<a_shift.m_particles.size(); j++)
    {
      m_particles[j].increment(a_shift.m_particles[j]);
    }
}
// Dimensionally Independent Deposition
void  ParticleSet::deposit(RectMDArray<double>& a_Charge)
{
  array<double,DIM> pos;
  double interpcoeff;
  for (int it = 0; it <m_particles.size(); it++)
    {
      array<int,DIM> iposLow, iposHigh;
      for (int l = 0; l < DIM; l++)
        {
	  pos[l] = m_particles[it].m_x[l];
	  iposLow[l] = floor(pos[l]/m_dx);
	  iposHigh[l] = ceil(pos[l]/m_dx); 
	}

      Point Shift = getUnitv(0);
      Shift *= 0;
      for (int j =0; j<DIM; j++)
	{
	  Shift += getUnitv(j);
	}
      Shift *= (m_W.supportSize() - 1);
      Point HighCorner(iposHigh);
      Point LowCorner(iposLow);
      Point LC = LowCorner - Shift;
      Point HC = HighCorner+ Shift;
      DBox SupportBox(LC, HC); 
      interpcoeff = 1/m_dx*m_particles[it].strength;
      for (Point s = SupportBox.getLowCorner(); SupportBox.notDone(s); SupportBox.increment(s))
	{
	  double KernelProduct = 1.0;
	  for (int j =0; j<DIM; j++)
	    {
	      int region = min(abs(s[j] - LowCorner[j]), abs(s[j] - HighCorner[j]));
	      double val = abs(s[j]*m_dx - pos[j])/m_dx; 
	      KernelProduct *= m_W.apply(val, region);
	    }
	  a_Charge[s] += interpcoeff*KernelProduct;
	}
    }
}
//Dimensionally Independent Deposition
void  ParticleSet::InterpolateForce(RectMDArray<double>& a_Field)
{
  
  for (int it = 0; it <m_particles.size(); it++)
    {
      array<int,DIM> iposLow, iposHigh;
      array<double,DIM> pos;
      for (int l = 0; l < DIM; l++)
        {
	  pos[l] = m_particles[it].m_x[l];
	  iposLow[l] = floor(pos[l]/m_dx);
	  iposHigh[l] = ceil(pos[l]/m_dx); 
	}

      Point Shift = getUnitv(0);
      Shift *= 0;
      for (int j =0; j<DIM; j++)
	{
	  Shift += getUnitv(j);
	}
      Shift *= (m_W.supportSize() - 1);
      Point HighCorner(iposHigh);
      Point LowCorner(iposLow);
      Point LC = LowCorner - Shift;
      Point HC = HighCorner+ Shift;
      DBox SupportBox(LC, HC); 
      for (Point s = SupportBox.getLowCorner(); SupportBox.notDone(s); SupportBox.increment(s))
	{
	  double KernelProduct = 1.0;
	  for (int j =0; j<DIM; j++)
	    {
	      int region = min(abs(s[j] - LowCorner[j]), abs(s[j] - HighCorner[j]));
	      double val = abs(s[j]*m_dx - pos[j])/m_dx; 
	      KernelProduct *= m_W.apply(val, region);
	    }
	  for (int j = 0; j<DIM; j++)
	    {
	      m_particles[it].EField[j] += a_Field(s, j)*KernelProduct;
	    }
	}
    }
}
// This function needs to be done carefully and might require changing how the boundary is computed.
void ParticleSet::getGhost(RectMDArray<double >& a_phi )
{
  for (int k = 0; k < 2*DIM; k+=2)
    {
      DBox bx = m_bdry[k];
      for (Point pt=bx.getLowCorner(); bx.notDone(pt);bx.increment(pt))
        {
          int image[DIM];
          for (int dir = 0; dir < DIM; dir++)
            {
              image[dir] = (pt[dir] + m_domainSize)%m_domainSize;
            }
          Point ptimage(image);
          a_phi[pt] += a_phi[ptimage];
	  a_phi[ptimage] = a_phi[pt];
        }
    }
};
