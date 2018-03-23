#include "ParticleSet.H"


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
}
void ParticleSet::increment(const ParticleShift& a_shift)
{
  for (int j = 0; j<a_shift.m_particles.size(); j++)
    {
      m_particles[j].increment(a_shift.m_particles[j]);
    }
}
void  ParticleSet::deposit(RectMDArray<double>& a_Charge)
{
  //loop over particles
  //I would like to write this in a dimensionally independent way some day.
  int ipos[2];
  double sk, pos, interpcoeff;
  for (int it = 0; it <m_particles.size(); it++)
    {
      // cout<<"it ="<<it<<endl;
      pos = m_particles[it].m_x[0];
      ipos[0] = floor(pos/m_dx);
      ipos[1] = 0;
      sk = (pos - ipos[0]*m_dx)/m_dx;
      interpcoeff = 1/m_dx*m_particles[it].strength;
      Point e0 = getUnitv(0);
      Point ipoint(ipos);
      Point shift;
      // cout<<"ipoint[0] ="<<ipoint[0]<<", ipoint[1] ="<<ipoint[1]<<endl;
      for (int s = 0; s < m_W.supportSize(); s++)
	{
	  shift = e0;
	  shift *= s;
	  // cout<<"shift[0] ="<<shift[0]<<", shift[1] ="<<shift[1]<<endl;
	  // cout<<"s ="<< s<<endl;
	  // cout<<"interpcoeff ="<< interpcoeff<<",W(sk) = "<< interpcoeff*m_W.apply(sk, s)<<", W(1 - sk) = "<< m_W.apply(1-sk, s)<<endl;
	  a_Charge[ipoint - shift] += interpcoeff*m_W.apply(sk, s);
	  a_Charge[ipoint + e0 + shift] += interpcoeff*m_W.apply(1-sk, s);
	  // cout<<"ChargeL = "<< a_Charge[ipoint - shift]<<", ChargeR = "<<a_Charge[ipoint + e0 + shift]<<endl;
	}
    }
}
void ParticleSet::InterpolateForce(RectMDArray<double>& a_Field)
{
   //loop over particles
  //I would like to write this in a dimensionally independent way some day.
  int ipos[2];
  double sk, pos;
  for (int it = 0; it <m_particles.size(); it++)
    {
      // cout<<"it ="<<it<<endl;
      pos = m_particles[it].m_x[0];
      ipos[0] = floor(pos/m_dx);
      ipos[1] = 0;
      sk = (pos - ipos[0]*m_dx)/m_dx;
      Point e0 = getUnitv(0);
      Point ipoint(ipos);
      Point shift;
      // cout<<"ipoint[0] ="<<ipoint[0]<<", ipoint[1] ="<<ipoint[1]<<endl;
      for (int s = 0; s < m_W.supportSize(); s++)
	{
	  shift = e0;
	  shift *= s;
	  // cout<<"shift[0] ="<<shift[0]<<", shift[1] ="<<shift[1]<<endl;
	  // cout<<"s ="<< s<<endl;
	  // cout<<"interpcoeff ="<< interpcoeff<<",W(sk) = "<< interpcoeff*m_W.apply(sk, s)<<", W(1 - sk) = "<< m_W.apply(1-sk, s)<<endl;
	   m_particles[it].EField += a_Field[ipoint - shift]*m_W.apply(sk, s);
	   m_particles[it].EField += a_Field[ipoint + e0 + shift]*m_W.apply(1-sk, s);
	  // cout<<"ChargeL = "<< a_Charge[ipoint - shift]<<", ChargeR = "<<a_Charge[ipoint + e0 + shift]<<endl;
	}
    }
}
