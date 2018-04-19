#include "ParticleSet.H"

//Constructors
// articleSet::ParticleSet(
//               DBox& a_box,
//               double& a_dx, 
//               array<double, DIM>& a_lowCorner,
//               int a_M, int order, int smoothness):
//   m_particles{},
//   m_dx{a_dx},
//   m_box{a_box},
//   m_lowCorner{a_lowCorner},
//   m_W{order, smoothness}
// {
//   //  m_tempArray.define(a_box.grow(m_W.supportSize()));
// }
ParticleSet::ParticleSet(
	      DBox& a_box,
              double& a_dx, 
              array<double, DIM>& a_lowCorner,
              int a_M, double a_L, int order, int smoothness):
  m_particles{},
  m_dx{a_dx},
  m_box{a_box},
  m_lowCorner{a_lowCorner},
  m_W{order, smoothness},
  m_L{a_L}
{
  //  m_tempArray.define(a_box.grow(m_W.supportSize()));
}
// Functions
void ParticleSet::increment(const ParticleShift& a_shift)
{
  for (int j = 0; j<a_shift.m_particles.size(); j++)
    {
      m_particles[j].increment(a_shift.m_particles[j]);
    }
}
void ParticleSet::incrementPositionandVelocity(const ParticleShift& a_shift)
{
  for (int j = 0; j<a_shift.m_particles.size(); j++)
    {
      m_particles[j].increment(a_shift.m_particles[j]);
    }
}
// Dimensionally Independent Deposition
void  ParticleSet::deposit(RectMDArray<double>& a_Charge, vector<Particle>& t_particles)
{
  array<double,DIM> pos;
  double interpcoeff;
  // for (auto it= t_particles.begin(); it!= t_particles.end(); ++it)
  //   {
  //     it->print();
  //   }
  // cout<<"m_dx = "<< m_dx<<endl;
  // m_box.print();
  // a_Charge.getDBox().print();
  for (int it = 0; it <t_particles.size(); it++)
    {
      array<int,DIM> iposLow, iposHigh;
      for (int l = 0; l < DIM; l++)
        {
	  pos[l] = t_particles[it].m_x[l];
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
      Point LC = (LowCorner - Shift);
      Point HC = (HighCorner+ Shift);
      DBox SupportBox(LC, HC);
      //SupportBox.print();
      //Only place that it seems like the scaling needs to be taken into account.
      interpcoeff = 1/pow(m_dx*m_L, DIM)*t_particles[it].strength;
      for (Point s = SupportBox.getLowCorner(); SupportBox.notDone(s); SupportBox.increment(s))
	{
	  double KernelProduct = 1.0;
	  for (int j =0; j<DIM; j++)
	    {
	      int region = min(abs(s[j] - LowCorner[j]), abs(s[j] - HighCorner[j]));
	      double val = abs(s[j]*m_dx - pos[j])/m_dx;
	      //cout<<"s[j] = "<<s[j]<<endl;
	      //cout<<"val = "<<val<<endl;
	      KernelProduct *= m_W.apply(val, region);
	    }
	  a_Charge[s] += interpcoeff*KernelProduct;
	}
    }
}
void  ParticleSet::deposit(RectMDArray<double>& a_Charge)
{
  deposit( a_Charge, m_particles);
}
// void  ParticleSet::deposit(RectMDArray<double>& a_Charge)
// {
//   array<double,DIM> pos;
//   double interpcoeff;
//   for (int it = 0; it <m_particles.size(); it++)
//     {
//       array<int,DIM> iposLow, iposHigh;
//       for (int l = 0; l < DIM; l++)
//         {
// 	  pos[l] = m_particles[it].m_x[l];
// 	  iposLow[l] = floor(pos[l]/m_dx);
// 	  iposHigh[l] = ceil(pos[l]/m_dx); 
// 	}

//       Point Shift = getUnitv(0);
//       Shift *= 0;
//       for (int j =0; j<DIM; j++)
// 	{
// 	  Shift += getUnitv(j);
// 	}
//       Shift *= (m_W.supportSize() - 1);
//       Point HighCorner(iposHigh);
//       Point LowCorner(iposLow);
//       Point LC = LowCorner - Shift;
//       Point HC = HighCorner+ Shift;
//       DBox SupportBox(LC, HC); 
//       interpcoeff = 1/m_dx*m_particles[it].strength;
//       for (Point s = SupportBox.getLowCorner(); SupportBox.notDone(s); SupportBox.increment(s))
// 	{
// 	  double KernelProduct = 1.0;
// 	  for (int j =0; j<DIM; j++)
// 	    {
// 	      int region = min(abs(s[j] - LowCorner[j]), abs(s[j] - HighCorner[j]));
// 	      double val = abs(s[j]*m_dx - pos[j])/m_dx; 
// 	      KernelProduct *= m_W.apply(val, region);
// 	    }
// 	  a_Charge[s] += interpcoeff*KernelProduct;
// 	}
//     }
// }

//Dimensionally Independent Interpolation
void  ParticleSet::InterpolateForce(RectMDArray<double, DIM>& a_Field, vector<Particle>& t_particles)
{
  for (int it = 0; it <t_particles.size(); it++)
    {
      array<int,DIM> iposLow, iposHigh;
      array<double,DIM> pos;
      for (int l = 0; l < DIM; l++)
        {
	  pos[l] = t_particles[it].m_x[l];
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
	      t_particles[it].EField[j] += a_Field(s, j)*KernelProduct;
	    }
	}
    }
}
void  ParticleSet::InterpolateForce(RectMDArray<double, DIM>& a_Field, vector<DX>& t_particles)
{
  for (int it = 0; it <t_particles.size(); it++)
    {
      array<int,DIM> iposLow, iposHigh;
      array<double,DIM> pos;
      for (int l = 0; l < DIM; l++)
        {
	  pos[l] = t_particles[it].m_x[l];
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
	      t_particles[it].EField[j] += a_Field(s, j)*KernelProduct;
	    }
	}
    }
}
void  ParticleSet::InterpolateForce(RectMDArray<double, DIM>& a_Field)
{
  InterpolateForce(a_Field, m_particles);
}
// void  ParticleSet::InterpolateForce(RectMDArray<double>& a_Field)
// {
  
//   for (int it = 0; it <m_particles.size(); it++)
//     {
//       array<int,DIM> iposLow, iposHigh;
//       array<double,DIM> pos;
//       for (int l = 0; l < DIM; l++)
//         {
// 	  pos[l] = m_particles[it].m_x[l];
// 	  iposLow[l] = floor(pos[l]/m_dx);
// 	  iposHigh[l] = ceil(pos[l]/m_dx); 
// 	}

//       Point Shift = getUnitv(0);
//       Shift *= 0;
//       for (int j =0; j<DIM; j++)
// 	{
// 	  Shift += getUnitv(j);
// 	}
//       Shift *= (m_W.supportSize() - 1);
//       Point HighCorner(iposHigh);
//       Point LowCorner(iposLow);
//       Point LC = LowCorner - Shift;
//       Point HC = HighCorner+ Shift;
//       DBox SupportBox(LC, HC); 
//       for (Point s = SupportBox.getLowCorner(); SupportBox.notDone(s); SupportBox.increment(s))
// 	{
// 	  double KernelProduct = 1.0;
// 	  for (int j =0; j<DIM; j++)
// 	    {
// 	      int region = min(abs(s[j] - LowCorner[j]), abs(s[j] - HighCorner[j]));
// 	      double val = abs(s[j]*m_dx - pos[j])/m_dx; 
// 	      KernelProduct *= m_W.apply(val, region);
// 	    }
// 	  for (int j = 0; j<DIM; j++)
// 	    {
// 	      m_particles[it].EField[j] += a_Field(s, j)*KernelProduct;
// 	    }
// 	}
//     }
// }

//Wraps Particles that have crossed boundary.
void ParticleSet::wrapParticles(vector<Particle>& t_particles)
{
  double highBoundary[DIM];
  double lowBoundary[DIM];
  double domainSpecs[DIM];
  for (int j=0; j<DIM; j++)
    {  
      highBoundary[j] = (m_box.getHighCorner()[j]+1)*m_dx;
      lowBoundary[j] = (m_box.getLowCorner()[j])*m_dx;
      domainSpecs[j] = highBoundary[j] - lowBoundary[j];
    }
  //The ordering of this should be explored. Since it is a vector of particles this should be best.
  for (int it = 0; it < t_particles.size(); it++)
    {
      for (int j =0; j<DIM; j++)
	{
	  //Determines if the marticles left the spatial domain.
	  //Might be good to implement this for different 
	  if (t_particles[it].m_x[j]>highBoundary[j])
	    {
	      t_particles[it].m_x[j] = t_particles[it].m_x[j] - domainSpecs[j];
	    }
	  else if (t_particles[it].m_x[j]<lowBoundary[j])
	    {
	      t_particles[it].m_x[j] = t_particles[it].m_x[j] + domainSpecs[j];
	    }
	}
    }
}
void ParticleSet::wrapParticles()
{
  wrapParticles(m_particles);
}
// void ParticleSet::wrapParticles()
// {
//   double highBoundary[DIM];
//   double lowBoundary[DIM];
//   double domainSpecs[DIM];
//   for (int j=0; j<DIM; j++)
//     {  
//       highBoundary[j] = (m_box.getHighCorner()[j]+1)*m_dx;
//       lowBoundary[j] = (m_box.getLowCorner()[j])*m_dx;
//       domainSpecs[j] = highBoundary[j] - lowBoundary[j];
//     }
//   //The ordering of this should be explored. Since it is a vector of particles this should be best.
//   for (int it = 0; it < m_particles.size(); it++)
//     {
//       for (int j =0; j<DIM; j++)
// 	{
// 	  //Determines if the marticles left the spatial domain.
// 	  //Might be good to implement this for different 
// 	  if (m_particles[it].m_x[j]>highBoundary[j])
// 	    {
// 	      m_particles[it].m_x[j] = m_particles[it].m_x[j] - domainSpecs[j];
// 	    }
// 	  else if (m_particles[it].m_x[j]<lowBoundary[j])
// 	    {
// 	      m_particles[it].m_x[j] = m_particles[it].m_x[j] + domainSpecs[j];
// 	    }
// 	}
//     }
// }
