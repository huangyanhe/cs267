#ifndef _PARTICLE_H_
#define _PARTICLE_H_
#include <cstdio>
#include <cmath>
#include <cassert>
#include <array>
#include <iostream>
using namespace std;
/// Class representing the increment of a single point.
/**
ParticleShift is a vector of these, one for each particle.
*/
class DX
{
public:
  /// Default contructor sets values to zero.
  DX()
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]=0;
      }
  }
  /// Data members are a DIM-tuple. 
  array<double, DIM> m_x;
  array<double, DIM> m_v;
  array<double, DIM> EField;
  /// Increment is used to implement ParticleShift::increment. 
  inline void increment(double a_scale, const DX& a_rhs)
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]+=a_scale*a_rhs.m_x[i];
      }
  }
  inline void addVelocity(double a_scale)
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]+=a_scale*m_v[i];
      }
  }
  inline void incrementVelocity(double a_scale, const DX& a_rhs)
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_v[i]+=a_scale*a_rhs.m_v[i];
      }
  }
  inline void incrementPositionandVelocity(double a_scaleX, double a_scaleV,const DX& a_rhs)
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]+=a_scaleX*a_rhs.EField[i];
	m_v[i]+=a_scaleV*a_rhs.EField[i];
      }
  }
  /// Used to implement ParticleShift::scale.
  /**
     Does not have a data member corresponding to strength since the vorticity doesn't change.
   */ 
  inline void operator*=(double a_scale)
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	EField[i]*=a_scale;
      }
  }
  void print()
  {
    cout<< "x= ";
      for (int k=0; k<DIM; k++)
        {
      cout<<m_x[k] <<", ";
        }
      cout<<"v=";
      for (int k=0; k<DIM; k++)
        {
	  cout<<m_v[k] <<", ";
        }
      cout<<"E=";
      for (int k=0; k<DIM; k++)
        {
	  cout<<EField[k] <<", ";
        }
      cout<<" "<<std::endl;
  }
};




//////////////////////////////////////////////////////////////////////////////////////////////
/// Represents a single particle in a ParticleSet.
//////////////////////////////////////////////////////////////////////////////////////////////

class Particle
{
public:
  /** m_x stores the spatial location of the particle,
      m_v stores the velocity of the particle, 
      EField stores the Field for each particle,
      strength stores the charge for each particle.
   */
  array<double, DIM> m_x;
  array<double, DIM> m_v;
  array<double, DIM> EField;
  double strength;
  Particle()
  {
    for( int k=0; k<DIM; k++)
      {
	m_x[k] = 0.0;
	m_v[k] = 0.0;
	EField[k] = 0.0;
      }
  }
  /// used to implement ParticleSet::increment.
  inline void increment(const DX& a_shift)
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]+=a_shift.EField[i];
      }
  }
  inline void incrementDelta(const DX& a_shift)
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]+=a_shift.m_x[i];
	m_v[i]+=a_shift.m_v[i];
      }
  }
  inline void addVelocity(double a_scale)
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]+=a_scale*m_v[i];
      }
  }
  inline void incrementPositionandVelocity(const DX& a_shift)
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]+=a_shift.m_x[i];
	m_v[i]+=a_shift.m_v[i];
      }
  }
  void print()
  {
    cout<<"strength= "<< strength<<", ";
    cout<< "x= ";
      for (int k=0; k<DIM; k++)
        {
      cout<<m_x[k] <<", ";
        }
      cout<<"v=";
      for (int k=0; k<DIM; k++)
        {
	  cout<<m_v[k] <<", ";
        }
      cout<<"E=";
      for (int k=0; k<DIM; k++)
        {
	  cout<<EField[k] <<", ";
        }
      cout<<" "<<std::endl;
  }
};

#endif
