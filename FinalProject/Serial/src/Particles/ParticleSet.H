#ifndef _PARTICLESET_H_
#define _PARTICLESET_H_

#include <vector>
#include "Particle.H"
#include "RectMDArray.H"
#include "InterpolationKernel.H"

using namespace std;

class ParticleSet;
/// Represents an increment in a ParticleSet conforming to the dX class in RK4. 
class ParticleShift
{
public:
  /// vector of increments of particles.
  vector<DX> m_particles;
  ParticleShift(){};
  ~ParticleShift(){};
  /// Sets length of m_particles, initializes to zero via default constructor of DX.
  void init(const ParticleSet& a_particles);
  void initDelta(const ParticleSet& a_particles);
  /// m_particles[k] += a_rhs.m_particles[k]*a_scale.
  void increment(
                 double a_scale, 
                 const ParticleShift& a_rhs);
  void incrementVelocity(
			 double a_scale,
			 const ParticleShift& a_rhs);
  void incrementPositionandVelocity(
				    double a_scaleX,
				    double a_scaleV,
				    const ParticleShift& a_rhs);
  /// m_particles[k] *= a_scale
  void operator*=(double a_scale);
  /// reinitializes the values m_particles[k] to zero. Not used in RK4.
  void setToZero();
  void zeroEField();
};
/// Represents a collection of particles, conforming to the X class in RK4. Also contains the persistent data required to compute the right-hand side in ParticlesVelocities.
class ParticleSet
{
public:
  /// constructor. Called by main.
  /**
     a_kerptr: required to construct Hockney.
     a_box: DBox over which convolution is performed.
     a_dx: mesh spacing of PIC grid.
     a_lowCorner: DIM-tuple giving the location in physical space of the lowCorner of the PIC grid.
     a_M: N = 2^a_M is the number of grid points in the PIC grid.
   */
  // ParticleSet(
  //             DBox& a_box,
  //             double& a_dx, 
  //             array<double, DIM>& a_lowCorner,
  //             int a_M, int order=2, int smoothness=0);
  ParticleSet(
              DBox& a_box,
              double& a_dx, 
              array<double, DIM>& a_lowCorner,
              int a_M, double a_L, int order=2, int smoothness=0);
  ParticleSet(){};
  ~ParticleSet(){};

  //Data
  /// Main container for particle data.
  vector<Particle> m_particles;
  /// Mesh spacing.
  double m_dx;
  double m_L;
  /// Spatial grid.
  DBox m_box;
  /// Location in physical space of m_box.m_lowCorner;
  array<double,DIM> m_lowCorner;
  /// Interpolation Kernel used for deposition and interpolation
  InterpolationKernel m_W;
  /// RMDA used as temporary for deposition and interpolation. It is defined on Dbox.grow(m_W.support)
  //  RectMDArray<double, DIM> m_tempArray; 

  //Functions
  /// inplements increment function required for RK$.
  void incrementPositionandVelocity(const ParticleShift& a_shift);
  void increment(const ParticleShift& a_shift);
  void incrementDelta(const ParticleShift& a_shift, double a_dt);
  /// Bin the particles by spatial location 
  //void bin();
  /// Deposit the particles on the spatial grid
  void deposit(RectMDArray<double>& a_Charge, vector<Particle>& t_particles);
  void deposit(RectMDArray<double>& a_Charge);
  /// Interpolates the force from the spatial grid onto the particles
  void InterpolateForce(RectMDArray<double, DIM>& a_Field, vector<Particle>& t_particles);
  void InterpolateForce(RectMDArray<double, DIM>& a_Field, vector<DX>& t_particles);
  void InterpolateForce(RectMDArray<double, DIM>& a_Field);
  /// Wraps Particles that have (left) domain/crossed periodic boundary
  void wrapParticles(vector<Particle>& t_particles);
  void wrapParticles();
  
  
};
#endif
