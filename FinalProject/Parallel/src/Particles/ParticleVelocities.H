#ifndef _PARTICLEVELOCITIES_H_
#define _PARTICLEVELOCITIES_H_
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "DBox.H"
#include "RectMDArray.H"
#include "ParticleSet.H"
#include "PoissonSolver.H"
#include "mpi.h"
using namespace std;
/// Class for computing the RHS in RK4, corresponding to the template parameter F.
class ParticleVelocities
{
public:
  /// Member function that computes increment to solution.Conforms to the m_f(...) interface in RK4. 
  /** Implements an operator that evaluates a_dt*F(t,a_X+a_k), given the inputs (a_k, a_time, a_dt, a_X ). It first computes a temporary a_X + a_k, then evaluates the right-hand side, scales the result by  a_dt and stores it in  a_k. 
   */
  /// Maximum support size (assumed to be set by the interpolation kernel)
  int m_supportSize;
  /// DomainSize (= highCorner+1). 0 is associated with N+1th point.
  int m_domainSize[DIM];
  /// Mesh spacing.
  double m_dx;
  double m_L;
  /// Spatial grid (unenlarged).
  DBox m_box;
  /// Boundary intersection with ghost cells
  DBox m_bdry[2*DIM];
  PoissonSolver PS;
  ParticleVelocities();
  ParticleVelocities(ParticleSet& a_state);
  void define(ParticleSet& a_state);
  /// Function that deals with Ghost Cells for deposition, fd, and interpolation 
  void getGhostDeposition(RectMDArray<double>& enlargedGrid);
  //void setGhost(RectMDArray<double>& enlargedGrid);
  void setGhost(RectMDArray<double>& enlargedGrid);
  void setGhostMD(RectMDArray<double, DIM>& enlargedGrid);
  /// Main operator called by RK4
  void operator()(ParticleShift& a_k, 
                     const double& a_time, const double& dt, 
                     ParticleSet& a_state);
};
#endif
