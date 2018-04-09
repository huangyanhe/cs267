#ifndef _POISSONSOLVER_
#define _POISSONSOLVER_
#include "RectMDArray.H"
#include "FFTMD.H"
#include <array>
#include <memory>
using namespace std;
class PoissonSolver
{
public:
  
  PoissonSolver();
  PoissonSolver(const double& a_h, int a_M);
  void PoissonSolver(const double& a_h, int a_M);
  void PoissonSolver(RectMDArray<double>& a_rhs);
  //Called by constructor for building eigenvalues for spectral/fd solver 
  void buildEigenvalues();
  void Solve( RectMDArray<double>& a_rhs)
  ~PoissonSolver(){};
protected:
  double m_h[DIM];
  int m_M[DIM], m_N[DIM];
  FFTMD m_fftmd;
  DBox m_box;
  RectMDArray<complex<double> > m_eigenvalues;
};
#endif
