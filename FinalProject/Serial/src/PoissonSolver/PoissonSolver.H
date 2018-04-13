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
  PoissonSolver(double a_h, int a_M, DBox a_box);
  void define(double a_h,int a_M,  DBox a_box);
  //Called by constructor for building eigenvalues for spectral/fd solver 
  void buildEigenvalues();
  //Used to Solve Poisson's equation on a periodic domain
  void Solve( RectMDArray<double>& a_rhs);
  ~PoissonSolver(){};
protected:
  double m_h;
  int m_M, m_N;
  FFTMD m_fftmd;
  DBox m_box;
  RectMDArray<double > m_eigenvalues;
};
#endif
