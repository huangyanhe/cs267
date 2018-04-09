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
  ~PoissonSolver(){};
protected:
  double m_h;
  int m_M,m_N;
  FFTMD m_fftmd;
};
#endif
