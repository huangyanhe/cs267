#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <memory>
#include "PowerItoI.H"
#include "RectMDArray.H"
#include "WriteRectMDArray.H"
#include "DBox.H"
#include "FFT1DW.H"
#include "FFT1D.H"
#include "FFTMD.H"
#include "PoissonSolver.H"

using namespace std;
void copyReal(RectMDArray<complex<double> >& a_cxarray,RectMDArray<double >& a_real)
{
  DBox d = a_real.getDBox()&a_cxarray.getDBox();
  for (Point pt=d.getLowCorner();d.notDone(pt);d.increment(pt))
    {
      a_real[pt] = real(a_cxarray[pt]);
    }
};
PoissonSolver::PoissonSolver():
  m_h{},
  m_M{},
  m_N{},
  m_fftmd{},
{
};
PoissonSolver::PoissonSolver(const double& a_h,int a_M)
{ 
  shared_ptr<FFT1DW> p_fftw1d = shared_ptr<FFT1DW>(new FFT1DW(a_M+1));
  shared_ptr<FFT1D> p_fft = dynamic_pointer_cast<FFT1D>(p_fftw1d);
  m_fftmd.define(p_fft);
  m_h = a_h;
  m_M = a_M;
  m_N = Power(2,a_M);
};
void PoissonSolver::define(const double& a_h,int a_M)
{ 
  shared_ptr<FFT1DW> p_fftw1d = shared_ptr<FFT1DW>(new FFT1DW(a_M+1));
  shared_ptr<FFT1D> p_fft = dynamic_pointer_cast<FFT1D>(p_fftw1d);
  m_fftmd.define(p_fft);
  m_h = a_h;
  m_M = a_M;
  m_N = Power(2,a_M);
};
void PoissonSolver::Solve( RectMDArray<double>& a_rhs)
{
  DBox rhsDomain = a_rhs.getDBox();
  Point low = rhsDomain.getLowCorner();
  Point high = rhsDomain.getHighCorner();

  assert(low == getZeros());
  assert(high == getOnes()*m_N);

  //low = high*(-1);
  DBox ddomain(low,high);
  RectMDArray<complex<double> > rhsDouble(ddomain);
  complex<double> zero(0.,0.);
  rhsDouble.setVal(zero);
  double scale = 1./pow(m_N*1.,DIM*2)/4;
  for (Point pt = rhsDomain.getLowCorner();rhsDomain.notDone(pt);rhsDomain.increment(pt))
    {
      rhsDouble[pt].real(a_rhs[pt]);
    }
  //RectMDArray<complex<double> > kernel(ddomain);
  RectMDArray<double > realOut(ddomain);
  //m_kerPtr->getKernel(kernel,m_h);
  //need to determine if this is the correct transform to be using.
  m_fftmd.forwardCC(rhsDouble);
  for (Point pt = ddomain.getLowCorner();ddomain.notDone(pt);ddomain.increment(pt))
    {
      rhsDouble[pt] *= eigenvalue[pt];
    }
  m_fftmd.inverseCCcen(rhsDouble);
  a_rhs.setVal(0.);
  DBox bx(rhsDomain.getLowCorner(),rhsDomain.getHighCorner() - getOnes());
  for (Point pt = bx.getLowCorner();bx.notDone(pt);bx.increment(pt))
    {
      a_rhs[pt] = real(rhsDouble[pt])*scale;
    }
}
