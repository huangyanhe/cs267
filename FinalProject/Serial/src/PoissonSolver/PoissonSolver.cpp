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
PoissonSolver::buildEigenvalues(RectMDArray<complex<double > >& a_cxarray,RectMDArray<double >& a_real)
{
  double scale[DIM];
  Point highCorner = m_box.getHighCorner();
  for (int j=0; j<DIM; j++)
    {
      //scale[j]= 2*M_PI/((highCorner[j]+1)*m_h[j]);
      scale[j]= 2*M_PI/(highCorner[j]*m_h[j]);
    }
  for (Point pt= m_box.getLowCorner(); m_box.notDone(pt); m_box.increment(pt))
    {
      double k = 0.0;
      for (int j=0; j<DIM; j++)
	{
	  if (2*pt[j] <= m_N[j] )
	    {
	      k += pow(scale*(pt[j]), 2);
	    }
	  else
	    {
	      k += pow(scale*(pt[j] - m_N[j]), 2);
	    }
	}
      complex<double> K (1/k, 0.0);
      m_eigenvalues[pt] = K;
    }
};
PoissonSolver::PoissonSolver():
  m_h{},
  m_N{},
  m_fftmd{},
  m_box{},
  m_eigenvalues{};
{
};
PoissonSolver::PoissonSolver(const double& a_h[DIM],int a_M[DIM], Dbox a_box)
{ 
  shared_ptr<FFT1DW> p_fftw1d = shared_ptr<FFT1DW>(new FFT1DW(a_M+1));
  shared_ptr<FFT1D> p_fft = dynamic_pointer_cast<FFT1D>(p_fftw1d);
  m_fftmd.define(p_fft);
  m_h = a_h;
  m_M = a_M;
  for (int j = 0; j<DIM; j++)
    {
      m_N[j] = Power(2,a_M[j]);
    }
  m_box = a_box;
  buildEigenvalues();
};
void PoissonSolver::define(const double& a_h[DIM],int a_M[DIM],  Dbox a_box[DIM])
{ 
  shared_ptr<FFT1DW> p_fftw1d = shared_ptr<FFT1DW>(new FFT1DW(a_M+1));
  shared_ptr<FFT1D> p_fft = dynamic_pointer_cast<FFT1D>(p_fftw1d);
  m_fftmd.define(p_fft);
  m_h = a_h;
  m_M = a_M;
  for (int j = 0; j<DIM; j++)
    {
      m_N[j] = Power(2,a_M[j]);
    }
  m_box = a_box;
  buildEigenvalues();
};
PoissonSolver::Solve( RectMDArray<double>& a_rhs)
{
  //DBox for grid
  //DBox rhsDomain = a_rhs.getDBox();
  Point low = m_box.getLowCorner();
  Point high = m_box.getHighCorner();

  assert(low == getZeros());
  assert(high == getOnes()*m_N);

  //low = high*(-1);
  //  DBox ddomain(low,high);
  RectMDArray<complex<double> > rhsDouble(m_box);
  complex<double> zero(0.,0.);
  rhsDouble.setVal(zero);
  //double scale = 1./pow(m_N*1.,DIM*2)/4;
  //Copy the double RMDA into a complex RMDA
  for (Point pt = m_box.getLowCorner(); m_box.notDone(pt); m_box.increment(pt))
    {
      rhsDouble[pt].real(a_rhs[pt]);
    }
  //RectMDArray<complex<double> > kernel(ddomain);
  RectMDArray<double > realOut(m_box);
  //m_kerPtr->getKernel(kernel,m_h);
  //need to determine if this is the correct transform to be using.
  m_fftmd.forwardCC(rhsDouble);
  for (Point pt = m_box.getLowCorner(); m_box.notDone(pt); m_box.increment(pt))
    {
      rhsDouble[pt] *= m_eigenvalues[pt];
    }
  m_fftmd.inverseCC(rhsDouble);
  a_rhs.setVal(0.);
  DBox bx(rhsDomain.getLowCorner(),rhsDomain.getHighCorner() - getOnes());
  for (Point pt = bx.getLowCorner();bx.notDone(pt);bx.increment(pt))
    {
      a_rhs[pt] = real(rhsDouble[pt]);
    }
}
