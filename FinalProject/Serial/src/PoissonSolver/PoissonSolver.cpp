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
void PoissonSolver::buildEigenvalues()
{
  cout<<"m_N = "<<m_N<<endl;
  double scale[DIM];
  Point highCorner = m_box.getHighCorner();
  for (int j=0; j<DIM; j++)
    {
      //scale[j]= 2*M_PI/((highCorner[j]+1)*m_h);
      scale[j]= 2*M_PI/(m_N*m_h);
    }
  
  for (Point pt= m_box.getLowCorner(); m_box.notDone(pt); m_box.increment(pt))
    {
      cout<<"Point: ";
      pt.print();
      double k = 0.0;
      for (int j=0; j<DIM; j++)
	{
	  if (2*pt[j] <= m_N )
	    {
	      k += pow(scale[j]*(pt[j]), 2);
	    }
	  else
	    {
	      k += pow(scale[j]*(pt[j] - m_N), 2);
	    }
	}

      double K;
      if (k <= pow(10.0, -14.0))
	{
	  cout<<"k = 0"<<endl;
	  K =0.0;
	}
      else
	{
	  cout<<"in other case"<<endl;
	  K = 1/k;
	  cout<<real(K)<<endl;
	}
      m_eigenvalues[pt] = K;
    }
};
PoissonSolver::PoissonSolver():
  m_h{},
  m_N{},
  m_fftmd{},
  m_box{},
  m_eigenvalues{}
{};
PoissonSolver::PoissonSolver(double a_h, int a_M, DBox a_box)
{
  //This line was used because it was on a grid twice the size I think in hockney
  //shared_ptr<FFT1DW> p_fftw1d = shared_ptr<FFT1DW>(new FFT1DW(a_M+1));
  shared_ptr<FFT1DW> p_fftw1d = shared_ptr<FFT1DW>(new FFT1DW(a_M));
  shared_ptr<FFT1D> p_fft = dynamic_pointer_cast<FFT1D>(p_fftw1d);
  m_fftmd.define(p_fft);
  // for (int j=0; j<DIM; j++)
  //   {
  //     m_h[j] = a_h[j];
  //     m_M[j] = a_M[j];
  //   }
  // for (int j = 0; j<DIM; j++)
  //   {
  //     m_N[j] = Power(2,a_M[j]);
  //   }
  m_h = a_h;
  m_M = a_M;
  m_N = Power(2, a_M);
  m_box = a_box;
  m_eigenvalues.define(m_box);
  buildEigenvalues();
};
void PoissonSolver::define(double a_h,int a_M,  DBox a_box)
{ 
  //shared_ptr<FFT1DW> p_fftw1d = shared_ptr<FFT1DW>(new FFT1DW(a_M+1));
  shared_ptr<FFT1DW> p_fftw1d = shared_ptr<FFT1DW>(new FFT1DW(a_M));
  shared_ptr<FFT1D> p_fft = dynamic_pointer_cast<FFT1D>(p_fftw1d);
  m_fftmd.define(p_fft);
  // for (int j=0; j<DIM; j++)
  //   {
  //     m_h[j] = a_h[j];
  //     m_M[j] = a_M[j];
  //   }
  // for (int j = 0; j<DIM; j++)
  //   {
  //     m_N[j] = Power(2,a_M[j]);
  //   }
  m_h = a_h;
  m_M = a_M;
  m_N = Power(2, a_M);
  m_box = a_box;
  m_eigenvalues.define(m_box);
  buildEigenvalues();
};
void PoissonSolver::Solve( RectMDArray<double>& a_rhs)
{
  //DBox for grid
  //DBox rhsDomain = a_rhs.getDBox();
  Point low = m_box.getLowCorner();
  Point high = m_box.getHighCorner();

  //assert(low == getZeros());
  //assert(high == getOnes()*m_N);

  //low = high*(-1);
  //  DBox ddomain(low,high);
  RectMDArray<complex<double> > rhsDouble(m_box);
  complex<double> zero(0.,0.);
  rhsDouble.setVal(zero);
  //See FFTW site but applying a forward and inverse give the original scaled by the size of the array.
  //Hence the scale factor is 1/array size.
  double scale = 1./pow(m_N*1.,DIM);
  //Copy the double RMDA into a complex RMDA
  cout<<"made it to writing rhs double"<<endl;
  for (Point pt = m_box.getLowCorner(); m_box.notDone(pt); m_box.increment(pt))
    {
      rhsDouble[pt].real(a_rhs[pt]);
      pt.print();
      cout<<"rhsdouble[^] = "<<rhsDouble[pt]<<endl;
    }
  //RectMDArray<complex<double> > kernel(ddomain);
  RectMDArray<double > realOut(m_box);
  //m_kerPtr->getKernel(kernel,m_h);
  //need to determine if this is the correct transform to be using.
  cout<<"made it to fft"<<endl;
  m_fftmd.forwardCC(rhsDouble);
  cout<<"made it to multiplying eigenvalues"<<endl;
  for (Point pt = m_box.getLowCorner(); m_box.notDone(pt); m_box.increment(pt))
    {
      cout<<"eval["<< pt[0]<<"]= "<<m_eigenvalues[pt]<<endl;
      rhsDouble[pt] *= m_eigenvalues[pt];
      pt.print();
      cout<<"rhsdouble[^] = "<<rhsDouble[pt]<<endl;
    }
  cout<<"made it to inverse FFT"<<endl;
  m_fftmd.inverseCC(rhsDouble);
  a_rhs.setVal(0.);
  // DBox bx(m_box.getLowCorner(),m_box.getHighCorner() - getOnes());
  for (Point pt = m_box.getLowCorner(); m_box.notDone(pt); m_box.increment(pt))
    {
      a_rhs[pt] = real(rhsDouble[pt])*scale;
    }
}
