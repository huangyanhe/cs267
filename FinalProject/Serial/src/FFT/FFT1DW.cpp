#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include "PowerItoI.H"
#include "FFT1D.H"
#include "FFT1DW.H"
#include <fftw3.h>
FFT1DW::FFT1DW(const unsigned int& a_M):
  FFT1D(a_M)
{
};
FFT1DW::~FFT1DW(){}; 
void FFT1DW::forwardFFTCC(vector<complex<double> > & a_fHat, 
                          const vector<complex<double> >& a_f) const
{
  fftw_plan p;
  //This is super sloppy but there isn't a version off fftw_plan that takes a type const as input. 
  vector<complex<double> > in = a_f;
  p = fftw_plan_dft_1d(m_N, reinterpret_cast<fftw_complex*>(&(in[0])), reinterpret_cast<fftw_complex*>(&(a_fHat[0])), FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
};
void FFT1DW::inverseFFTCC(vector<complex<double> > & a_f, 
                          const vector<complex<double> > & a_fHat) const
{
  fftw_plan p;
  vector<complex<double> > in = a_fHat;
  p = fftw_plan_dft_1d(m_N, reinterpret_cast<fftw_complex*>(&(in[0])), reinterpret_cast<fftw_complex*>(&(a_f[0])), FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
};
