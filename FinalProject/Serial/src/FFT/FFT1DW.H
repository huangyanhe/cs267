#ifndef _FFT1DW_H_
#define _FFT1DW_H_
#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include "PowerItoI.H"
#include "FFT1D.H"
#include <fftw3.h>
//#include "fftw3.h"
class FFT1DW:public FFT1D
{
public:
  FFT1DW(const unsigned int& a_M);
  virtual ~FFT1DW();
  virtual void forwardFFTCC(vector<complex<double> > & a_fHat,
			    const vector<complex<double> >& a_f) const;
  virtual void inverseFFTCC(vector<complex<double> > & a_f,
                            const vector<complex<double> > & a_fHat) const;
private:
  //Need to think about if there should be any member data.
  complex<double>* in_fwd;
  complex<double>* out_fwd;
  complex<double>* in_bck;
  complex<double>* out_bck;
  fftw_plan forward;
  fftw_plan backward;
};
#endif
