#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include "PowerItoI.H"
#include "FFT1D.H"
#include "FFT1DW.H"
#include <fftw3.h>
FFT1DW::FFT1DW(const unsigned int& a_M):
  FFT1D(a_M)
{
    in_fwd = (complex<double>*)fftw_malloc(sizeof(complex<double>)*m_N);
    out_fwd = (complex<double>*)fftw_malloc(sizeof(complex<double>)*m_N);
    in_bck = (complex<double>*)fftw_malloc(sizeof(complex<double>)*m_N);
    out_bck = (complex<double>*)fftw_malloc(sizeof(complex<double>)*m_N);

    forward = fftw_plan_dft_1d(m_N, reinterpret_cast<fftw_complex*>(in_fwd),
            reinterpret_cast<fftw_complex*>(out_fwd),
            FFTW_FORWARD, FFTW_ESTIMATE);
    backward = fftw_plan_dft_1d(m_N, reinterpret_cast<fftw_complex*>(in_bck),
            reinterpret_cast<fftw_complex*>(out_bck),
            FFTW_BACKWARD, FFTW_ESTIMATE);
}
FFT1DW::~FFT1DW(){
    fftw_destroy_plan(forward);
    fftw_destroy_plan(backward);
    fftw_free(in_fwd);
    fftw_free(out_fwd);
    fftw_free(in_bck);
    fftw_free(out_bck);
    fftw_cleanup();
}

void FFT1DW::forwardFFTCC(vector<complex<double> > & a_fHat,
                          const vector<complex<double> >& a_f) const
{
  //fftw_plan p;
  //This is super sloppy but there isn't a version off fftw_plan that takes a type const as input.
  //vector<complex<double> > in = a_f;
  //p = fftw_plan_dft_1d(m_N, reinterpret_cast<fftw_complex*>(&in[0]), reinterpret_cast<fftw_complex*>(&a_fHat[0]), FFTW_FORWARD, FFTW_ESTIMATE);
  //fftw_execute(p);
  //fftw_destroy_plan(p);
    copy(a_f.begin(), a_f.end(), in_fwd);
    fftw_execute(forward);
    copy(out_fwd, out_fwd+m_N, a_fHat.begin());
};
void FFT1DW::inverseFFTCC(vector<complex<double> > & a_f,
                          const vector<complex<double> > & a_fHat) const
{
  //fftw_plan p;
  //vector<complex<double> > in = a_fHat;
  //p = fftw_plan_dft_1d(m_N, reinterpret_cast<fftw_complex*>(&in[0]), reinterpret_cast<fftw_complex*>(&a_f[0]), FFTW_BACKWARD, FFTW_ESTIMATE);
  //fftw_execute(p);
  //fftw_destroy_plan(p);
    copy(a_fHat.begin(), a_fHat.end(), in_bck);
    fftw_execute(backward);
    copy(out_bck, out_bck+m_N, a_f.begin());
};
