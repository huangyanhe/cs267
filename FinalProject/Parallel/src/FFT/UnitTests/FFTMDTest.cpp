#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include <memory>
#include <sys/time.h>
#include "PowerItoI.H"
#include "FFT1D.H"
#include "FFT1DW.H" 
#include "RectMDArray.H"
#include "DBox.H"
#include "FFTMD.H"

using namespace std;
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }

    gettimeofday( &end, NULL );

    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}
double test(const FFTMD& a_fftmd,double& a_seconds)
{
  int N = a_fftmd.getN();
  int low[DIM],high[DIM];

  for (int dir = 0;dir <DIM;dir++)
    {
      low[dir] = 0;
      high[dir] = N-1;
    }
  DBox b(low,high);
     
  RectMDArray<complex<double> > f(b);
  RectMDArray<complex<double> > fSave(b);
  RectMDArray<complex<double> > fHat(b);
  double h = 1./N;

  for (Point pt = b.getLowCorner();b.notDone(pt);b.increment(pt))
    {
      double x = 1.;
      for (int dir = 0;dir < DIM;dir++)
        {
         
          double y = (pt[dir]*h-.5);
          x += y*y*32*32;
        }
      f[pt] = complex<double >(exp(-x),0);
      fSave[pt] = f[pt];
    }
  double maxramp = 0.;
  double maxiamp = 0.;
  double maxamp = 0.;
  Point ptmax = getZeros();
  double time = read_timer();
  a_fftmd.forwardCC(f);
  a_seconds = a_seconds + read_timer() - time;
  for (Point pt = b.getLowCorner();b.notDone(pt);b.increment(pt))
    {
      if (real(f[pt])*real(f[pt]) + imag(f[pt])*imag(f[pt])  > maxamp )
        {
          maxramp = real(f[pt]);
          maxiamp = imag(f[pt]);
          maxramp = real(f[pt])*real(f[pt]) + imag(f[pt])*imag(f[pt]);         
          ptmax = pt;
        }
    }
  time = read_timer();
  a_fftmd.inverseCC(f);
  a_seconds = a_seconds + read_timer() - time;
  double maxerr = 0.;
  double minerr = 10000.;
    int normalize = Power(N,DIM);
  
  for (Point pt = b.getLowCorner();b.notDone(pt);b.increment(pt))
      {
        if (fabs(real(f[pt])/normalize) > maxerr )
          {
            maxerr = fabs(real(f[pt])/normalize - real(fSave[pt]));
          }
      }
  return maxerr;
}
int main(int argc, char* argv[])
{
  int M;
  double time;
  cout << "input log_2(N), N = number of points" << endl;
  cin >> M ;
  
  shared_ptr<FFT1D> p_fft;
  

  p_fft = dynamic_pointer_cast<FFT1D >(shared_ptr<FFT1DW>(new FFT1DW(M)));
   
  FFTMD fftmd(p_fft);
  double error = test(fftmd,time);
  cout << "test 1: error in Gaussian  = " << error << endl;
  cout << "time in FFTMD = " << time << " seconds" << endl;
  
};
