#ifndef _INTERPOLATIONKERNEL_H_
#define _INTERPOLATIONKERNEL_H_

#include <memory>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
using namespace std;

class InterpolationKernel
{
public:
  /// Default Constructor
  InterpolationKernel();
  /// Destructor
  ~InterpolationKernel();
  /// Constructs polynomial interpolant from desired order and smoothness
  InterpolationKernel(int a_Order, int a_Smoothness = 0);
  /// Computes the value of the polynomial interpolant given x
  double apply(double a_x, int supportRegion);
  /// returns m_supportSize()
  int supportSize();
  
private:
  //supportSize is the number of points in the support minus 1.
  int m_supportSize;
  vector<vector<double> > m_coefficients;
  vector<vector<int> > m_powers;
  
};
#endif
