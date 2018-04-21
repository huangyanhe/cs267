#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <memory>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <cstring>
#include <memory>
#include <xmmintrin.h>
#include "InterpolationKernel.H"
#include "VisitWriter.H"
int main(int argc, char* argv[])
{
  cout << "input order of interpolating function" << endl;
  int order;
  cin >> order;

  cout << "input smoothness of interpolating function" << endl;
  int smoothness;
  cin >> smoothness;

  InterpolationKernel W(order, smoothness);

  double dx = 1.0/25.0;
  double supportBoundary = 3;
  int num_points = 150;

  double x[num_points];
  double out[num_points];

  ofstream myfile;
  //Filenaming
  string filename = "plotInterpolationKernelW" + to_string(order) + "r" + to_string(smoothness);
  myfile.open(filename);
  
  for (int j = 0; j<num_points; j++)
    {
      x[j] = -supportBoundary + j*dx;
      if (abs(x[j])<1)
	{
	  out[j] = W.apply(x[j], 0);
	}
      else
	{
	  out[j] = W.apply(x[j], 1);
	}
     
      myfile<< x[j]<<" "<<out[j]<<endl;
    }

  
  
  myfile.close();
  
  
  //const char* a_filename = filename.c_str();
  //const char* varnames = filename.c_str();
  
  //write_regular_mesh( a_filename, 0, *{num_points, 0, 0}, {1}, {1}, {1}, {varnames}, {out});
}
