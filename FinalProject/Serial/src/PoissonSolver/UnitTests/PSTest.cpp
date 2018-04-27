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
#include "RectMDArray.H"
#include "ParticleSet.H"
#include "DBox.H"
#include "Point.H"
#include "PoissonSolver.H"
#include "VisitWriter.H"
int main(int argc, char* argv[])
{  
  // // Grid Setup 
  // int intlowCorner[2] = {0, 0};
  // int inthighCorner[2] = {7, 7};
  // int inthighCornerX[2] = {4, 0};
  // Point lowCorner(intlowCorner);
  // Point highCorner(inthighCorner);
  // Point highCornerX(inthighCornerX);
  // DBox phaseGridBox(lowCorner, highCorner);
  // DBox xGridBox(lowCorner, highCornerX);
  // RectMDArray<double> Charge(xGridBox);
  // Charge.setVal(0.0);
  
  // // ParticleSet build
  // int M =3; //2^M on PIC grid
  // array<double,2> lowCorn;
  // lowCorn[0] = 0.0;
  // lowCorn[1] = 0.0;
  // double dx = 1.0/4.0;
  // ParticleSet TestSet( phaseGridBox, dx, lowCorn, M, order, smoothness);
  
  
  // Build grid
  int M =2; //2^M on PIC grid
  double N = pow(2.0, M);
  double h = 2*M_PI*1.0/N;

  int intlowCorner[DIM];
  int inthighCorner[DIM];
  for (int j = 0; j<DIM; j++)
    {
      intlowCorner[j] = 0;
      inthighCorner[j] = N-1;
    }
  
  Point lowCorner(intlowCorner);
  Point highCorner(inthighCorner);
  DBox GridBox(lowCorner, highCorner);


  
  PoissonSolver PS(h, M, 1.0, GridBox);
  cout<<"made it here"<<endl;
  RectMDArray<double > RHS(GridBox);
  
  for (Point p=GridBox.getLowCorner(); GridBox.notDone(p); GridBox.increment(p))
    {
      RHS[p] = 1.0;
      for (int j=0; j<DIM; j++)
	{
	  RHS[p] *= sin(p[j]*h);
	}
    }
  cout<<"Initialization"<<endl;
  for (Point p=GridBox.getLowCorner(); GridBox.notDone(p); GridBox.increment(p))
    {
      cout<< "x = ";
      for (int j =0; j<DIM; j++)
	{
	  cout<<p[j]<<", ";
	}
      cout<<" RHS = "<<RHS[p]<<endl;
    }
   
  PS.Solve(RHS);

  cout<<"Solution"<<endl;
  for (Point p=GridBox.getLowCorner(); GridBox.notDone(p); GridBox.increment(p))
    {
      cout<< "x = ";
      for (int j =0; j<DIM; j++)
	{
	  cout<<p[j]<<", ";
	}
      cout<<" RHS = "<<RHS[p]<<endl;
    }

  cout<<"Exact"<<endl;
  for (Point p=GridBox.getLowCorner(); GridBox.notDone(p); GridBox.increment(p))
    {
      cout<< "x = ";
      for (int j =0; j<DIM; j++)
	{
	  cout<<p[j]<<", ";
	}
      double exact = 1.0;
      for (int j =0; j<DIM; j++)
	{
	  exact *= -sin(2.0*M_PI*p[j]/N)/pow(2*M_PI/N,2)*DIM;
	}
      cout<<" exact = "<<exact<<endl;
    }
  
}
