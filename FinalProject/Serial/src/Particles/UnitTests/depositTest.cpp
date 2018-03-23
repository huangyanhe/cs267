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
#include "VisitWriter.H"
int main(int argc, char* argv[])
{
  // interpolation =2,
  cout << "Test (deposition = 1, deposit and then intepolate back with force = 2, etc.)" << endl;
  int test;
  cin >> test;
  cout << "input order of interpolating function(lowest = 2)" << endl;
  int order;
  cin >> order;
  cout << "input smoothness of interpolating function(lowest = 0)" << endl;
  int smoothness;
  cin >> smoothness;  
  
  // Grid Setup 
  int intlowCorner[2] = {0, 0};
  int inthighCorner[2] = {7, 7};
  int inthighCornerX[2] = {4, 0};
  Point lowCorner(intlowCorner);
  Point highCorner(inthighCorner);
  Point highCornerX(inthighCornerX);
  DBox phaseGridBox(lowCorner, highCorner);
  DBox xGridBox(lowCorner, highCornerX);
  RectMDArray<double> Charge(xGridBox);
  Charge.setVal(0.0);
  
  // ParticleSet build
  int M =3; //2^M on PIC grid
  array<double,2> lowCorn;
  lowCorn[0] = 0.0;
  lowCorn[1] = 0.0;
  double dx = 1.0/4.0;
  ParticleSet TestSet( phaseGridBox, dx, lowCorn, M, order, smoothness);

  
  int num_points = phaseGridBox.size(0);

  double x[num_points];
  double out[num_points];

  for (int j =0; j<num_points; j++)
    {
      Particle part;
      x[j] = j*dx/2.0 + 1.0/16.0;
      part.m_x[0] = j*dx/2.0 + 1.0/16.0;
      part.m_x[1] = 0.0;
      part.strength = 1.0;
      TestSet.m_particles.push_back(part);
      TestSet.m_particles[j].EField = 0.0;
    }

  for (int j =0; j<TestSet.m_particles.size(); j++)
    {
      cout<< "xparticle = "<< TestSet.m_particles[j].m_x[0]<< "particle charge = "<< TestSet.m_particles[j].strength<<endl;
    }
  
  TestSet.deposit(Charge);
  
  for (Point p=Charge.getDBox().getLowCorner(); Charge.getDBox().notDone(p); Charge.getDBox().increment(p))
     {
       cout<< "x = "<< p[0]<< ", Charge = "<<Charge[p]<<endl;
       out[p[0]] = Charge[p];
     }

  if (test == 1)
    {
      ofstream myfile;
      //Filenaming
      string filename = "plotDepositionW" + to_string(order) + "r" + to_string(smoothness);
      myfile.open(filename);
  
      for (int j = 0; j < num_points; j++)
	{
	  myfile<< x[j]<<" "<<out[j]<<endl;
	}
    }

  if (test == 2)
    {
      TestSet.InterpolateForce(Charge);

       ofstream myfile;
      //Filenaming
      string filename = "plotDepositionForceW" + to_string(order) + "r" + to_string(smoothness);
      myfile.open(filename);
      myfile << "Charge"<<endl;
      for (Point p=Charge.getDBox().getLowCorner(); Charge.getDBox().notDone(p); Charge.getDBox().increment(p))
     {
       myfile<<  p[0]<< " "<< Charge[p]<<endl;
     }


      myfile << "Force"<<endl;

      for (int j = 0; j < TestSet.m_particles.size(); j++)
	{
	  myfile<< TestSet.m_particles[j].m_x[0]<<" "<< TestSet.m_particles[j].EField<<endl;
	}
      
      
      
    }

  
}
