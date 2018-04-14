#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "RectMDArray.H"
#include "ParticleSet.H"
#include "ParticleVelocities.H"
#include "CutoffKernel.H"
#include "WriteRectMDArray.H" 
#include "VisitWriter.H"
#include "RK4.H"
void outField(ParticleSet& p, int a_coarsenFactor)
{
  int coarsenFactor = a_coarsenFactor;
  DBox bx = p.m_box.coarsen(coarsenFactor);
  double h = p.m_dx*coarsenFactor;
  RectMDArray<double> outVort(bx);
  array<int,DIM> ipos;
  array<double,DIM> xpos;
  double weight;
  Point e0 = getUnitv(0);
  Point e1 = getUnitv(1);
  
  outVort.setVal(0.);
  for (int k = 0; k < p.m_particles.size(); k++)
    {
      for (int l = 0; l < DIM; l++)
        {
          double newpos = p.m_particles[k].m_x[l]; 
          ipos[l] = newpos/h;
          xpos[l] = (newpos - ipos[l]*h)/h;
        }
      Point pt(ipos);
      assert(p.m_box.contains(pt));
      for (int l0=0; l0 < DIM;l0++)
        {
          for (int l1=0;l1 < DIM ; l1++)
            {
              outVort[pt+e0*l0 + e1*l1] += 
                (1.-xpos[0] + (2*xpos[0] - 1.)*l0)*
                (1.- xpos[1] + (2*xpos[1] - 1.)*l1)*p.m_particles[k].strength/coarsenFactor/coarsenFactor;
            }
        }
    }
  const char* foo = MDWrite(&outVort);
};
int main(int argc, char* argv[])
{
  unsigned int M;
  unsigned int N;
  cout << "input test = 1 (Linear Landau Damping), 2, other" << endl;
  int test;
  cout << "input log_2(number of grid points)" << endl; 
  cin >> M;
  cin >> test;
  cout << "input order of interpolating function(lowest = 2)" << endl;
  int order;
  cin >> order;
  cout << "input smoothness of interpolating function(lowest = 0)" << endl;
  int smoothness;
  cin >> smoothness;  
  cout << "input particle refinement factor" << endl;
  unsigned int cfactor;
  cin >> cfactor;
  cout << "enter stopping time" << endl;
  double timeStop;
  cin >> timeStop;
  
  N = Power(2,M);
  double h = 1./N;
  double hp = h/cfactor; //pow(h,4./3.);
  int Np = 1./hp;
  hp = 1./Np;
  double delta = h;
  int pcfactor = 4/cfactor;
  if (pcfactor < 1 ) pcfactor = 1;
  cout << "number of particles per cell = " << h*h/hp/hp << endl;

  int intlowCorner[DIM];
  int inthighCorner[DIM];
  for (int j = 0; j<DIM; j++)
    {
      intlowCorner[j] = 0;
      inthighCorner[j] = Np-1;
    }
  array<double, DIM> lowCorn;
  for (int j =0; j<DIM; j++)
    {
      lowCorn[j] = 0.0;
    }
  
  Point lowCorner(intlowCorner);
  Point highCorner(inthighCorner);
  DBox domain(lowCorner, highCorner);
  
  ParticleSet p(domain, hp, lowCorn, M, order, smoothness );

  //Assumes equally spaced grids in x and v.
  //Changing this is pretty easy.
  int intPhaselowCorner[2*DIM];
  int intPhasehighCorner[2*DIM];
  for (int j = 0; j<2*DIM; j++)
    {
      intPhaselowCorner[j] = 0;
      inthighCorner[j] = N-1;
    }
  Point PhaselowCorner(intPhaselowCorner);
  Point PhasehighCorner(intPhasehighCorner);
  DBox PhaseSpace(PhaselowCorner, PhasehighCorner);
  //Still need to sort out what values this should be.
  double k[DIM];
  for (int j=0; j<DIM; j++)
    {
      k[j] = 1.0/2.0;
    }
  double alpha = 0.05;
  
  if (test == 1)
  {
    p.m_particles.resize(Power(N, 2));
    int j = 0;
    for (Point pt=PhaseSpace.getLowCorner(); PhaseSpace.notDone(pt); PhaseSpace.increment(pt), j++)
      {
	for(int k=0; k<DIM; k++)
	  {
	    p.m_particles[j].m_x[k] = pt[k]*h;
	    p.m_particles[j].m_v[k] = pt[DIM-1+k]*h;
	    p.m_particles[j].EField[k] = 0.0;
	  }
	double fFirstTerm, fSecondTerm;
	for(int k=0; k<DIM; k++)
	  {
	    fFirstTerm *= exp(-p.m_particles[j].m_v[k]*p.m_particles[j].m_v[k]/2.0);
	    fSecondTerm *= exp(-p.m_particles[j].m_v[k]*p.m_particles[j].m_v[k]/2.0)*cos(k[dim]*p.m_particles[j].m_x[k]);
	  }
	p.m_particles[j].strength = 1/(2.0*M_PI)*(fFirstTerm + alpha*fSecondTerm)*pow(h,2*DIM);
      }
    
  }
  //Need to check on some of these parameters.
  double dx = 1./N;
  cout << "number of particles = " << p.m_particles.size() << endl;
  ParticleShift kIn,kOut;
  kIn.init(p);
  kOut.init(p);
  kIn.setToZero();
  ParticleVelocities pv; 
  double time = 0.;
  double dt = 140*.025/N;
  int m = 5000;
  
  RK4<ParticleSet,ParticleVelocities,ParticleShift> integrator;
// #if ANIMATION
//   outField(p,pcfactor);
//   PWrite(&p);
// #endif 
  for(int i=0; i<m; i++)
    {
      integrator.advance(time, dt, p);
      time = time + dt;
      cout << "time = " << time << "  dt " << dt << endl;
// #if ANIMATION
//       outField(p,pcfactor);
//       PWrite(&p);
// #endif
      if (time >= timeStop) 
        {
          break;
        }
    }
  // if (!((test == 1) || (test == 2) || (test == 3)))
  //   {
  //     outField(p,pcfactor);
  //     PWrite(&p);
  //   }
}
