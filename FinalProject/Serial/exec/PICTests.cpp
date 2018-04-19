#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "RectMDArray.H"
#include "ParticleSet.H"
#include "ParticleVelocities.H"
#include "WriteRectMDArray.H" 
#include "VisitWriter.H"
#include "RK4.H"
auto removeParticle = [](Particle p) -> bool
{
  double minStrength = pow(10.0, -9);
  bool removeIfFalse =  (p.strength < minStrength);
  return removeIfFalse;
};
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
  cin >> test;
  cout << "input log_2(number of grid points)" << endl; 
  cin >> M;
  cout << "input order of interpolating function(lowest = 2)" << endl;
  int order;
  cin >> order;
  cout << "input smoothness of interpolating function(lowest = 0)" << endl;
  int smoothness;
  cin >> smoothness;  
  //cout << "input particle refinement factor" << endl;
  //unsigned int cfactor;
  //cin >> cfactor;
  cout << "enter stopping time" << endl;
  double timeStop;
  cin >> timeStop;

  //N = 2*Grid_Size in Matlab
  N = Power(2,M);
  double alpha = 0.01;
  double vmax = 10.0;
  double L = 4*M_PI;
  int cfactor = 2;
  double h = 1./N;
  double hp = h/cfactor; //pow(h,4./3.);
  int Np = 1./hp;
  hp = 1./Np;
  //cout<<"hp = "<<hp<<endl;
  double delta = h;
  int pcfactor = 4/cfactor;
  //if (pcfactor < 1 ) pcfactor = 1;
  //cout << "number of particles per cell = " << h*h/hp/hp << endl;

  int intlowCorner[DIM];
  int inthighCorner[DIM];
  for (int j = 0; j<DIM; j++)
    {
      intlowCorner[j] = 0;
      inthighCorner[j] = N-1;
    }
  array<double, DIM> lowCorn;
  for (int j =0; j<DIM; j++)
    {
      lowCorn[j] = 0.0;
    }
  
  Point lowCorner(intlowCorner);
  Point highCorner(inthighCorner);
  DBox domain(lowCorner, highCorner);
  domain.print();
  ParticleSet p(domain, h, lowCorn, M, L, order, smoothness );

  //Assumes equally spaced grids in x and v.
  //Changing this is pretty easy.
  int intPhaseXlowCorner[DIM];
  int intPhaseXhighCorner[DIM];
  int intPhaseVlowCorner[DIM];
  int intPhaseVhighCorner[DIM];
  for (int j = 0; j<DIM; j++)
    {
      intPhaseXlowCorner[j] = 0;
      intPhaseXhighCorner[j] = Np-1;
      intPhaseVlowCorner[j] = -Np;
      intPhaseVhighCorner[j] = Np-1;
    }
  Point PhaseXlowCorner(intPhaseXlowCorner);
  Point PhaseXhighCorner(intPhaseXhighCorner);
  Point PhaseVlowCorner(intPhaseVlowCorner);
  Point PhaseVhighCorner(intPhaseVhighCorner);
  DBox PhaseXSpace(PhaseXlowCorner, PhaseXhighCorner);
  DBox PhaseVSpace(PhaseVlowCorner, PhaseVhighCorner);
  //PhaseXSpace.print();
  //PhaseVSpace.print();
  //Still need to sort out what values this should be.
  double Modes[DIM];
  for (int j=0; j<DIM; j++)
    {
      Modes[j] = 1.0/2.0;
    }
  
  cout<<"Num Particles = "<<Power(Np, DIM)*Power(2*Np, DIM) <<endl;
  if (test == 1)
  {
    p.m_particles.resize(Power(Np, DIM)*Power(2*Np, DIM));
    int j = 0;
    for (Point ptV=PhaseVSpace.getLowCorner(); PhaseVSpace.notDone(ptV); PhaseVSpace.increment(ptV))
      {
	for (Point ptX=PhaseXSpace.getLowCorner(); PhaseXSpace.notDone(ptX); PhaseXSpace.increment(ptX), j++)
	  {
	    //p.m_particles.push_back();
	    for(int k=0; k<DIM; k++)
	      {
		//p.m_particles[j].m_x[k] = ptX[k]*hp*L;
		//p.m_particles[j].m_v[k] = ptV[k]*hp*vmax;
		p.m_particles[j].m_x[k] = ptX[k]*hp;
		p.m_particles[j].m_v[k] = ptV[k]*hp;
		p.m_particles[j].EField[k] = 0.0;
	      }
	    double fFirstTerm= 1.0;
	    double fSecondTerm =1.0;
	    for(int k=0; k<DIM; k++)
	      {
		//fFirstTerm *= exp(-p.m_particles[j].m_v[k]*p.m_particles[j].m_v[k]/2.0);
		//fSecondTerm *= exp(-p.m_particles[j].m_v[k]*p.m_particles[j].m_v[k]/2.0)*cos(Modes[k]*p.m_particles[j].m_x[k]);
		fFirstTerm *= exp(-p.m_particles[j].m_v[k]*vmax*p.m_particles[j].m_v[k]*vmax/2.0);
		fSecondTerm *= exp(-p.m_particles[j].m_v[k]*vmax*p.m_particles[j].m_v[k]*vmax/2.0)*cos(Modes[k]*p.m_particles[j].m_x[k]*L);
	      }
	    p.m_particles[j].strength = 1/sqrt(2.0*M_PI)*(fFirstTerm + alpha*fSecondTerm)*pow(hp*L,DIM)*pow(hp*vmax,DIM);
	  }
      }
  }

  p.m_particles.erase(remove_if(p.m_particles.begin(), p.m_particles.end(), removeParticle), p.m_particles.end());

  // for (auto it=p.m_particles.begin(); it!=p.m_particles.end(); ++it)
  //   {
  //     it->print();
  //   }
  
  //Need to check on some of these parameters.
  double dx = 1./N;
  cout << "number of particles = " << p.m_particles.size() << endl;
  ParticleShift kIn,kOut;
  kIn.init(p);
  kOut.init(p);
  kIn.setToZero();
  //cout<<"In Particle Velocities"<<endl;
  //ParticleVelocities pv(p); 
  //cout<<"Out Particle Velocities"<<endl;
  double time = 0.;
  double dt = 2.0/N;
  int m = 5000;

  RK4<ParticleSet,ParticleVelocities,ParticleShift> integrator;
  integrator.define(p);
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
