#include <cstdio>
#include <iostream>
#include "WriteRectMDArray.H"
#include "DBox.H"
#include "RectMDArray.H"
#include "VisitWriter.H"
//#include "Stack.H"
//#include "ParticleSet.H"
using namespace std;


const char* MDWrite(RectMDArray<double,1>* a_array)
{
  static int fileCount = 0;
  static char nameBuffer[10];

  if (a_array == NULL)
    {
      return nameBuffer;
    }

  sprintf(nameBuffer, "md%06d",fileCount);
  MDWrite(nameBuffer, a_array);

  fileCount++;

  return nameBuffer;
};

void MDWrite(const char*           a_filename,
             RectMDArray<double>* a_array)
{
  FILE* fp = vtk_open_file(a_filename);
  double origin[3]={0,0,0};
  double dx=1.0;
  char* vars[1];
  vars[0] = "0";
  MDWrite(fp, *a_array, vars, origin, dx);
  }
const char* PWrite(const ParticleSet* a_array)
{
  static int fileCount = 0;
  static char nameBuffer[10];
  if(a_array == NULL)
    {
      return nameBuffer;
    }
  sprintf(nameBuffer, "PART.%d",fileCount);
  PWrite(nameBuffer, a_array);
  fileCount++;
  return nameBuffer;
}

void PWrite(const char* a_filename, const ParticleSet* a_p)
{
  if(a_filename == NULL || a_p == NULL)
    {
      return;
    }
  unsigned int size = a_p->m_particles.size();
  std::vector<double> x(3*size);
  for(unsigned int i=0; i<size; i++)
    {
      const Particle& p = a_p->m_particles[i];
      x[i*3] = p.m_x[0];
      x[i*3+1] = p.m_x[1];
#if DIM==3
      x[i*3+2] = p.m_x[2];
#else
      x[i*3+2] = 0.0;
#endif
    }

  write_point_mesh(a_filename, 0, size,
		   &(x[0]), 0, 0,
		   0, 0);
}
