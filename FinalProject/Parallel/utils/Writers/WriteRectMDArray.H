#ifndef _WRITERECTMDARRAY_H_
#define _WRITERECTMDARRAY_H_
#include <cstring>
#include <cstdio>
#include <array>
#include "RectMDArray.H"
#include "ParticleSet.H"

using namespace std;

// first two are handy when debugging code.  Can call both from the debugger

const char* MDWrite(RectMDArray<double>* a_array);

void MDWrite(const char* a_filename,
             RectMDArray<double>* a_array);

template<class T, unsigned int C>
void MDWrite(const string& a_filename,
             RectMDArray<T, C>& a_array,
             string a_varnames[C],
             const array<double,DIM>& a_corner,
             const double& a_h);

template<class T, unsigned int C>
void MDWrite(FILE*                      a_fp,
             RectMDArray<T,C>&      a_array,
             const char* const          a_varnames[C],
             const double*              a_corner,
             const double&              a_h);

void PWrite(const char* filename, const ParticleSet* a_particles);

const char* PWrite(const ParticleSet* a_particles);

#include "WriteRectMDArrayImplem.H"
#endif 
