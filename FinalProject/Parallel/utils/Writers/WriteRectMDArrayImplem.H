#include <cstdio>
#include <iostream>


#include "VisitWriter.H"
//#include "Stack.H"

template<class T, unsigned int C>
void MDWrite(const string&  a_filename,
             RectMDArray<T, C>& a_array,
             string a_varnames[C],
             const array<double,DIM>&  a_corner,
             const double& a_h)
{
  /* if (a_filename == NULL)
    {
      return;
      }*/
  char* filename = new char[a_filename.length()+1];
  strcpy(filename,a_filename.c_str());
  char* varnames[C];
  for (int icomp = 0; icomp < C ; icomp++)
    {
      varnames[icomp] = new char[a_varnames[icomp].length()+1];
      strcpy(varnames[icomp],a_varnames[icomp].c_str());
    }

  FILE* fp = vtk_open_file(filename);
  MDWrite(fp, a_array, varnames, &a_corner[0], a_h);
  vtk_close_file(fp);
  for (int icomp = 0; icomp <C;icomp++)
    {
      delete [] varnames[icomp];
    }
  delete [] filename; 
}


template<class T, unsigned int C>
void MDWrite(FILE*                      a_fp,
             RectMDArray<T,C>&      a_array,
             const char* const          a_varnames[C],
             const double*              a_corner,
             const double&              a_h)
{


  int dim[3] = {1,1,1};

  int vardims[C];
  int centering[C];
  double* vars[C];
  double* x[3] {nullptr,nullptr,nullptr};

  Point lo, hi;
  const DBox& box = a_array.getDBox();
  lo = box.getLowCorner();
  hi = box.getHighCorner();
  for (int i = 0; i < DIM; i++)
    {
      dim[i] = hi[i] - lo[i] + 2;
      x[i] = new double[dim[i]+1];
      //x[i] = (double*)StackAlloc::alloc((dim[i]+1)*sizeof(double));
      for (int k = 0; k < dim[i]+1; k++)
      {
        x[i][k] = (double)(a_corner[i] + a_h * k);
      }
    }

  if (DIM == 2)
  {
    //x[2] = (double*)StackAlloc::alloc(2*sizeof(double));
    x[2] = new double[2];
    x[2][0] = 0.0;
    x[2][1] = 1.0;
  }

  for (int i = 0; i < C; i++)
  {
    vardims  [i] = 1;
    centering[i] = 0;

    vars[i] = &((a_array)(a_array.getDBox().getLowCorner(),i));
  }

  write_rectilinear_mesh(a_fp,
                         1,
                         dim,
                         x[0],x[1],x[2],
                         C,
                         vardims,
                         centering,
                         a_varnames,
                         vars);

  //StackAlloc::free(x[2]); // it's stack, have to do it backwards :-)
  //StackAlloc::free(x[1]);
  // StackAlloc::free(x[0]);
  delete [] x[2]; // it's stack, have to do it backwards :-)
  delete [] x[1];
  delete [] x[0];
}


