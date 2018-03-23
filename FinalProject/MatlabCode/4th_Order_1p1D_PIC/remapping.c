/*
 * deposition.c - speed up for depositions
 *
 *This script handles the nested loops for the deposition and remapping 
 *steps of the PIC algorithm.  
 *
 * The calling syntax is:
 *
 *		newWeights = arrayProduct(oldWeights, oldPositions, newPositions)
 *
 * This is a MEX file for MATLAB.
*/

#include "mex.h"
#include "math.h"

/* Remapping Function */
double W3(double y)
{
    double x = fabs(y);
    double z;
    if (x <= 1)
    {
        z = 1-2.5*pow(x,2)+1.5*pow(x,3);
    }
    else if (x <= 2)
    {
        z = 2-4*x+2.5*pow(x,2)-0.5*pow(x,3);
    }
    else
    {
        z = 0;
    }
    return z;
    
}

double W4(double y)
{
    double x = fabs(y);
    double z;
    if (x <= 1)
    {
        z = 1-0.5*x-pow(x,2)+0.5*pow(x,3);
    }
    else if (x <= 2)
    {
        z = 1- 1.8333*x+pow(x,2)-0.1667*pow(x,3);
    }
    else
    {
        z = 0;
    }
    return z;
    
}

double W6(double y)
{
    double x = fabs(y);
    double z;
    if (x <= 1)
    {
        z = 1-0.3333*x - 1.2500*pow(x,2)+0.4166667*pow(x,3) + 0.2500*pow(x,4) - 0.0833*pow(x,5);
    }
    else if (x <= 2)
    {
        z = 1-1.0833*x - 0.6250*pow(x,2)+1.04166667*pow(x,3) - 0.3750*pow(x,4) + 0.04166667*pow(x,5);
    }
    else if (x <= 3)
    {
        z = 1-2.28333333*x + 1.8750*pow(x,2)-0.7083333*pow(x,3) + 0.125*pow(x,4) - 0.00833333*pow(x,5);
    }
    else
    {
        z = 0;
    }
    return z;
    
}


/* actual deposition algorithm */
void remapping(double *w, double *x, double *v, double *y, double *u, double *z, int n, int m)
{
    mwSize j;
    mwSize i;
    
    for (i=0; i<m; i++)
    {
        for (j=0; j<n; j++) 
        {
             z[i] += w[j]*W6(x[j]-y[i])*W6(v[j]-u[i]);
        }
    }
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
/* variable declarations here */
double *oldWeights;       /* 1xN input matrix */
double *oldPositions;       /* 1xN input matrix */
double *oldVelocities;       /* 1xN input matrix */
double *newPositions;       /* 1xN input matrix */
double *newVelocities;       /* 1xN input matrix */
mwSize newSize;           /* size of new matrix */
mwSize oldSize;           /* size of old matrix */
double *newWeights;      /* output matrix */
/* code here */

/* create a pointer to the real data in oldWeights  */
oldWeights = mxGetPr(prhs[0]);
/* create a pointer to the real data in oldPositions  */
oldPositions = mxGetPr(prhs[1]);
/* create a pointer to the real data in the oldVelocities  */
oldVelocities = mxGetPr(prhs[2]);
/* create a pointer to the real data in the newPositions  */
newPositions = mxGetPr(prhs[3]);
/* create a pointer to the real data in the oldVelocities  */
newVelocities = mxGetPr(prhs[4]);
/* get dimensions of the input matrix */
oldSize = mxGetN(prhs[1]);
/* get dimensions of the input matrix */
newSize = mxGetN(prhs[3]);
/* create the output matrix */
plhs[0] = mxCreateDoubleMatrix(1, newSize, mxREAL);
/* get a pointer to the real data in the output matrix newWeights */
newWeights = mxGetPr(plhs[0]);

/* call the computational routine */
remapping( oldWeights, oldPositions, oldVelocities, newPositions, newVelocities,newWeights, oldSize, newSize);

}