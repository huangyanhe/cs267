#ifndef _DECOMPOSITION_H_
#define _DECOMPOSITION_H_

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"


typedef struct{
    int* ind;
    int Num;
    int Capacity;
}region;

//typedef struct{
    //region* region_list;
    //int M, N, Num_region;
//}decomp;

//void malloc_decomp(decomp* a_decomp);
//void initial_decomp(double num_particle, particle_t* particles, decomp* a_decomp);
//void free_decomp(decomp* a_decomp);

//int region_indexing(int i, int j, decomp* a_decomp);

//void add_particle(int particle_ind, int region_ind, decomp* a_decomp);
//void remove_particle(int particle_ind, int region_ind, decomp* a_decomp);

//void check(int num_particle, decomp* a_decomp);
//
class decomp{
    public:
        region* region_list;
        int M, N, Num_region;

        decomp(double num_particle, particle_t* particles);
        ~decomp();

        region& operator()(int i, int j);
        void add_particle(int particle_ind, int i, int j);
        void remove_particle(int particle_ind, int i, int j);

        void check(int num_particle);
};

#endif
