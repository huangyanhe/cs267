#ifndef _DECOMPOSITION_H_
#define _DECOMPOSITION_H_

#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"


class region{
    public:
        int* ind;
        int Num;
        int Capacity;
};

class decomp{
    public:
        region* region_list;
        int* region_length;
        int M, N, Num_region;

        decomp(int num_particle, particle_t* particles);
        void init(int num_particle, particle_t* particles);
        ~decomp();

        std::vector<int > grid_M, grid_N;
        int num_sub_M, num_sub_N, Num_sub;

        region& operator()(int i, int j);

        void add_particle(int particle_ind, int i, int j);
        void delete_particle(int particle_ind, int i, int j);
        void check(int num_particle, particle_t* particles);
        void sub_decomp(int numthreads);
};

#endif
