#ifndef _DECOMPOSITION_H_
#define _DECOMPOSITION_H_

#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"

class decomp{
    public:
        //region* region_list;
        std::vector< std::vector<int> > region_list;
        std::vector< int >region_length;
        int M, N, Num_region;

        std::vector<int> grid_M, grid_N;
        int num_sub_M, num_sub_N, Num_sub;

        decomp(double num_particle, particle_t* particles);
        void init(double num_particle, particle_t* particles);
        //~decomp();

        std::vector<int>& operator()(int i, int j);

        void add_particle(int particle_ind, int i, int j);
        void delete_particle(int particle_ind, int i, int j);
        void check(int num_particle, particle_t* particles);

        void malloc_sub_decomp(int numthreads);
};

#endif
