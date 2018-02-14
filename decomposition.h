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

class sub_decomp{
    public:
        int ul_m, ul_n;
        int rd_m, rd_n;
        std::vector<std::vector<int> > particle_removed;
};

class decomp{
    public:
        region* region_list;
        int M, Num_region;

        sub_decomp* sub_decomp_list;
        std::vector<int> grid_M, grid_N;
        int num_sub_M, num_sub_N, Num_sub;

        decomp(double num_particle, particle_t* particles);
        ~decomp();

        region& operator()(int i, int j);

        void add_particle(int particle_ind, int i, int j);
        void delete_particle(int particle_ind, int i, int j);
        void check(int num_particle);

        void malloc_sub_decomp(int numthreads);
        void free_sub_decomp();
};

#endif
