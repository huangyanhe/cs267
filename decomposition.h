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
        int M, N, Num_region;

        decomp(double num_particle, particle_t* particles);
        void init(double num_particle, particle_t* particles);
        ~decomp();

        region& operator()(int i, int j);

        void add_particle(int particle_ind, int i, int j);
        void delete_particle(int particle_ind, int i, int j);
        void check(int num_particle, particle_t* particles);
};

#endif
