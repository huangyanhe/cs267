#include "decomposition.h"

void malloc_decomp(decomp* a_decomp){
    int M = ceil(size/RegionSize);
    a_decomp->M = M;
    a_decomp->N = M;
    a_decomp->Num_region = M*M;
    a_decomp->region_list = (region*)malloc(a_decomp->Num_region*sizeof(region));
    for(int i = 0; i < a_decomp->Num_region; i++){
        a_decomp->region_list[i].Num = 0;
        a_decomp->region_list[i].Capacity = InitialCapacity;
        a_decomp->region_list[i].ind = (int*)malloc(InitialCapacity*sizeof(int));
    }
}

void initial_decomp(double num_particle, particle_t* particles, decomp* a_decomp){
    for(int i = 0; i < a_decomp->Num_region; i++){
        a_decomp->region_list[i].Num = 0;
    }
    for(int i = 0; i < num_particle; i++){
        int m = particles[i].x/RegionSize, n = particles[i].y/RegionSize;
        add_particle(i, region_indexing(m,n,a_decomp), a_decomp);
    }
}

void free_decomp(decomp* a_decomp){
    for(int i = 0; i < a_decomp->Num_region; i++){
        free(a_decomp->region_list[i].ind);
    }
    free(a_decomp->region_list);
}

int region_indexing(int i, int j, decomp* a_decomp){
    if(i < 0 || j < 0 || i >= a_decomp->M || j>= a_decomp->N){
        printf("indices excess the boundary\n");
        printf("i = %d, j = %d\n", i,j);
        abort();
    }
    return i+j*a_decomp->M;
}

void add_particle(int particle_ind, int region_ind, decomp* a_decomp){
    region& temp = a_decomp->region_list[region_ind];
    if(temp.Num == temp.Capacity){
        temp.Capacity+=InitialCapacity;
        temp.ind = (int*)realloc(temp.ind, temp.Capacity*sizeof(int));
    }
    temp.ind[temp.Num] = particle_ind;
    temp.Num++;
}

void remove_particle(int particle_ind, int region_ind, decomp* a_decomp){
    region& temp = a_decomp->region_list[region_ind];
    for(int i = 0; i < temp.Num; i++){
        if(temp.ind[i] == particle_ind){
            for(int j = i; j < temp.Num-1; j++){
                temp.ind[j] = temp.ind[j+1];
            }
            break;
        }
        if(i == temp.Num-1){
            printf("Fail to remove particle %d from region %d\n", particle_ind, region_ind);
            abort();
        }
    }
    temp.Num--;
}

void check(int num_particle, decomp* a_decomp){
    int test = 0;
    for(int i = 0;i < a_decomp->Num_region; i++){
        test += a_decomp->region_list[i].Num;
    }
    if(test != num_particle){
        printf("wrong total number");
        abort();
    }
}
