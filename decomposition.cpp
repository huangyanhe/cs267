#include "decomposition.h"

extern double size;
extern double RegionSize;
extern int InitialCapacity;

decomp::decomp(double num_particle, particle_t* particles){
    int required = ceil(size/RegionSize);
    M = required;
    N = required;
    Num_region = M*N;
    region_list = (region*)malloc(Num_region*sizeof(region));
    for(int i = 0; i < Num_region; i++){
        region_list[i].Num = 0;
        region_list[i].Capacity = InitialCapacity;
        region_list[i].ind = (int*)malloc(InitialCapacity*sizeof(int));
    }

    for(int i = 0; i < num_particle; i++){
        int m = particles[i].x/RegionSize, n = particles[i].y/RegionSize;
        add_particle(i, m, n);
    }
}

decomp::~decomp(){
    for(int i = 0; i < Num_region; i++){
        free(region_list[i].ind);
    }
    free(region_list);
}

region& decomp::operator()(int i, int j){
    if(i < 0 || j < 0 || i >= M || j>= N){
        printf("indices excess the boundary\n");
        printf("i = %d, j = %d\n", i,j);
        abort();
    }
    return region_list[i+j*M];
}

void decomp::add_particle(int particle_ind, int i, int j){
    region& temp = region_list[i+j*M];
    if(temp.Num == temp.Capacity){
        temp.Capacity+=InitialCapacity;
        temp.ind = (int*)realloc(temp.ind, temp.Capacity*sizeof(int));
    }
    temp.ind[temp.Num] = particle_ind;
    temp.Num++;
}

void decomp::remove_particle(int particle_ind, int i, int j){
    region& temp = region_list[i+j*M];
    for(int i = 0; i < temp.Num; i++){
        if(temp.ind[i] == particle_ind){
            for(int j = i; j < temp.Num-1; j++){
                temp.ind[j] = temp.ind[j+1];
            }
            break;
        }
        if(i == temp.Num-1){
            printf("Fail to remove particle %d from region (%d, %d)\n", particle_ind,i,j);
            abort();
        }
    }
    temp.Num--;
}

void decomp::check(int num_particle){
    int test = 0;
    for(int i = 0;i < Num_region; i++){
        test += region_list[i].Num;
    }
    if(test != num_particle){
        printf("wrong total number");
        abort();
    }
}
