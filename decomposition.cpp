#include "decomposition.h"
#include <math.h>

extern double size;
extern int InitialCapacity;
extern double RegionSize;

decomp::decomp(double num_particle, particle_t* particles){
    int required = ceil(size/RegionSize);
    M = required;
    Num_region = M*M;
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
    if(i < 0 || j < 0 || i >= M || j>= M){
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

void decomp::delete_particle(int particle_ind, int i, int j){
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

void decomp::malloc_sub_decomp(int numthreads){
    Num_sub = numthreads;
    int num_sub_M = sqrt(numthreads);
    for(int i = num_sub_M; i > 0; i--){
        if(Num_sub%i == 0){
            num_sub_M = i;
            break;
        }
    }
    num_sub_N = M/num_sub_M;

    int size_M = M/num_sub_M, size_N = M/num_sub_N;
    int remain_M = M%num_sub_M, remain_N = M%um_sub_N;
    grid_M.clear(); grid_N.clear();
    int temp_grid = 0;
    for(int i = 0; i < num_sub_M; i++){
        grid_M.push_back(temp_grid);
        temp_grid += (i < remain_M)? size_M+1:size_M;
    }
    grid_M.push_back(M);
    temp_grid = 0;
    for(int i = 0; i < num_sub_N; i++){
        grid_N.push_back(temp_grid);
        temp_grid += (i < remain_N)? size_N+1:size_N;
    }
    grid_N.push_back(M);
    particle_removed.resize(Num_sub*Num_sub);
}

