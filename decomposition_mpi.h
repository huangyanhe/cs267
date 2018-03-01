#ifndef _DECOMPOSITION_MPI_H_
#define _DECOMPOSITION_MPI_H_

#include <vector>
#include <mpi.h>
#include "common.h"

class decomp_proc{
    public:
        int proc_m, proc_n;
        int local_rank;
        int rank_m, rank_n;
        //double block_size_m, block_size_n;
        //int num_region_m, num_region_n;
        double region_size_m, region_size_n;
        int interior_l_m, interior_u_m, interior_l_n, interior_u_n;
        int closure_l_m, closure_u_m, closure_l_n, closure_u_n;
        int lda;
        int max_recvcount;
        std::vector<int> local_ind;
        std::vector<std::vector<int> > region_list;
        //std::vector<std::vector<int> > out_ind;
        double* sendbufs[8];
        double* recvbufs[8];
        MPI_Request requests[16];
        std::vector<bool> neighbor_exist;
        std::vector<int> neighbor_shift;

        decomp_proc(int num_proc, int rank, int num_particle, particle_t* particles);
        ~decomp_proc();

        std::vector<int>& operator()(int i, int j);
        void clean_boundary();
        void delete_particle(int m_ind, int n_ind, int index);
        void local_delete_particle(int local_index);
        void pre_sync(int ind, int m_ind, int n_ind, particle_t* particles);
        void synchronization(particle_t* particles, int time_step);
        void sendbuf_push_back(int neighbor_ind, int ind, particle_t* particles);
        //void write_sendbuf(int rank, int step, double* sendbuf, std::vector<int>& send_ind, particle_t* particles);
        void read_recvbuf(double* recvbuf, particle_t* particles);
};

#endif
