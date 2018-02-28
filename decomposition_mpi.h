#ifndef _DECOMPOSITION_MPI_H_
#define _DECOMPOSITION_MPI_H_

#include <vector>
#include <mpi.h>
#include "common.h"

class decomp_proc{
    public:
        int proc_m, proc_n; // arrangement of processors
        int local_rank;
        int rank_m, rank_n; // position of this processor
        double block_size_m, block_size_n; // size of the domain in a single processor
        int num_region_m, num_region_n; // number of small bins in a single processor
        double region_size_m, region_size_n; // size of a single bin
        int interior_l_m, interior_u_m, interior_l_n, interior_u_n; // bins that are updated by this proc
        int closure_l_m, closure_u_m, closure_l_n, closure_u_n; // bins that are saved in this proc
        int lda; // leading dimension of regoin_list(bins)
        int max_recvcount; // maximum value of the recv count ======================
        std::vector<int> local_ind; // all the particles inside the interoir
        std::vector<std::vector<int> > region_list; // the particles in every bins
        std::vector<std::vector<int> > out_ind; // the particles that should be sent out
        double* sendbufs[8]; // sending buffers
        double* recvbufs[8]; // receiving buffers
        MPI_Status stats[16]; // dummy variable
        MPI_Request requests[16]; // request needed by non-blocking and waitall

        decomp_proc(int num_proc, int rank, int num_particle, particle_t* particles);
        ~decomp_proc();

        std::vector<int>& operator()(int i, int j);
        void clean_boundary();
        void delete_particle(int m_ind, int n_ind, int index);
        void local_delete_particle(int local_index);
        void pre_sync(int ind, int m_ind, int n_ind);
        void synchronization(particle_t* particles, int time_step);
        void write_sendbuf(int rank, int step, double* sendbuf, std::vector<int>& send_ind, particle_t* particles);
        void read_recvbuf(double* recvbuf, particle_t* particles);
        //void check(int num_particle, particle_t* particles);
};

#endif
