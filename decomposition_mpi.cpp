#include "decomposition_mpi.h"
#include <cmath>
#include <algorithm>
#include <cassert>

#define cutoff  0.01
extern double size;

decomp_proc::decomp_proc(int num_proc, int rank, int num_particle, particle_t* particles){
    local_rank = rank;
    for(int i = sqrt(num_proc); i > 0; i--){
        if(num_proc%i == 0){
            proc_m = i;
            break;
        }
    }
    proc_n = num_proc / proc_m;

    rank_m = rank % proc_m; rank_n = rank / proc_m;

    double block_size_m = size/proc_m, block_size_n = size/proc_n;

    int num_region_m = floor(block_size_m/cutoff), num_region_n = floor(block_size_n/cutoff);

    region_size_m = block_size_m/num_region_m;
    region_size_n = block_size_n/num_region_n;

    interior_l_m = rank_m*num_region_m;
    interior_u_m = (rank_m+1)*num_region_m;
    interior_l_n = rank_n*num_region_n;
    interior_u_n = (rank_n+1)*num_region_n;

    closure_l_m = (rank_m == 0)? interior_l_m: interior_l_m-1;
    closure_u_m = (rank_m == proc_m-1)? interior_u_m: interior_u_m+1;
    closure_l_n = (rank_n == 0)? interior_l_n: interior_l_n-1;
    closure_u_n = (rank_n == proc_n-1)? interior_u_n: interior_u_n+1;

    lda = closure_u_m - closure_l_m;

    max_recvcount = 5*num_particle+1;

    region_list.resize((closure_u_m - closure_l_m)*(closure_u_n - closure_l_n));
    for(int i = 0; i < num_particle; i++){
        int ind_m = particles[i].x/region_size_m,
            ind_n = particles[i].y/region_size_n;
        if(ind_m >= closure_l_m && ind_m < closure_u_m
                && ind_n >= closure_l_n && ind_n < closure_u_n){
            (*this)(ind_m, ind_n).push_back(i);
            if(ind_m >= interior_l_m && ind_m < interior_u_m
                    && ind_n >= interior_l_n && ind_n < interior_u_n){
                local_ind.push_back(i);
            }
        }
    }

    for(int i = 0; i < 8; i++){
        sendbufs[i] = (double*) malloc(sizeof(double)*max_recvcount);
        sendbufs[i][0] = 0.0;
        recvbufs[i] = (double*) malloc(sizeof(double)*max_recvcount);
    }

    neighbor_exist.resize(8);
    neighbor_exist[0] = (rank_m < proc_m-1);
    neighbor_exist[1] = (rank_m < proc_m-1 && rank_n >0);
    neighbor_exist[2] = (rank_n > 0);
    neighbor_exist[3] = (rank_m > 0 && rank_n > 0);
    neighbor_exist[4] = (rank_m > 0);
    neighbor_exist[5] = (rank_m > 0 && rank_n < proc_n-1);
    neighbor_exist[6] = (rank_n < proc_n-1);
    neighbor_exist[7] = (rank_n < proc_n-1 && rank_m < proc_m-1);

    neighbor_shift.resize(16);
    neighbor_shift[0] = 1;
    neighbor_shift[1] = 0;
    neighbor_shift[2] = 1;
    neighbor_shift[3] = -1;
    neighbor_shift[4] = 0;
    neighbor_shift[5] = -1;
    neighbor_shift[6] = -1;
    neighbor_shift[7] = -1;
    neighbor_shift[8] = -1;
    neighbor_shift[9] = 0;
    neighbor_shift[10] = -1;
    neighbor_shift[11] = 1;
    neighbor_shift[12] = 0;
    neighbor_shift[13] = 1;
    neighbor_shift[14] = 1;
    neighbor_shift[15] = 1;
}

decomp_proc::~decomp_proc(){
    for(int i = 0; i < 8; i++){
        free(sendbufs[i]);
        free(recvbufs[i]);
    }
}

std::vector<int>& decomp_proc::operator()(int i, int j){
    int ind_m = i - closure_l_m, ind_n = j - closure_l_n;
    return region_list[ind_m+ind_n*lda];
}

void decomp_proc::clean_boundary(){
    for(int m = closure_l_m; m < closure_u_m; m++){
        for(int n = closure_l_n; n < closure_u_n; n++){
            if(m < interior_l_m || m >= interior_u_m ||
                    n < interior_l_n || n >= interior_u_n){
                (*this)(m,n).clear();
            }
        }
    }
}

void decomp_proc::delete_particle(int m_ind, int n_ind, int index){
    std::vector<int>::iterator it = find((*this)(m_ind, n_ind).begin(), (*this)(m_ind, n_ind).end(), index);
    if(it == (*this)(m_ind, n_ind).end()){
        printf("fail to find the index in delete_particle\n"); fflush(stdout);
        exit(1);
    }
    (*this)(m_ind,n_ind).erase(it);
}

void decomp_proc::local_delete_particle(int local_index){
    local_ind[local_index] = local_ind.back();
    local_ind.pop_back();
}

void decomp_proc::pre_sync(int ind, int m_ind, int n_ind, particle_t* particles){
    if(m_ind >= interior_l_m+1 && m_ind < interior_u_m-1){
        if(n_ind >= interior_l_n+1 && n_ind < interior_u_n-1){
        }
        else if( n_ind < interior_l_n+1){
            sendbuf_push_back(2, ind, particles);
        }
        else{
            sendbuf_push_back(6, ind, particles);
        }
    }
    else if(n_ind >= interior_l_n+1 && n_ind < interior_u_n-1){
        if(m_ind < interior_l_m+1){
            sendbuf_push_back(4, ind, particles);
        }
        else{
            sendbuf_push_back(0, ind, particles);
        }
    }
    else if(n_ind < interior_l_n+1 && n_ind >= interior_l_n-1){
        if(m_ind < interior_l_m-1){
            sendbuf_push_back(3, ind, particles);
            sendbuf_push_back(4, ind, particles);
        }
        else if(m_ind >= interior_l_m-1 && m_ind < interior_l_m+1){
            sendbuf_push_back(2, ind, particles);
            sendbuf_push_back(3, ind, particles);
            sendbuf_push_back(4, ind, particles);
        }
        else if(m_ind >= interior_u_m-1 && m_ind < interior_u_m+1){
            sendbuf_push_back(0, ind, particles);
            sendbuf_push_back(1, ind, particles);
            sendbuf_push_back(2, ind, particles);
        }
        else{
            sendbuf_push_back(0, ind, particles);
            sendbuf_push_back(1, ind, particles);
        }
    }
    else if(n_ind >= interior_u_n-1 && n_ind < interior_u_n+1){
        if(m_ind < interior_l_m-1){
            sendbuf_push_back(4, ind, particles);
            sendbuf_push_back(5, ind, particles);
        }
        else if(m_ind >= interior_l_m-1 && m_ind < interior_l_m+1){
            sendbuf_push_back(4, ind, particles);
            sendbuf_push_back(5, ind, particles);
            sendbuf_push_back(6, ind, particles);
        }
        else if(m_ind >= interior_u_m-1 && m_ind < interior_u_m+1){
            sendbuf_push_back(6, ind, particles);
            sendbuf_push_back(7, ind, particles);
            sendbuf_push_back(0, ind, particles);
        }
        else{
            sendbuf_push_back(7, ind, particles);
            sendbuf_push_back(0, ind, particles);
        }
    }
    else if(n_ind >= interior_u_n+1){
        if(m_ind < interior_l_m-1){
            sendbuf_push_back(5, ind, particles);
        }
        else if(m_ind >= interior_l_m-1 && m_ind < interior_l_m+1){
            sendbuf_push_back(5, ind, particles);
            sendbuf_push_back(6, ind, particles);
        }
        else if(m_ind >= interior_u_m-1 && m_ind < interior_u_m+1){
            sendbuf_push_back(6, ind, particles);
            sendbuf_push_back(7, ind, particles);
        }
        else{
            sendbuf_push_back(7, ind, particles);
        }
    }
    else{
         if(m_ind < interior_l_m-1){
            sendbuf_push_back(3, ind, particles);
        }
        else if(m_ind >= interior_l_m-1 && m_ind < interior_l_m+1){
            sendbuf_push_back(3, ind, particles);
            sendbuf_push_back(2, ind, particles);
        }
        else if(m_ind >= interior_u_m-1 && m_ind < interior_u_m+1){
            sendbuf_push_back(2, ind, particles);
            sendbuf_push_back(1, ind, particles);
        }
        else{
            sendbuf_push_back(1, ind, particles);
        }
    }
}

void decomp_proc::synchronization(particle_t* particles, int time_step){
    int indexcount = 0;
    for(int i = 0; i < 8; i++){
        if(neighbor_exist[i]){
            MPI_Isend(sendbufs[i], int(sendbufs[i][0])*5+1, MPI_DOUBLE,
                    (rank_m+neighbor_shift[2*i])+(rank_n+neighbor_shift[2*i+1])*proc_m,
                    time_step, MPI_COMM_WORLD, &requests[indexcount]);
            MPI_Irecv(recvbufs[i], max_recvcount, MPI_DOUBLE,
                    (rank_m+neighbor_shift[2*i])+(rank_n+neighbor_shift[2*i+1])*proc_m,
                    time_step, MPI_COMM_WORLD, &requests[indexcount+1]);
            indexcount += 2;
        }
    }
    MPI_Waitall(indexcount, requests, MPI_STATUSES_IGNORE);

    for(int i = 0; i < 8; i++){
        if(neighbor_exist[i]){
            read_recvbuf(recvbufs[i], particles);
        }
    }

    for(int i = 0; i < 8; i++){
        sendbufs[i][0] = 0.0;
    }
}

void decomp_proc::sendbuf_push_back(int neighbor_ind, int ind, particle_t* particles){
    int total_num = sendbufs[neighbor_ind][0];
    sendbufs[neighbor_ind][5*total_num+1] = double(ind);
    sendbufs[neighbor_ind][5*total_num+2] = particles[ind].x;
    sendbufs[neighbor_ind][5*total_num+3] = particles[ind].y;
    sendbufs[neighbor_ind][5*total_num+4] = particles[ind].vx;
    sendbufs[neighbor_ind][5*total_num+5] = particles[ind].vy;
    sendbufs[neighbor_ind][0] += 1.0;
}

//void decomp_proc::write_sendbuf(int rank, int step, double* sendbuf, std::vector<int>& send_ind, particle_t* particles){
    //sendbuf[0] = double(send_ind.size());
    //for(int i = 0; i < send_ind.size(); i++){
        //sendbuf[5*i+1] = double(send_ind[i]);
        //particle_t& temp = particles[send_ind[i]];
        //sendbuf[5*i+2] = temp.x;
        //sendbuf[5*i+3] = temp.y;
        //sendbuf[5*i+4] = temp.vx;
        //sendbuf[5*i+5] = temp.vy;
    //}
//}

void decomp_proc::read_recvbuf(double* recvbuf, particle_t* particles){
    int recvcount = recvbuf[0];
    for(int i = 0; i < recvcount; i++){
        particle_t& temp = particles[int(recvbuf[5*i+1])];
        temp.x = recvbuf[5*i+2];
        temp.y = recvbuf[5*i+3];
        temp.vx = recvbuf[5*i+4];
        temp.vy = recvbuf[5*i+5];
        int m = temp.x/region_size_m, n = temp.y/region_size_n;
        (*this)(m,n).push_back(int(recvbuf[5*i+1]));
        if(m >= interior_l_m && m < interior_u_m
                && n >= interior_l_n && n < interior_u_n){
            local_ind.push_back(recvbuf[5*i+1]);
        }
    }
}
