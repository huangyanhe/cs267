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
    //printf("rank = %d, rank_m = %d, rank_n = %d\n", rank, rank_m, rank_n); fflush(stdout);

    block_size_m = size/proc_m; block_size_n = size/proc_n;

    num_region_m = floor(block_size_m/cutoff);
    num_region_n = floor(block_size_n/cutoff);

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
    //printf("rank %d, closure_l_m = %d, closure_u_m = %d, closure_l_n = %d, closure_u_n = %d\n", rank, closure_l_m, closure_u_m, closure_l_n, closure_u_n); fflush(stdout);

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

    out_ind.resize(8);

    for(int i = 0; i < 8; i++){
        sendbufs[i] = (double*) malloc(sizeof(double)*max_recvcount);
        //sendbufs[i][0] = 0.0;
        //if(rank == 0){
            //printf("sendbufs[i][0] = %f\n", sendbufs[i][0]); fflush(stdout);
        //}
        recvbufs[i] = (double*) malloc(sizeof(double)*max_recvcount);
    }
}
//checked

decomp_proc::~decomp_proc(){
    for(int i = 0; i < 8; i++){
        free(sendbufs[i]);
        free(recvbufs[i]);
    }
}
//checked

std::vector<int>& decomp_proc::operator()(int i, int j){
    int ind_m = i - closure_l_m, ind_n = j - closure_l_n;
    return region_list[ind_m+ind_n*lda];
}
//checked

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
//checked

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

void decomp_proc::pre_sync(int ind, int m_ind, int n_ind){
    if(m_ind >= interior_l_m+1 && m_ind < interior_u_m-1){
        if(n_ind >= interior_l_n+1 && n_ind < interior_u_n-1){
        }
        else if( n_ind < interior_l_n+1){
            out_ind[2].push_back(ind);
        }
        else{
            out_ind[6].push_back(ind);
        }
    }
    else if(n_ind >= interior_l_n+1 && n_ind < interior_u_n-1){
        if(m_ind < interior_l_m+1){
            out_ind[4].push_back(ind);
        }
        else{
            out_ind[0].push_back(ind);
        }
    }
    else if(n_ind < interior_l_n+1 && n_ind >= interior_l_n-1){
        if(m_ind < interior_l_m-1){
            out_ind[3].push_back(ind);
            out_ind[4].push_back(ind);
        }
        else if(m_ind >= interior_l_m-1 && m_ind < interior_l_m+1){
            out_ind[2].push_back(ind);
            out_ind[3].push_back(ind);
            out_ind[4].push_back(ind);
        }
        else if(m_ind >= interior_u_m-1 && m_ind < interior_u_m+1){
            out_ind[0].push_back(ind);
            out_ind[1].push_back(ind);
            out_ind[2].push_back(ind);
        }
        else{
            out_ind[0].push_back(ind);
            out_ind[1].push_back(ind);
        }
    }
    else if(n_ind >= interior_u_n-1 && n_ind < interior_u_n+1){
        if(m_ind < interior_l_m-1){
            out_ind[4].push_back(ind);
            out_ind[5].push_back(ind);
        }
        else if(m_ind >= interior_l_m-1 && m_ind < interior_l_m+1){
            out_ind[4].push_back(ind);
            out_ind[5].push_back(ind);
            out_ind[6].push_back(ind);
        }
        else if(m_ind >= interior_u_m-1 && m_ind < interior_u_m+1){
            out_ind[6].push_back(ind);
            out_ind[7].push_back(ind);
            out_ind[0].push_back(ind);
        }
        else{
            out_ind[7].push_back(ind);
            out_ind[0].push_back(ind);
        }
    }
    else if(n_ind >= interior_u_n+1){
        if(m_ind < interior_l_m-1){
            out_ind[5].push_back(ind);
        }
        else if(m_ind >= interior_l_m-1 && m_ind < interior_l_m+1){
            out_ind[5].push_back(ind);
            out_ind[6].push_back(ind);
        }
        else if(m_ind >= interior_u_m-1 && m_ind < interior_u_m+1){
            out_ind[6].push_back(ind);
            out_ind[7].push_back(ind);
        }
        else{
            out_ind[7].push_back(ind);
        }
    }
    else{
         if(m_ind < interior_l_m-1){
            out_ind[3].push_back(ind);
        }
        else if(m_ind >= interior_l_m-1 && m_ind < interior_l_m+1){
            out_ind[3].push_back(ind);
            out_ind[2].push_back(ind);
        }
        else if(m_ind >= interior_u_m-1 && m_ind < interior_u_m+1){
            out_ind[2].push_back(ind);
            out_ind[1].push_back(ind);
        }
        else{
            out_ind[1].push_back(ind);
        }
    }
}
// checked

void decomp_proc::synchronization(particle_t* particles, int time_step){
    //int see_rank = 0;
    //if(rank_m + rank_n*proc_m == see_rank && time_step == 1){
        //printf("rank %d step 1 synchronization starts\n", see_rank); fflush(stdout);
    //}
    int indexcount = 0;
    //if(rank_m + rank_n*proc_m == see_rank && time_step == 1){
        //printf("0\n"); fflush(stdout);
            //for(int i = 0; i < out_ind[0].size(); i++){
                //printf("%d\n", out_ind[0][i]);
            //}
    //}
    if(rank_m < proc_m-1){
        //if(rank_m + rank_n*proc_m == see_rank && time_step == 1){
            //printf("oh\n"); fflush(stdout);
            //sendbufs[0][0] = double(out_ind[0].size());
            //printf("before write, max_recvcount = %d, out_ind[0].size() = %d\n", max_recvcount, out_ind[0].size()); fflush(stdout);
        //}
        write_sendbuf(local_rank, time_step, sendbufs[0], out_ind[0], particles);
        //if(rank_m + rank_n*proc_m == 2 && time_step == 1){
            //printf("after write\n"); fflush(stdout);
        //}
        MPI_Isend(sendbufs[0], out_ind[0].size()*5+1, MPI_DOUBLE,
                (rank_m+1)+rank_n*proc_m, time_step, MPI_COMM_WORLD, &requests[indexcount]);
        MPI_Irecv(recvbufs[0], max_recvcount, MPI_DOUBLE,
                (rank_m+1)+rank_n*proc_m, time_step, MPI_COMM_WORLD, &requests[indexcount+1]);
        indexcount += 2;
    }
    //if(rank_m + rank_n*proc_m == see_rank && time_step == 1){
        //printf("1\n"); fflush(stdout);
    //}
    if(rank_m < proc_m-1 && rank_n > 0){
        write_sendbuf(local_rank, time_step, sendbufs[1], out_ind[1], particles);
        MPI_Isend(sendbufs[1], out_ind[1].size()*5+1, MPI_DOUBLE,
                (rank_m+1)+(rank_n-1)*proc_m, time_step, MPI_COMM_WORLD, &requests[indexcount]);
        MPI_Irecv(recvbufs[1], max_recvcount, MPI_DOUBLE,
                (rank_m+1)+(rank_n-1)*proc_m, time_step, MPI_COMM_WORLD, &requests[indexcount+1]);
        indexcount += 2;
    }
    //if(rank_m + rank_n*proc_m == see_rank && time_step == 1){
        //printf("2\n"); fflush(stdout);
    //}
    if(rank_n > 0){
        write_sendbuf(local_rank, time_step, sendbufs[2], out_ind[2], particles);
        MPI_Isend(sendbufs[2], out_ind[2].size()*5+1, MPI_DOUBLE,
                rank_m+(rank_n-1)*proc_m, time_step, MPI_COMM_WORLD, &requests[indexcount]);
        MPI_Irecv(recvbufs[2], max_recvcount, MPI_DOUBLE,
                rank_m+(rank_n-1)*proc_m, time_step, MPI_COMM_WORLD, &requests[indexcount+1]);
        indexcount += 2;
    }
    //if(rank_m + rank_n*proc_m == see_rank && time_step == 1){
        //printf("3\n"); fflush(stdout);
    //}
    if(rank_m>0 && rank_n>0){
        write_sendbuf(local_rank, time_step, sendbufs[3], out_ind[3], particles);
        MPI_Isend(sendbufs[3], out_ind[3].size()*5+1, MPI_DOUBLE,
                (rank_m-1)+(rank_n-1)*proc_m, time_step, MPI_COMM_WORLD, &requests[indexcount]);
        MPI_Irecv(recvbufs[3], max_recvcount, MPI_DOUBLE,
                (rank_m-1)+(rank_n-1)*proc_m, time_step, MPI_COMM_WORLD, &requests[indexcount+1]);
        indexcount += 2;
    }
    //if(rank_m + rank_n*proc_m == see_rank && time_step == 1){
        //printf("4\n"); fflush(stdout);
    //}
    if(rank_m>0){
        write_sendbuf(local_rank, time_step, sendbufs[4], out_ind[4], particles);
        MPI_Isend(sendbufs[4], out_ind[4].size()*5+1, MPI_DOUBLE,
                (rank_m-1)+rank_n*proc_m, time_step, MPI_COMM_WORLD, &requests[indexcount]);
        MPI_Irecv(recvbufs[4], max_recvcount, MPI_DOUBLE,
                (rank_m-1)+rank_n*proc_m, time_step, MPI_COMM_WORLD, &requests[indexcount+1]);
        indexcount += 2;
    }
    //if(rank_m + rank_n*proc_m == see_rank && time_step == 1){
        //printf("5\n"); fflush(stdout);
    //}
    if(rank_m>0 && rank_n<proc_n-1){
        write_sendbuf(local_rank, time_step, sendbufs[5], out_ind[5], particles);
        MPI_Isend(sendbufs[5], out_ind[5].size()*5+1, MPI_DOUBLE,
                (rank_m-1)+(rank_n+1)*proc_m, time_step, MPI_COMM_WORLD, &requests[indexcount]);
        MPI_Irecv(recvbufs[5], max_recvcount, MPI_DOUBLE,
                (rank_m-1)+(rank_n+1)*proc_m, time_step, MPI_COMM_WORLD, &requests[indexcount+1]);
        indexcount += 2;
    }
    //if(rank_m + rank_n*proc_m == see_rank && time_step == 1){
        //printf("6\n"); fflush(stdout);
    //}
    if(rank_n<proc_n-1){
        write_sendbuf(local_rank, time_step, sendbufs[6], out_ind[6], particles);
        MPI_Isend(sendbufs[6], out_ind[6].size()*5+1, MPI_DOUBLE,
                rank_m+(rank_n+1)*proc_m, time_step, MPI_COMM_WORLD, &requests[indexcount]);
        MPI_Irecv(recvbufs[6], max_recvcount, MPI_DOUBLE,
                rank_m+(rank_n+1)*proc_m, time_step, MPI_COMM_WORLD, &requests[indexcount+1]);
        indexcount += 2;
    }
    //if(rank_m + rank_n*proc_m == see_rank && time_step == 1){
        //printf("7\n"); fflush(stdout);
    //}
    if(rank_m<proc_m-1 && rank_n<proc_n-1){
        write_sendbuf(local_rank, time_step, sendbufs[7], out_ind[7], particles);
        MPI_Isend(sendbufs[7], out_ind[7].size()*5+1, MPI_DOUBLE,
                (rank_m+1)+(rank_n+1)*proc_m, time_step, MPI_COMM_WORLD, &requests[indexcount]);
        MPI_Irecv(recvbufs[7], max_recvcount, MPI_DOUBLE,
                (rank_m+1)+(rank_n+1)*proc_m, time_step, MPI_COMM_WORLD, &requests[indexcount+1]);
        indexcount += 2;
    }
    //if(rank_m + rank_n*proc_m == see_rank && time_step == 1){
        //printf("8\n"); fflush(stdout);
    //}

    //if(rank_m + rank_n*proc_m == see_rank && time_step == 1){
        //printf("rank %d step 1 synchronization before waitall\n", see_rank); fflush(stdout);
    //}

    MPI_Waitall(indexcount, requests, MPI_STATUSES_IGNORE);
    //if(rank_m + rank_n*proc_m == see_rank && time_step == 1){
        //printf("rank %d step 1 synchronization after waitall\n", see_rank); fflush(stdout);
    //}

    if(rank_m < proc_m-1){
        read_recvbuf(recvbufs[0], particles);
    }
    if(rank_m < proc_m-1 && rank_n > 0){
        read_recvbuf(recvbufs[1], particles);
    }
    if(rank_n > 0){
        read_recvbuf(recvbufs[2], particles);
    }
    if(rank_m>0 && rank_n>0){
        read_recvbuf(recvbufs[3], particles);
    }
    if(rank_m>0){
        read_recvbuf(recvbufs[4], particles);
    }
    if(rank_m>0 && rank_n<proc_n-1){
        read_recvbuf(recvbufs[5], particles);
    }
    if(rank_n<proc_n-1){
        read_recvbuf(recvbufs[6], particles);
    }
    if(rank_m<proc_m-1 && rank_n<proc_n-1){
        read_recvbuf(recvbufs[7], particles);
    }

    //if(rank_m + rank_n*proc_m == see_rank && time_step == 1){
        //printf("rank %d step 1 synchronization after read recvbufs\n", see_rank); fflush(stdout);
    //}

    for(int i = 0; i < 8; i++){
        out_ind[i].clear();
    }
}
// checked

void decomp_proc::write_sendbuf(int rank, int step, double* sendbuf, std::vector<int>& send_ind, particle_t* particles){
    //if(rank == 0 && step == 1){
        //printf("rank 0 step 1 write_sendbuf starts\n"); fflush(stdout);
        //printf("send_ind.size() = %d\n", send_ind.size()); fflush(stdout);
        //printf("one more = %f\n", sendbuf[0]); fflush(stdout);
    //}
    sendbuf[0] = double(send_ind.size());
    //if(rank == 0 && step == 1){
        //printf("rank 0 step 1 after set the total num\n"); fflush(stdout);
    //}
    for(int i = 0; i < send_ind.size(); i++){
        //if(rank == 0 && step == 1){
            //printf("the %d iteration in write_sendbuf, before get the index\n", i); fflush(stdout);
        //}
        sendbuf[5*i+1] = double(send_ind[i]);
        //if(rank == 0 && step == 1){
            //printf("the %d iteration in write_sendbuf, get the index, send_ind[%d] = %d\n", i, i, send_ind[i]); fflush(stdout);
        //}
        particle_t& temp = particles[send_ind[i]];
        sendbuf[5*i+2] = temp.x;
        sendbuf[5*i+3] = temp.y;
        sendbuf[5*i+4] = temp.vx;
        sendbuf[5*i+5] = temp.vy;
    }
}
//checked

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
//checked
