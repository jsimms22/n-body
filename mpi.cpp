#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <signal.h>
#include <unistd.h>
#include "common.h"

using std::vector;
using std::map;
using std::set;

#define _cutoff 0.01    //Value copied from common.cpp
#define _density 0.0005

double bin_size, grid_size;
int bin_count;

inline void build_bins(vector<bin_t>& bins, particle_t* particles, int n)
{
    grid_size = sqrt(n * _density);
    bin_size = _cutoff;  
    bin_count = int(grid_size / bin_size) + 1; // Should be around sqrt(N/2)

    // printf("Grid Size: %.4lf\n", grid_size);
    // printf("Number of Bins: %d*%d\n", bin_count, bin_count);
    // printf("Bin Size: %.2lf\n", bin_size);
    // Increase\Decrease bin_count to be something like 2^k?
    
    bins.resize(bin_count * bin_count);

    for (int i = 0; i < n; i++)
    {
        int x = int(particles[i].x / bin_size);
        int y = int(particles[i].y / bin_size);
        bins[x*bin_count + y].push_back(particles[i]);
    }
}

inline void compute_forces_for_bin(vector<bin_t>& bins, int i, int j, double& dmin, double& davg, int& navg)
{
    bin_t& vec = bins[i * bin_count + j];

    for (int k = 0; k < vec.size(); k++)
        vec[k].ax = vec[k].ay = 0;

    for (int dx = -1; dx <= 1; dx++)   //Search over nearby 8 bins and itself
    {
        for (int dy = -1; dy <= 1; dy++)
        {
            if (i + dx >= 0 && i + dx < bin_count && j + dy >= 0 && j + dy < bin_count)
            {
                bin_t& vec2 = bins[(i+dx) * bin_count + j + dy];
                for (int k = 0; k < vec.size(); k++)
                    for (int l = 0; l < vec2.size(); l++)
                        apply_force( vec[k], vec2[l], &dmin, &davg, &navg);
            }
        }
    }
}

void bin_particle(particle_t& particle, vector<bin_t>& bins)
{
    int x = particle.x / bin_size;
    int y = particle.y / bin_size;
    bins[x*bin_count + y].push_back(particle);
}


inline void get_neighbors(int i, int j, vector<int>& neighbors)
{
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            if (dx == 0 && dy == 0)
                continue;
            if (i + dx >= 0 && i + dx < bin_count && j + dy >= 0 && j + dy < bin_count) {
                int index = (i + dx) * bin_count + j + dy;
                neighbors.push_back(index);
            }
        }
    }
}


//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    //signal(SIGSEGV, sigsegv);

    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg; 
 
    //
    //  process command line parameters
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;


    particle_t *particles = new particle_t[n];
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    
    set_size( n );
    if( rank == 0 )
        init_particles( n, particles );

    MPI_Bcast(particles, n, PARTICLE, 0, MPI_COMM_WORLD);

    vector<bin_t> bins;
    build_bins(bins, particles, n);

    delete[] particles;
    particles = NULL;

    int x_bins_per_proc = bin_count / n_proc;

    // although each worker has all particles, we only access particles within
    // my_bins_start, my_bins_end.

    int my_bins_start = x_bins_per_proc * rank;
    int my_bins_end = x_bins_per_proc * (rank + 1);

    if (rank == n_proc - 1)
        my_bins_end = bin_count;
    
    //  simulate a number of time steps
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;

        // compute local forces
        for (int i = my_bins_start; i < my_bins_end; ++i) {
            for (int j = 0; j < bin_count; ++j) {
                compute_forces_for_bin(bins, i, j, dmin, davg, navg);
            }
        }

        if (find_option( argc, argv, "-no" ) == -1) {
          MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
          if (rank == 0){
            if (rnavg) {
              absavg +=  rdavg/rnavg;
              nabsavg++;
            }
            if (rdmin < absmin) 
                absmin = rdmin;
          }
        }

        // move, but not rebin
        bin_t local_move;
        bin_t remote_move;

        for (int i = my_bins_start; i < my_bins_end; ++i) {
            for (int j = 0; j < bin_count; ++j) {
                bin_t& bin = bins[i * bin_count + j];
                int tail = bin.size(), k = 0;
                for (; k < tail; ) {
                    move(bin[k]);
                    int x = int(bin[k].x / bin_size);
                    int y = int(bin[k].y / bin_size);
                    if (my_bins_start <= x && x < my_bins_end) {
                        if (x == i && y == j)
                            ++k;
                        else {
                            local_move.push_back(bin[k]);
                            bin[k] = bin[--tail];
                        }
                    } else {
                        remote_move.push_back(bin[k]);
                        bin[k] = bin[--tail];
                    }
                }
                bin.resize(k);
            }
        }

        for (int i = 0; i < local_move.size(); ++i) {
            bin_particle(local_move[i], bins);
        }

        if (rank != 0) {
            for (int i = my_bins_start - 1, j = 0; j < bin_count; ++j) {
                bin_t& bin = bins[i * bin_count + j];
                bin.clear();
            }
            for (int i = my_bins_start, j = 0; j < bin_count; ++j) {
                bin_t& bin = bins[i * bin_count + j];
                remote_move.insert(remote_move.end(), bin.begin(), bin.end());
                bin.clear();
            }
        }

        if (rank != n_proc - 1) {
            for (int i = my_bins_end, j = 0; j < bin_count; ++j) {
                bin_t& bin = bins[i * bin_count + j];
                bin.clear();
            }
            for (int i = my_bins_end - 1, j = 0; j < bin_count; ++j) {
                bin_t& bin = bins[i * bin_count + j];
                remote_move.insert(remote_move.end(), bin.begin(), bin.end());
                bin.clear();
            }
        }

        bin_t incoming_move;
        int send_count = remote_move.size();
        int recv_counts[n_proc];

        MPI_Gather(&send_count, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

        int displs[n_proc];
        int total_num = 0;

        if (rank == 0) {
            displs[0] = 0;
            for (int i = 1; i < n_proc; ++i) {
                displs[i] = displs[i-1] + recv_counts[i-1];
            }
            total_num = recv_counts[n_proc-1] + displs[n_proc-1];
            incoming_move.resize(total_num);
        }

        // now root knows total_num.

        MPI_Gatherv(remote_move.data(), send_count, PARTICLE, 
            incoming_move.data(), recv_counts, displs, PARTICLE, 
            0, MPI_COMM_WORLD);

        vector<bin_t> scatter_particles;
        scatter_particles.resize(n_proc);

        if (rank == 0) {
            for (int i = 0; i < incoming_move.size(); ++i) {
                int x = int(incoming_move[i].x / bin_size);

                assert(incoming_move[i].x >= 0 && incoming_move[i].y >= 0 &&
                    incoming_move[i].x <= grid_size && incoming_move[i].y <= grid_size);

                int who = min(x / x_bins_per_proc, n_proc-1);
                scatter_particles[who].push_back(incoming_move[i]);

                int row = x % x_bins_per_proc;
                if (row == 0 && who != 0)
                    scatter_particles[who - 1].push_back(incoming_move[i]);
                if (row == x_bins_per_proc-1 && who != n_proc-1)
                    scatter_particles[who + 1].push_back(incoming_move[i]);
            }
            for (int i = 0; i < n_proc; ++i) {
                recv_counts[i] = scatter_particles[i].size();
            }
            displs[0] = 0;
            for (int i = 1; i < n_proc; ++i) {
                displs[i] = displs[i-1] + recv_counts[i-1];
            }
        }

        send_count = 0;
        MPI_Scatter(recv_counts, 1, MPI_INT, &send_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        bin_t outgoing_move;
        outgoing_move.resize(send_count);

        bin_t scatter_particles_flatten;
        for (int i = 0; i < scatter_particles.size(); ++i) {
            scatter_particles_flatten.insert(scatter_particles_flatten.end(),
                scatter_particles[i].begin(), scatter_particles[i].end());
        }

        MPI_Scatterv(scatter_particles_flatten.data(), recv_counts, displs, PARTICLE, 
            outgoing_move.data(), send_count, PARTICLE, 0, MPI_COMM_WORLD);

        // printf("worker: %d. Bin.\n", rank);
        for (int i = 0; i < send_count; ++i) {
            particle_t &p = outgoing_move[i];
            assert(p.x >= 0 && p.y >= 0 && p.x <= grid_size && p.y <= grid_size);
            bin_particle(p, bins);
        }
    }
    simulation_time = read_timer( ) - simulation_time;
  
    if (rank == 0) {  
      printf( "n = %d, simulation time = %g seconds", n, simulation_time);

      if( find_option( argc, argv, "-no" ) == -1 )
      {
        if (nabsavg) absavg /= nabsavg;
          // 
          //  -The minimum distance absmin between 2 particles during the run of the simulation
          //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
          //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
          //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
          //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
          //
          printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
          if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
          if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
      }
      printf("\n");     
        
      //  
      // Printing summary data
      //  
      if( fsum)
        fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
    }
  
    //
    //  release resources
    //
    if ( fsum )
        fclose( fsum );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize();
    
    return 0;
}
