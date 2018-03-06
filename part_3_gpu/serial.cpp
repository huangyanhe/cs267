#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <tuple>
#include <list>
#include "common_gsi.h"

// Ordered pair for debugging
// a's acceleration is affected by interaction with b
struct pair_t {
  int a, b;

  pair_t(int a, int b) {
    this->a = a;
    this->b = b;
  }

  bool operator==(const pair_t &y) const {
    if (this->a == y.a && this->b == y.b) {
      return true;
    } else {
      return false;
    }
  }
};

//
//  benchmarking program
//
#define cutoff 0.01
int main( int argc, char **argv )
{
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

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

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    double size = set_size(n);
    //double cutoff = get_cutoff();

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    init_particles( n, particles );

    double bin_size = cutoff;

    int bindim = ((int) (ceil(size / bin_size) + 0.5));
    int nbins = bindim*bindim;

    //printf("grid dimensions are %d -> %lf and cutoff is %lf, so %d bins of size %lf\n", 0, size, cutoff, nbins, bin_size);//fflush(stdout);

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    std::vector <std::list <int> > bins(nbins, std::list <int> ());

    // Bin all the particles
    for (int i = 0; i < n; i++) {
        int xbin = (int) (particles[i].x / bin_size);
        int ybin = (int) (particles[i].y / bin_size);
        bins[xbin + ybin*bindim].push_back(i);
    }

    for( int step = 0; step < NSTEPS; step++ )
    {
	navg = 0;
        davg = 0.0;
	dmin = 1.0;

        for (int i = 0; i < bindim; i++) {
            for (int j = 0; j < bindim; j++) {
                for (auto &p1 : bins[j + i*bindim]) {
                    particles[p1].ax = particles[p1].ay = 0;
                    for (int ii = i-1; ii <= i+1; ii++) {
                        for (int jj = j-1; jj <= j+1; jj++) {
                            if (ii >= 0 && ii < bindim && jj >= 0 && jj < bindim) {
                                for (auto &p2 : bins[jj + ii*bindim]) {
                                    apply_force(particles[p1], particles[p2], &dmin,&davg,&navg);
                                }
                            }
                        }
                    }
                }
            }
        }

/*
        for( int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
            for (int j = 0; j < n; j++ )
				apply_force( particles[i], particles[j],&dmin,&davg,&navg);
        }
*/

        //
        //  move particles
        //
        std::list <std::tuple <std::list <int>::iterator, std::list <int> *>> toRebin;
        for (int i = 0; i < bindim; i++) {
            for (int j = 0; j < bindim; j++) {
                for (std::list <int>::iterator pIt = bins[j + i*bindim].begin(); pIt != bins[j + i*bindim].end(); pIt++) {
                    int p = *pIt;
                    int xbin = (int) (particles[p].x / bin_size);
                    int ybin = (int) (particles[p].y / bin_size);
                    move(particles[p]);
                    int newxbin = (int) (particles[p].x / bin_size);
                    int newybin = (int) (particles[p].y / bin_size);
                    if (newxbin != xbin || newybin != ybin) {
                        toRebin.push_back(make_pair(pIt, &bins[j + i*bindim]));
                    }
                }
            }
        }

        for (auto &i : toRebin) {
            std::list <int>::iterator pIt = std::get <0> (i);
            std::list <int> *bin = std::get <1> (i);
            int p = *pIt;

            bin->erase(pIt);

            int xbin = (int) (particles[p].x / bin_size);
            int ybin = (int) (particles[p].y / bin_size);

            bins[xbin + ybin*bindim].push_back(p);
        }

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;

          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;

    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    //
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
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
        fprintf(fsum,"%d %g\n",n,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );
    free( particles );
    if( fsave )
        fclose( fsave );

    return 0;
}
