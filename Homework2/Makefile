#
# Edison - NERSC 
#
# Intel Compilers are loaded by default; for other compilers please check the module list
#
CC = CC
#CC = clang++
MPCC = CC
<<<<<<< HEAD:Homework2/Makefile
OPENMP = -openmp #Note: this is the flag for Intel compilers. Change this to -fopenmp for GNU compilers. See http://www.nersc.gov/users/computational-systems/edison/programming/using-openmp/
CFLAGS = -O3 -std=c++11
#CFLAGS = -g -std=c++11
=======
OPENMP = -fopenmp #Note: this is the flag for Intel compilers. Change this to -fopenmp for GNU compilers. See http://www.nersc.gov/users/computational-systems/edison/programming/using-openmp/
CFLAGS = -O3
>>>>>>> c_version:Makefile
LIBS =


TARGETS = serial openmp mpi autograder

all:	$(TARGETS)

serial: serial.o common.o decomposition.o
	$(CC) -o $@ $(LIBS) serial.o common.o decomposition.o
autograder: autograder.o common.o
	$(CC) -o $@ $(LIBS) autograder.o common.o
<<<<<<< HEAD:Homework2/Makefile
openmp: openmp.o common.o
	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o common.o
mpi: mpi.o common.o Particle.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common.o
=======
openmp: openmp.o common.o decomposition.o
	$(CC) -o $@ $(LIBS) $(OPENMP) openmp.o common.o decomposition.o
mpi: mpi.o common.o decomposition_mpi.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common.o decomposition_mpi.o
>>>>>>> c_version:Makefile

autograder.o: autograder.cpp common.h
	$(CC) -c $(CFLAGS) autograder.cpp
openmp.o: openmp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp -o openmp.o
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp
mpi.o: mpi.cpp common.h Particle.h
	$(MPCC) -c $(CFLAGS) mpi.cpp
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp
decomposition.o: decomposition.cpp decomposition.h
	$(CC) -c $(CFLAGS) decomposition.cpp
decomposition_mpi.o: decomposition_mpi.cpp decomposition_mpi.h
	$(CC) -c $(CFLAGS) decomposition_mpi.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt *.error
