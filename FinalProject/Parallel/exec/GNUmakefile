#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.
HOME = ./..
WRITERS = $(HOME)/utils/Writers
RMDA = $(HOME)/src/RectMDArray
FFT = $(HOME)/src/FFT
PS = $(HOME)/src/PoissonSolver
PARTICLES = $(HOME)/src/Particles
IK = $(HOME)/src/InterpolationKernel/
TIMER = $(HOME)/utils/timer
VPATH= . $(HOME) $(PARTICLES) $(RMDA) $(FFT) $(PS) $(TIMER) $(WRITERS) $(IK)
ANIMATION = 1
OPENMP = -qopenmp

DIM= 1
CXX = CC
MPCC = CC
FFTWDIR = /opt/cray/pe/fftw/3.3.4.6/haswell
#FFTWDIR = /usr/local/Cellar/fftw/3.3.6-pl2
CFLAGS = -g -Wall
CFLAGS += $(MPILIBS)
#CFLAGS = -O3
CFLAGS += -std=c++11 -I$(PARTICLES) -I$(RMDA) -I$(IK) -I$(PS) -I$(FFT) -I$(TIMER) -I$(WRITERS) -I$(FFTWDIR)/include

odir = ./o.$(DIM)d
ddir = ./d.$(DIM)d

CFLAGS += -DDIM=$(DIM) -D ANIMATION=$(ANIMATION)
LIBS:=  -L$(FFTWDIR)/lib -lfftw3

CFLAGS += $(OPENMP)

SRCFILES:= $(notdir $(wildcard $(TIMER)/*.cpp $(RMDA)/*.cpp $(IK)/*.cpp $(PS)/*.cpp $(WRITERS)/*.cpp $(FFT)/*.cpp $(PARTICLES)/*.cpp ./*.cpp))

OBJS:=$(patsubst %.cpp,$(odir)/%.o,$(SRCFILES))


DEPS:=$(patsubst $(odir)/%.o,$(ddir)/%.d, $(OBJS))

$(odir)/%.o:%.cpp GNUmakefile
	mkdir -p $(odir);$(CXX) -c $(CFLAGS) $< -o $@
#	mkdir -p $(ddir);$(CXX) -MM $(CFLAGS) $<  > $*.d;mv $*.d $(ddir)
	mkdir -p $(ddir);$(CXX) -MM $(CFLAGS) $< | sed '1s/^/o.$(DIM)d\//' > $*.d;mv $*.d $(ddir)

export OMP_NUM_THREADS=32
PICLLD: GNUmakefile $(OBJS) 
	$(CXX) $(CFLAGS) $(OBJS) $(LIBS) -o PIC$(DIM)p$(DIM)D.exe

-include $(DEPS)

clean:
	rm -r *.exe $(odir) $(ddir)

listsrc:
	@echo $(SRCFILES)
listobj:
	@echo $(OBJS)
listdep:
	@echo $(DEPS)

