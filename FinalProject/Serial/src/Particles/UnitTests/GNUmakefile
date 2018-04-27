#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.
HOME = ./../../..
IK = $(HOME)/src/InterpolationKernel/
RMDA = $(HOME)/src/RectMDArray/
PARTICLES = $(HOME)/src/Particles/
POISSON = $(HOME)/src/PoissonSolver/
CONV = $(HOME)/src/Hockney
FFT = $(HOME)/src/FFT
WRITERS = $(HOME)/utils/Writers
TIMER = $(HOME)/utils/timer
FFTWDIR = /opt/cray/pe/fftw/3.3.4.6/haswell
VPATH= . $(IK) $(WRITERS) $(RMDA) $(PARTICLES) $(POISSON) $(CONV) $(FFT) $(TIMER)


DIM=1
CXX=CC
CFLAGS = -g -Wall
CFLAGS += -std=c++11 -I$(IK) -I$(WRITERS) -I$(RMDA) -I$(PARTICLES) -I$(POISSON) -I$(CONV) -I$(FFT) -I$(TIMER) -I$(FFTWDIR)/include
CFLAGS += -DDIM=$(DIM)

odir = ./o.$(DIM)d
ddir = ./d.$(DIM)d

LIBS:= -L$(FFTWDIR)/lib -lfftw3 

SRCFILES:= $(notdir $(wildcard ./*.cpp $(IK)/*.cpp $(RMDA)/*.cpp $(PARTICLES)/*.cpp  $(POISSON)/*.cpp $(CONV)/*.cpp $(FFT)/*.cpp $(TIMER)/*.cpp))

OBJS:=$(patsubst %.cpp,$(odir)/%.o,$(SRCFILES))

DEPS:=$(patsubst $(odir)/%.o,$(ddir)/%.d, $(OBJS))

$(odir)/%.o:%.cpp GNUmakefile
	mkdir -p $(odir);$(CXX) -c $(CFLAGS) $< -o $@
	mkdir -p $(ddir);$(CXX) -MM $(CFLAGS) $< | sed '1s/^/o.1p1d\//' > $*.d;mv $*.d $(ddir)

testDeposition: GNUmakefile $(OBJS) 
	$(CXX) $(CFLAGS) $(OBJS) $(LIBS) -o deposition.exe

-include $(DEPS)

clean:
	rm -r *.exe $(odir) $(ddir)

listsrc:
	@echo $(SRCFILES)
listobj:
	@echo $(OBJS)
listdep:
	@echo $(DEPS)
