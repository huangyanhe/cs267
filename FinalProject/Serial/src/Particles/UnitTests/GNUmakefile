#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.
HOME = ./../../..
IK = $(HOME)/src/InterpolationKernel/
RMDA = $(HOME)/src/RectMDArray/
PARTICLES = $(HOME)/src/Particles/
PS = $(HOME)/src/PoissonSolver
FFT = $(HOME)/src/FFT
WRITERS = $(HOME)/utils/Writers
TIMER = $(HOME)/utils/timer
FFTWDIR = /usr/local/Cellar/fftw/3.3.6-pl2
VPATH= . $(IK) $(WRITERS) $(RMDA) $(PARTICLES) $(FFT) $(TIMER) $(PS)


DIM=1
CXX=clang++
CFLAGS = -g -Wall
CFLAGS += -std=c++11 -I$(IK) -I$(WRITERS) -I$(RMDA) -I$(PARTICLES) -I$(FFT) -I$(PS) -I$(TIMER) -I$(FFTWDIR)/include
CFLAGS += -DDIM=$(DIM)

odir = ./o.$(DIM)d
ddir = ./d.$(DIM)d

LIBS:= -L$(FFTWDIR)/lib -lfftw3 

SRCFILES:= $(notdir $(wildcard ./*.cpp $(IK)/*.cpp $(RMDA)/*.cpp $(PS)/*.cpp $(PARTICLES)/*.cpp $(FFT)/*.cpp $(TIMER)/*.cpp))

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
