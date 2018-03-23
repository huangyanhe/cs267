#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.
HOME = ./..
WRITERS = $(HOME)/utils/Writers
RMDA = $(HOME)/src/RectMDArray
FFT = $(HOME)/src/fftTools
CONV = $(HOME)/src/Hockney
PARTICLES = $(HOME)/src/Particles
TIMER = $(HOME)/utils/timer
VPATH= . $(HOME) $(PARTICLES) $(RMDA) $(FFT) $(CONV) $(TIMER) $(WRITERS)
ANIMATION = 1

DIM=2
#CXX=g++
CXX=clang++
FFTWDIR = /usr/local
#CFLAGS = -g -Wall
CFLAGS = -O3
CFLAGS += -std=c++11 -I$(PARTICLES) -I$(RMDA) -I$(FFT) -I$(CONV) -I$(TIMER) -I$(WRITERS) -I$(FFTWDIR)/include

odir = ./o.$(DIM)d
ddir = ./d.$(DIM)d

CFLAGS += -DDIM=$(DIM) -D ANIMATION=$(ANIMATION)
LIBS:=  -L$(FFTWDIR)/lib -lfftw3 

SRCFILES:= $(notdir $(wildcard $(TIMER)/*.cpp $(RMDA)/*.cpp $(WRITERS)/*.cpp $(FFT)/*.cpp ./*.cpp $(CONV)/*.cpp $(PARTICLES)/*.cpp))

OBJS:=$(patsubst %.cpp,$(odir)/%.o,$(SRCFILES))


DEPS:=$(patsubst $(odir)/%.o,$(ddir)/%.d, $(OBJS))

$(odir)/%.o:%.cpp GNUmakefile
	mkdir -p $(odir);$(CXX) -c $(CFLAGS) $< -o $@
#	mkdir -p $(ddir);$(CXX) -MM $(CFLAGS) $<  > $*.d;mv $*.d $(ddir)
	mkdir -p $(ddir);$(CXX) -MM $(CFLAGS) $< | sed '1s/^/o.$(DIM)d\//' > $*.d;mv $*.d $(ddir)

vortex2D: GNUmakefile $(OBJS) 
	$(CXX) $(CFLAGS) $(OBJS) $(LIBS) -o vortex$(DIM)D.exe

-include $(DEPS)

clean:
	rm -r *.exe $(odir) $(ddir)

listsrc:
	@echo $(SRCFILES)
listobj:
	@echo $(OBJS)
listdep:
	@echo $(DEPS)

