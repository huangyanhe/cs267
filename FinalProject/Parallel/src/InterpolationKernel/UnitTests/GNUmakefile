#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.
HOME = ./../../..
IK = $(HOME)/src/InterpolationKernel/
WRITERS = $(HOME)/utils/Writers
VPATH= . $(IK) $(WRITERS)

CXX=clang++
CFLAGS = -g -Wall
CFLAGS += -std=c++11 -I$(IK) -I$(WRITERS)

odir = ./o.$(DIM)d
ddir = ./d.$(DIM)d

LIBS = 

SRCFILES:= $(notdir $(wildcard ./*.cpp $(IK)/*.cpp))

OBJS:=$(patsubst %.cpp,$(odir)/%.o,$(SRCFILES))

DEPS:=$(patsubst $(odir)/%.o,$(ddir)/%.d, $(OBJS))

$(odir)/%.o:%.cpp GNUmakefile
	mkdir -p $(odir);$(CXX) -c $(CFLAGS) $< -o $@
	mkdir -p $(ddir);$(CXX) -MM $(CFLAGS) $< | sed '1s/^/o.1p1d\//' > $*.d;mv $*.d $(ddir)

testIK: GNUmakefile $(OBJS) 
	$(CXX) $(CFLAGS) $(OBJS) $(LIBS) -o interpKernel1p1D.exe

-include $(DEPS)

clean:
	rm -r *.exe $(odir) $(ddir)

listsrc:
	@echo $(SRCFILES)
listobj:
	@echo $(OBJS)
listdep:
	@echo $(DEPS)
