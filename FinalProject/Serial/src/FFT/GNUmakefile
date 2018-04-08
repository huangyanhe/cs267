#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.

HOME = ./..
include $(HOME)/GNUmakefile2

CFLAGS += -DDIM=$(DIM) -g
# Archiving tool on MACs.
ifeq ($(system),Darwin)
  LIBTOOL = libtool -static -s -o
endif

1DFFTSRC := $(wildcard $(FFT_HOME)/FFT1D*.cpp)
1DFFTSRC += $(FFT_HOME)/PowerItoI.cpp 
#$(FFT_HOME)/FFTRecursive.cpp
1DFFTOBJS = $(1DFFTSRC:.cpp=.o)

ifneq ($(DIM),1)
MDFFTOBJS = FFTMD$(osuffix)
libfft$(DIM)D.a: GNUmakefile $(MDFFTOBJS)
	$(LIBTOOL) libfft$(DIM)D.a $(MDFFTOBJS) 
	mkdir -p ../lib;mv libfft$(DIM)D.a $(LIBS_LOCAL);rm *.o
endif

libfft1D:GNUmakefile $(wildcard $(FFT_HOME)/*.cpp) $(wildcard $(FFT_HOME)/*.H)
	make $(1DFFTOBJS)
	$(LIBTOOL) libfft1D.a $(1DFFTOBJS)
	mkdir -p ../lib;mv libfft1D.a $(LIBS_LOCAL);rm *.o

-include $(1DFFTOBJS:.o=.d)
-include $(MDFFTOBJS:.o=.d)
