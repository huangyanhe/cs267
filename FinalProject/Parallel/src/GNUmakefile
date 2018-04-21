
HOME = ./..
include $(HOME)/GNUmakefile
DIM=3
CXX = clang++

test1d: FFT1DTest.cpp GNUmakefile $(LIBS_LOCAL)/libfft1D.a
	  $(CXX) $(CFLAGS) FFT1DTest.cpp $(LIB_FLAGS) -o test1d.exe 

testMD: GNUmakefile $(OBJS) FFTMDTest$(osuffix) $(LIBS_LOCAL)/libfft$(DIM)D.a $(LIBS_LOCAL)/libfft1D.a 
	$(CXX) $(CFLAGS) $(OBJS) FFTMDTest$(osuffix) $(LIB_FLAGS) -lfft$(DIM)D -lfft1D -o test$(DIM)D.exe

