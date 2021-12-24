#Directory for sourcecode:
IDIR =./sourcecode
ODIR=$(IDIR)

OMP=-fopenmp
OPTIMIZE=-Ofast #-01-3, -Ofast, or blank

WARNINGS=-Wall #-Wextra -Wpedantic

CXX=g++

CXXFLAGS=-I$(IDIR) -std=c++11 $(WARNINGS) $(OMP) $(OPTIMIZE) -llapacke -llapack -lblas

detected_OS := $(shell uname -s) #will return the Operating system name
$(info )
$(info Detected operating system: $(detected_OS))
$(info )

# ### Following for Andrei Max OS:
MACOSX=Darwin #NB: needs a trailing space!
ifeq ($(detected_OS),$(MACOSX))
  $(info Trying to use MAC OSX flags Andrei)
  $(info )
  CXX=/usr/local/opt/llvm/bin/clang++
  LDFLAGS += -L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib
  CXXFLAGS += -I/usr/local/opt/llvm/include -I/usr/local/opt/llvm/include/c++/v1/
endif
LIBS=-lgsl -lgslcblas -lm

################################################################################

#default tagret
all: convertmout.x dateConvert.x fetchAndCheckJPL30s.x impTwWindow.x \
	impTwWindowPattern.x plotClockData.x processNoise.x testMethods.x \
	stationNoise.x skewness.x excessPower.x
	$(info  Success!)

# #All files depend on these:
$(ODIR)/miscFunctions.o: $(IDIR)/miscFunctions.cpp $(IDIR)/miscFunctions.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

# #All files depend on these:
$(ODIR)/mathematicsFunctions.o: $(IDIR)/mathematicsFunctions.cpp \
$(IDIR)/mathematicsFunctions.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

# $(IDIR)/readJplBiasData_1s.o: $(IDIR)/readJplBiasData_1s.cpp
# 	$(CXX) -c -o $@ $< $(CXXFLAGS)

#All files depend on these:
$(ODIR)/JplGpsDataClass.o: $(IDIR)/JplGpsDataClass.cpp \
$(IDIR)/JplGpsDataClass.h \
$(IDIR)/miscFunctions.o $(IDIR)/mathematicsFunctions.o
	$(CXX) -c -o $@ $< $(CXXFLAGS)

# #All files depend on these:
OBJ = $(ODIR)/JplGpsDataClass.o $(IDIR)/miscFunctions.o $(IDIR)/mathematicsFunctions.o

#Dependencies/tagrets for each other program:

$(ODIR)/plotClockData.o: $(IDIR)/plotClockData.cpp $(IDIR)/plotClockData.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/convertmout.o: $(IDIR)/convertmout.cpp $(IDIR)/convertmout.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/noiseProfileFunctions.o: $(IDIR)/noiseProfileFunctions.cpp \
$(IDIR)/noiseProfileFunctions.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)


all: convertmout.x dateConvert.x fetchAndCheckJPL30s.x impTwWindow.x \
impTwWindowPattern.x plotClockData.x processNoise.x testMethods.x search.x

$(ODIR)/methods.o: $(IDIR)/methods.cpp $(IDIR)/methods.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/NumericCdfInverseClass.o: $(IDIR)/NumericCdfInverseClass.cpp \
$(IDIR)/NumericCdfInverseClass.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/SearchClass.o: $(IDIR)/SearchClass.cpp $(IDIR)/SearchClass.h 
	$(CXX) -c -o $@ $< $(CXXFLAGS)

## Targets:

COMP = $(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS) $(LDFLAGS)

convertmout.x: $(OBJ) $(IDIR)/convertmout.o
	$(COMP)

dateConvert.x: $(OBJ) $(IDIR)/dateConvert.o
	$(COMP)

fetchAndCheckJPL30s.x: $(OBJ) $(IDIR)/fetchAndCheckJPL30s.o
	$(COMP)

impTwWindow.x: $(OBJ) $(IDIR)/impTwWindow.o
	$(COMP)

impTwWindowPattern.x: $(OBJ) $(IDIR)/impTwWindowPattern.o
	$(COMP)

plotClockData.x: $(OBJ) $(IDIR)/plotClockData.o
	$(COMP)

processNoise.x: $(OBJ) $(IDIR)/processNoise.o $(ODIR)/noiseProfileFunctions.o
	$(COMP)


testMethods.x: $(OBJ) $(IDIR)/testMethods.o $(IDIR)/methods.o $(IDIR)/NumericCdfInverseClass.o
	$(COMP)

search.x: $(OBJ) $(IDIR)/search.o $(IDIR)/methods.o $(IDIR)/SearchClass.o $(IDIR)/NumericCdfInverseClass.o
	$(COMP)

stationNoise.x: $(OBJ) $(IDIR)/stationNoise.o $(ODIR)/noiseProfileFunctions.o
	$(COMP)

skewness.x: $(OBJ) $(IDIR)/skewness.o
	$(COMP)

excessPower.x: $(OBJ) $(IDIR)/excessPower.o
	$(COMP)


.PHONY: clean
clean:
	rm -f *.x *~ $(IDIR)/*.o $(IDIR)/*~
