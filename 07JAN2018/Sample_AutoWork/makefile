##// Compiler that supports C++11
CXX=g++ -std=c++0x -std=gnu++0x

##// Mycroft C++11 compiler
##CXX=/export/apps/gcc-4.9.1/bin/g++ -Wl,-rpath,/export/apps/gcc-4.9.1/lib64 -std=gnu++11

CXXFLAGS=-c
LDFLAGS=
SOURCES=./alglib/alglibinternal.cpp ./alglib/alglibmisc.cpp ./alglib/ap.cpp ./alglib/dataanalysis.cpp ./alglib/diffequations.cpp ./alglib/fasttransforms.cpp ./alglib/integration.cpp ./alglib/interpolation.cpp ./alglib/linalg.cpp ./alglib/optimization.cpp ./alglib/solvers.cpp ./alglib/specialfunctions.cpp ./alglib/statistics.cpp autowork.cpp varclass.cpp userinterface.cpp functions.cpp sysinfo.cpp structureclass.cpp moltemplatelmpdata.cpp lmpscripts.cpp workscripts.cpp amdatanalysis.cpp alglibfittingkernel.cpp fitdata.cpp theorytest.cpp fitdielectric.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=AUTOWORK

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -rf $(OBJECTS)
	rm -rf $(EXECUTABLE)
