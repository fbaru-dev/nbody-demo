CXX = icpc
COMPFLAGS =  -g -std=c++11 -O2
 INCLUDES = 

CXXFLAGS = $(COMPFLAGS) 

SOURCES = GSimulation.cpp main.cpp

ADVPRJ = "./adv-ver0"

ifeq ($(REPORT), yes)
	CXXFLAGS+=-qopt-report=5
ifeq ($(FILTER), yes) 
	CXXFLAGS+=-qopt-report-phase=vec -qopt-report-filter="GSimulation.cpp,125-175"
endif
endif

.SUFFIXES: .o .cpp

##########################################
OBJSC = $(SOURCES:.cpp=.o)
##########################################

EXEC=nbody.x

all: cpu

%.o: %.cpp
	$(info )
	$(info Compiling the object file for CPU: )
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@ 

cpu: $(OBJSC)
	$(info )
	$(info Linking the CPU executable:)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(EXEC) $(OBJSC)
	
run: 
	$(info )
	$(info Run the default test case on CPU: )
	./nbody.x 
	
clean: 
	rm -f $(OBJSC) nbody.x *.optrpt

#----------------------------------------------------------------
#---------- Intel Advisor Analysis ------------------------------
#----------------------------------------------------------------

survey:
	advixe-cl -collect survey -project-dir $(ADVPRJ) -- ./nbody.x

roofline:
	advixe-cl -collect survey -project-dir $(ADVPRJ) -- ./nbody.x
	advixe-cl -collect tripcounts -flop -project-dir $(ADVPRJ) -- ./nbody.x

open-gui:
	advixe-gui $(ADVPRJ)/$(ADVPRJ).advixeproj >/dev/null 2>&1 &

clean-results:
	rm -rf $(ADVPRJ)
 


