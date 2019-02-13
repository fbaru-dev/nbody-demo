CXX = icpc
COMPFLAGS = -g -std=c++11 -O2 -qopenmp
OPTFLAGS = -xCORE-AVX512 -qopt-zmm-usage=high

REPFLAGS = -qopt-report=5 -qopt-report-filter="GSimulation.cpp,125-175"
INCLUDES =

ADVPRJ = "./adv-ver8"

CXXFLAGS = $(COMPFLAGS) $(OPTFLAGS) $(REPFLAGS)

SOURCES = GSimulation.cpp main.cpp


##########################################
OBJSC = $(SOURCES:.cpp=.o)
##########################################

EXEC=nbody.x

all: clean cpu

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

asm:GSimulation.cpp
	$(info)
	$(info Generate assembly code: )
	$(CXX) $(COMPFLAGS) $(OPTFLAGS) $^ -S

clean:
	rm -f $(OBJSC) nbody.x *.optrpt *.s


#----------------------------------------------------------------
#---------- Intel Advisor Analysis ------------------------------
#----------------------------------------------------------------

survey:
	advixe-cl -collect survey -project-dir $(ADVPRJ) -- ./nbody.x

roofline:
	advixe-cl -collect survey -project-dir $(ADVPRJ) -- ./nbody.x
	advixe-cl -collect tripcounts -flop -project-dir $(ADVPRJ) -- ./nbody.x

map:
	advixe-cl -collect map -mark-up-list=1 -project-dir $(ADVPRJ) -- ./nbody.x

open-gui:
	advixe-gui $(ADVPRJ)/$(ADVPRJ).advixeproj >/dev/null 2>&1 &

clean-results:
	rm -rf $(ADVPRJ)

