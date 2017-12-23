# define the C++ compiler to use
CC = g++

# define any compile-time flags
CFLAGS = -Wall -g -std=c++11

# define any directories containing header files other than /usr/include
INCLUDES = -I. -Iclusters/ -I../dlib-18.8/ -I/home/vultor/novozerg/dlib-18.18/

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS = -lpthread -std=c++11 -lX11

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
LIBS = 

# define the C++ source files
SRCS = Random.cpp clusters/Similarity.cpp clusters/Hungarian.cpp clusters/MarquesEnantiomers.cpp AdministrateCreation.cpp AuxMath.cpp AuxMathGa.cpp BasicOperators.cpp Creation.cpp FuzzyAdministration.cpp GAPrinting.cpp GeneticAlgorithm.cpp Hooklaw.cpp main.cpp ParallelOptimization.cpp Predator.cpp UserFitness.cpp UserOperators.cpp clusters/ClustersFitness.cpp clusters/ClustersOperators.cpp clusters/Fitness.cpp clusters/InitializeAtoms.cpp clusters/ReadGaInput.cpp clusters/ReadQuantumOutput.cpp clusters/WriteQuantumInput.cpp clusters/Experiment.cpp clusters/Derivative.cpp Printing.cpp ../dlib-18.18/dlib/all/source.cpp


# define the C++ object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.cpp=.o)

# define the executable file 
MAIN = zerg.exe

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:    $(MAIN)
	@echo  zerg.exe foi compilado com sucesso:

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)
	$(RM) clusters/*.o

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
