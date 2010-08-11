###############################################################################
# This is an experimental Makefile for building guess both as a command line 
# version and as a parallel PBS version for Simba.
#
# Usage:
#
# make                     - Creates a command line version of guess
# make commandline	   - Same as above
# make pbs		   - Creates a version for parallel execution on Simba
# make clean               - Removes all generated files from the build
#
# The generated programs are placed in the 'output' directory, together with a
# submit script for PBS execution.

# Compiler to use
CC = pgCC

# Compiler options for C++ source files
CXXFLAGS = -w -O2 -Mnobuiltin

# Link flags for the PBS version for Simba
SIMBAPBSLINKFLAGS = -L/usr/lib/mpich/lib/shared -pgf90libs -lmpich

# Where to search for header files
INCLUDES = -I./libraries/include -I./modules -I./framework -I./cru/guessio \
		-I/usr/lib/mpich/include

# Descriptions of targets, paths and rules for building targets
# (nothing to change past this point)
###############################################################################

# Handle conditional compilation
ifeq ($(IO), DEMO)
	CXXFLAGS += -DUSE_DEMO_IO
endif

# Where to place the results of the build
OUTPUT = output
DEPS = $(OUTPUT)/deps
OBJS = $(OUTPUT)/objs

# Where to find source code
SOURCE_PATH := framework libraries/source modules cru/guessio

# Object files used in all versions
COMMONOBJFILES := $(foreach dir, $(SOURCE_PATH), \
	$(patsubst %.cpp, %.o, $(wildcard $(dir)/*.cpp)))

# Same as COMMONOBJFILES but with the correct path in the output directory
DESTOBJFILES := $(patsubst %.o, $(OBJS)/%.o, $(COMMONOBJFILES))

# The object files with main() for each version
COMMANDLINEMAIN = $(OBJS)/command_line_version/main.o
SIMBAPBSMAIN = $(OBJS)/parallel_version/simba/pbs/main.o

MAKESCRIPT = parallel_version/simba/pbs/makescript

SUBMITSCRIPT = $(OUTPUT)/submit.sh

# Determine correct compiler option for dependency generation
# -MM is standard, PGI 6.x however does the whole compilation
# if -MM is used. If -Mcpp=mm is used only the preprocessing
# is done.
ifeq ($(CC),pgCC)
	MM = -Mcpp=mm
else
	MM = -MM
endif

.PHONY : clean commandline pbs

commandline : $(OUTPUT)/guess

$(OUTPUT)/guess : $(DESTOBJFILES) $(COMMANDLINEMAIN) $(OUTPUT)
	$(CC) -o $(OUTPUT)/guess $(CXXFLAGS) $(INCLUDES)  \
		$(filter-out $(OUTPUT), $^)

pbs : $(OUTPUT)/guess_pbs $(SUBMITSCRIPT)

$(OUTPUT)/guess_pbs : $(DESTOBJFILES) $(SIMBAPBSMAIN) $(OUTPUT)
	$(CC) $(SIMBAPBSLINKFLAGS) -o $(OUTPUT)/guess_pbs \
		$(CXXFLAGS) $(INCLUDES) $(filter-out $(OUTPUT), $^)

$(SUBMITSCRIPT) : $(MAKESCRIPT)
	sh $(MAKESCRIPT) $(OUTPUT)/guess_pbs $(SUBMITSCRIPT)

$(OUTPUT) :
	mkdir -p $(OUTPUT)

# Pattern rule for building .o files from .cpp files
$(OBJS)/%.o : %.cpp
	@mkdir -p `dirname $@`
	$(CC) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

# Automatic dependency handling for headers
# (this is what creates the .d files in output/deps)
$(DEPS)/%.d : %.cpp
	@set -e; rm -f $@; \
	mkdir -p $(dir $@); \
        ( $(CC) $(MM) $(CXXFLAGS) $(INCLUDES) $< | tee $@.$$$$ ) &> /dev/null; \
        sed 's,.*\.o[ :]*,$*.o $@ : ,g' < $@.$$$$ > $@; \
        rm -f $@.$$$$

# Include the generated dependency files
-include $(foreach dep, $(COMMONOBJFILES:.o=.d), $(DEPS)/$(dep))

clean :
	rm -rf $(OUTPUT)
