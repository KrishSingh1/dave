##########################################################################
# Make file for LPJ-GUESS (Unix command line version)
# Written by Ben Smith, 2001-11-28
#
# Things you may have to change for your particular system:
#   1. Path to header (include) files for custom libraries (includepath)
#   2. Paths to library binaries (libraries) - you may have to create
#      these separately

# Name of executable file
exe= guess

# Include file locations
includepath= -I/home/ben/c++/lib/include

# Libraries to be linked to build executable (full path names)
libraries= /home/ben/c++/lib/gutil/gutil.o /home/ben/c++/lib/plib/plib.o

# Compiler options for building object files from C++ source files:
options= -w -O3 $(includepath)

# Compiler options for building executable from object files:
eoptions= -O3 $(includepath)

# Compiler
CC= pgCC

# Names of object, source and header files
# (nothing to change past this point)
##########################################################################

maino = main.o
mainc = main.cpp

guessframeo= guess.o
guessframec= guess.cpp
guessframeh= guess.h

guessioo= guessio_cru.o
guessioc= guessio_cru.cpp
guessioh= guessio.h

drivero= driver.o
driverc= driver.cpp
driverh= driver.h

canexcho= canexch.o
canexchc= canexch.cpp
canexchh= canexch.h

soilwatero= soilwater.o
soilwaterc= soilwater.cpp
soilwaterh= soilwater.h

growtho= growth.o
growthc= growth.cpp
growthh= growth.h

vegdynamo= vegdynam.o
vegdynamc= vegdynam.cpp
vegdynamh= vegdynam.h

somdynamo= somdynam.o
somdynamc= somdynam.cpp
somdynamh= somdynam.h

bvoco= bvoc.o
bvocc= bvoc.cpp
bvoch= bvoc.h

# List of target : dependency directives
# NB: executable must be first target

$(exe): $(maino) $(guessframeo) $(guessioo) $(drivero) $(canexcho) $(soilwatero) $(growtho) $(vegdynamo) $(somdynamo) $(bvoco)
	$(CC) $(eoptions) -o $(exe) $(maino) $(guessframeo) $(guessioo) $(drivero) $(canexcho) $(soilwatero) $(growtho) $(vegdynamo) $(somdynamo) $(bvoco) $(libraries) 

$(maino): $(mainc) $(guessframeh)
	$(CC) $(options) -c $(mainc)

$(guessframeo): $(guessframec) $(guessframeh) $(driverh) $(guessioh) $(canexchh) $(soilwaterh) $(somdynamh) $(growthh) $(vegdynamh)
	$(CC) $(options) -c $(guessframec)

$(guessioo): $(guessioc) $(guessframeh) $(guessioh) $(driverh)
	$(CC) $(options) -c $(guessioc)

$(drivero): $(driverc) $(guessframeh) $(driverh)
	$(CC) $(options) -c $(driverc)

$(canexcho): $(canexchc) $(guessframeh) $(canexchh)
	$(CC) $(options) -c $(canexchc)

$(soilwatero): $(soilwaterc) $(guessframeh) $(soilwaterh)
	$(CC) $(options) -c $(soilwaterc)

$(growtho): $(growthc) $(guessframeh) $(growthh)
	$(CC) $(options) -c $(growthc)

$(vegdynamo): $(vegdynamc) $(guessframeh) $(vegdynamh) $(growthh)
	$(CC) $(options) -c $(vegdynamc)

$(somdynamo): $(somdynamc) $(guessframeh) $(somdynamh)
	$(CC) $(options) -c $(somdynamc)

$(bvoco): $(bvocc) $(bvoch)
	$(CC) $(options) -c $(bvocc)
