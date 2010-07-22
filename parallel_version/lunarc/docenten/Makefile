##########################################################################
# Make file to build and set up LPJ-GUESS for parallel runs on
# docenten.lunarc.lu.se
#
# Written by Ben Smith, Lund University, 2005-01-21

# Name of executable file
exe= guess

# Path to include files for gutil and plib
include= -I/home/ben/c++/lib/include

# Paths to object files for gutil and plib
gutilo= /home/ben/c++/lib/gutil/bin-docenten/gutil.o
plibo= /home/ben/c++/lib/plib/bin-docenten/plib.o

# Compiler options for building object files from C++ source files:
options= -w -O3 -Mnobuiltin $(include)

# Compiler options for building executable from object files:
eoptions = -O3 -Mnobuiltin

# Compiler
CC= mpiCC 

# Names of object, source, header and script files
# (nothing to change past this point)
##########################################################################

makescript= makescript
submitscript= submit.sh

appende= append
appendc= append.cpp

splitgride= splitgrid
splitgridc= splitgrid.cpp

maino= main.o
mainc= main.cpp

guessframeo= guess.o
guessframec= guess.cpp
guessframeh= guess.h

guessioo= guessio.o
guessioc= guessio.cpp
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



# List of target : dependency directives
# NB: executable must be first target

$(submitscript): $(exe) $(appende) $(splitgride) $(makescript)
	sh makescript $(exe) $(submitscript)

$(appende): $(appendc)
	$(CC) $(eoptions) $(include) $(gutilo) -o $(appende) $(appendc)

$(splitgride): $(splitgridc)
	$(CC) $(eoptions) $(include) $(gutilo) -o $(splitgride) $(splitgridc)

$(exe): $(maino) $(guessframeo) $(guessioo) $(drivero) $(canexcho) $(soilwatero) $(growtho) $(vegdynamo) $(somdynamo)
	$(CC) $(eoptions) $(gutilo) $(plibo) -o $(exe) $(maino) $(guessframeo) $(guessioo) $(drivero) $(canexcho) $(soilwatero) $(growtho) $(vegdynamo) $(somdynamo) 

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
