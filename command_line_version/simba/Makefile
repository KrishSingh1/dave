##########################################################################
# Make file for LPJ-GUESS, command line version on Simba 
# Written by Ben Smith, 2005-01-26
#
# Updated following software updates to Simba 2007-05-02

# Name of executable file
exe= guess

# Include file locations
includepath=-I/usr/X11R6/include -I/home/ben/c++/lib/include

# GUESS custom libraries
gutilo= /home/ben/c++/lib/gutil/gutil.o
plibo= /home/ben/c++/lib/plib/plib.o

# Compiler options for building object files from C++ source files:
options= -w -O2 -Mnobuiltin $(includepath)

# Compiler options for building executable from object files:
eoptions= -O2 -lX11 -lpthread -Mnobuiltin 

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

$(exe): $(maino) $(guessframeo) $(guessioo) $(drivero) $(canexcho) $(soilwatero) $(growtho) $(vegdynamo) $(somdynamo)
	$(CC) $(eoptions) -o $(exe) $(maino) $(guessframeo) $(guessioo) $(drivero) $(canexcho) $(soilwatero) $(growtho) $(vegdynamo) $(somdynamo) $(gutilo) $(plibo) 

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
