######################################################################
# Make file for GUTIL and PLIB libraries
# Written by Ben Smith, 2001-11-28

# Name of Gutil archive file
gutila= gutil.a

# Name of Plib archive file
pliba= plib.a

# Include file location
includepath= ./include

# Compiler options
options= -w -fast -I$(includepath) -xar

# Nothing to change past this point
######################################################################

gutilc= gutil.cpp
gutilh= gutil.h

plibc= plib.cpp
plibh= plib.h

$(pliba): $(plibc) $(includepath)/$(plibh) $(gutila)
	CC $(options) -o $(pliba) $(plibc)

$(gutila): $(gutilc) $(includepath)/$(gutilh)
	CC $(options) -o $(gutila) $(gutilc)
