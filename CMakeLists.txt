################################################################################
# CMake configuration file for building LPJ-GUESS
#
# To build LPJ-GUESS with this build system cmake needs to be installed.
# If it's not installed it can be downloaded for free from www.cmake.org.
#

cmake_minimum_required(VERSION 2.8)
project(guess)

# Compiler flags for building on UNIX
if (UNIX)
  list(APPEND CMAKE_CXX_FLAGS "-Mnobuiltin -w")
endif()

# The IO variable decides which IO module to use. Rather than editing
# it here, you should pass it as a parameter to cmake, or edit it
# in the graphical cmake tool (cmake-gui on Windows, ccmake on Unix).
set(IO "cru" CACHE STRING "The IO module to use (cru/demo)")

# Where the compiler should search for header files
include_directories(${guess_SOURCE_DIR}/framework ${guess_SOURCE_DIR}/libraries/gutil ${guess_SOURCE_DIR}/libraries/plib ${guess_SOURCE_DIR}/modules)

# The following directories contain source code and 
# additional CMakeLists.txt files
add_subdirectory(framework)
add_subdirectory(modules)
add_subdirectory(cru)
add_subdirectory(libraries)

# Add the command line program's target
if (WIN32)
  # Let the exe be called guesscmd so it doesn't collide with the dll target
  set(guess_command_name "guesscmd")
else()
  # On Unix we don't have the dll so the command line binary can be called guess
  set(guess_command_name "guess")
endif()

add_executable(${guess_command_name} ${guess_sources} ${gutil_sources} ${plib_sources} command_line_version/main.cpp)

if (WIN32)
  # Create guess.dll (used with the graphical Windows shell)
  add_library(guess SHARED ${guess_sources} ${gutil_sources} ${plib_sources} windows_version/dllmain.cpp)
endif (WIN32)

# The custom build rule for generating the PBS submit script from an m4 template
# The submit script is generated each time the command line binary is built.
if (UNIX)
  # Choose template file depending on system
  if (DEFINED SYSTEM AND SYSTEM STREQUAL "milleotto")
    set(submit_template ${guess_SOURCE_DIR}/parallel_version/lunarc/milleotto/submit.m4)
  else()
    set(submit_template ${guess_SOURCE_DIR}/parallel_version/simba/pbs/submit.m4)
  endif()
  
  add_custom_command(TARGET ${guess_command_name}
                     POST_BUILD
		     COMMAND m4 -DBINARY=${CMAKE_CURRENT_BINARY_DIR}/${guess_command_name} -DDATE="`date`" ${submit_template} > ${CMAKE_BINARY_DIR}/submit.sh
		     COMMENT "Generating submit script")
endif (UNIX)

if (IO STREQUAL "demo")
  # Define the USE_DEMO_IO macro when compiling all sources
  set_property(
    DIRECTORY
    PROPERTY COMPILE_DEFINITIONS USE_DEMO_IO
    )
endif(IO STREQUAL "demo")

if (UNIX)
   # pgCC 6 doesn't seem to recognize -rdynamic, so remove it
   # (we shouldn't need it anyway)
   list(REMOVE_ITEM CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-rdynamic")
endif(UNIX)
