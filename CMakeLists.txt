################################################################################
# CMake configuration file for building LPJ-GUESS
#
# To build LPJ-GUESS with this build system cmake needs to be installed.
# If it's not installed it can be downloaded for free from www.cmake.org.
#

cmake_minimum_required(VERSION 3.10...4.4.4)

project(guess)

# C++ standard 17 for std::filesystem and other features.
set(CMAKE_CXX_STANDARD 17)

# Output compression is only possible if zlib is available.
set(COMPRESS_OUTPUT 0 BOOL "Compress output or not")
find_package(ZLIB QUIET)
if(ZLIB_FOUND)
  set(COMPRESS_OUTPUT 1)
  add_compile_definitions(COMPRESS_OUTPUT)
  set(LIBS ${LIBS} z)

  message(STATUS "Output compression has been enabled")
else()
  message(STATUS "find_package: not found: ZLIB. See LPJ-GUESS documentation." )
  message(STATUS "Output compression disabled." )
endif()

# Compiler flags for building with Microsoft Visual C++
if (MSVC)
  # Disable warnings about using secure functions like sprintf_s instead of
  # regular sprintf etc.
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D _CRT_SECURE_NO_WARNINGS /D _SCL_SECURE_NO_WARNINGS")
else()
  # Add a temporary flag for Linux and macOS
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
endif()

# The following are configuration variables for the build. 
# Rather than editing them here, you should either edit them
# in the graphical cmake tools (cmake-gui on Windows, ccmake
# on Unix) or pass them in as parameters to cmake.

# A variable controlling whether or not to include unit tests
# Unit tests are disabled in old VS 6.0, since CATCH doesn't compile in such 
# an old compiler.
if (NOT MSVC_VERSION EQUAL "1200")
  set(UNIT_TESTS "OFF" CACHE BOOL "Whether to include unit tests")
endif()

if (UNIX)
  # Setup the SYSTEM variable, currently only used to choose which 
  # submit.sh to generate (for submitting to job queue)

  # Figure out what value it should have initially, based on the
  # environment variable ARCH if it's set.
  if (NOT $ENV{ARCH} STREQUAL "")
	set(DEFAULT_SYSTEM $ENV{ARCH})
  else()
	set(DEFAULT_SYSTEM "")
  endif()

  set(SYSTEM ${DEFAULT_SYSTEM} CACHE STRING "System to build for (empty (=simba), gimle, lunarc, levante or multicore)")
  
  # A temporary helper for Lund users lingering with $ARCH=="aurora" instead of "lunarc".
  if (SYSTEM STREQUAL "aurora")
    set(SYSTEM "lunarc")
  endif()
  
endif (UNIX)

# Where to search for cmake modules
# (used by cmake when the include() command is used in a cmake file)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Detect the platform-specific stack trace API.
include(CheckIncludeFileCXX)
include(CheckCXXCompilerFlag)
include(CheckCXXSourceCompiles)
check_include_file_cxx("execinfo.h" HAVE_EXECINFO_H)
if (HAVE_EXECINFO_H)
  check_cxx_source_compiles("#include <execinfo.h>\nint main() { void *frames[1]; return backtrace(frames, 1) < 0; }" HAVE_EXECINFO_BACKTRACE)

  if (NOT HAVE_EXECINFO_BACKTRACE)
    find_library(EXECINFO_LIBRARY NAMES execinfo)
    if (EXECINFO_LIBRARY)
      set(CMAKE_REQUIRED_LIBRARIES ${EXECINFO_LIBRARY})
	  unset(HAVE_EXECINFO_BACKTRACE CACHE)
      check_cxx_source_compiles("#include <execinfo.h>\nint main() { void *frames[1]; return backtrace(frames, 1) < 0; }" HAVE_EXECINFO_BACKTRACE)
      unset(CMAKE_REQUIRED_LIBRARIES)
    endif()
  endif()
endif()

if (HAVE_EXECINFO_BACKTRACE)
  add_definitions(-DHAVE_EXECINFO_BACKTRACE)
  if (EXECINFO_LIBRARY)
    set(LIBS ${LIBS} ${EXECINFO_LIBRARY})
  endif()

  if (UNIX AND NOT APPLE)
    check_cxx_compiler_flag("-rdynamic" HAVE_RDYNAMIC)
    if (HAVE_RDYNAMIC)
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -rdynamic")
    endif()
  endif()
  message(STATUS "Stack traces enabled with execinfo")
elseif (NOT WIN32)
  message(STATUS "Stack traces unavailable: execinfo was not found")
endif()

if (WIN32)
  find_library(DBGHELP_LIBRARY NAMES dbghelp)
  if (DBGHELP_LIBRARY)
    add_definitions(-DHAVE_DBGHELP)
    set(LIBS ${LIBS} ${DBGHELP_LIBRARY})
    message(STATUS "Windows stack traces enabled with symbol names")
  else()
    message(STATUS "Windows stack traces enabled with addresses only")
  endif()
endif()

# NetCDF - used if found.
find_package(NetCDF QUIET)
if (NETCDF_FOUND)
  include_directories(${NETCDF_INCLUDE_DIRS})
  set(LIBS ${LIBS} ${NETCDF_LIBRARIES})
  add_definitions(-DHAVE_NETCDF)
else()
  message("find_package: not found: NETCDF_C_LIBRARY and NETCDF_INCLUDE_DIR. See LPJ-GUESS documentation." )
  message("NETCDF disabled." )
endif()

# MPI - used if found (not needed on Windows).
if (NOT CMAKE_HOST_WIN32)
  find_package(MPI QUIET)
endif()

# These are deprecated according to documentation in the FindMPI module,
# but for some reason not hidden. Let's not display them for the typical
# LPJ-GUESS user who doesn't have MPI installed.
mark_as_advanced(MPI_LIBRARY MPI_EXTRA_LIBRARY)

if (MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  set(LIBS ${LIBS} ${MPI_LIBRARIES})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS}")
  add_definitions(-DHAVE_MPI)
endif()

# MSVC doesn't seem to support variable-length arrays, so let's disable them on
# all compilers. Note that the MSVC compiler doesn't support -Werror=vla.
if (!MSVC)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=vla")

  # Let's also ignore deprecated declations warnings; without this we see a lot
  # of noise about the usage of std::auto_ptr, which is something that's being
  # worked on but not something we need to be informed of repeatedly.
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
endif()

# Use CCache if it's available.
find_program(CCACHE_FOUND ccache)
if (CCACHE_FOUND)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
endif(CCACHE_FOUND)

# AddressSanitizer (ASan) option (opt-in)
option(ENABLE_ASAN "Enable AddressSanitizer" OFF)
if (ENABLE_ASAN)
  if (UNIX AND NOT MSVC)
    # Only enable for GCC/Clang toolchains
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
      # Use conservative optimization level for better stack traces
      set(ASAN_CFLAGS "-fsanitize=address -fno-omit-frame-pointer")
      set(ASAN_LDFLAGS "-fsanitize=address")

      set(CMAKE_C_FLAGS           "${CMAKE_C_FLAGS} ${ASAN_CFLAGS}")
      set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} ${ASAN_CFLAGS}")
      set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${ASAN_LDFLAGS}")
      set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ASAN_LDFLAGS}")

      message(STATUS "AddressSanitizer enabled")
    else()
      message(WARNING "ENABLE_ASAN requested, but compiler ${CMAKE_CXX_COMPILER_ID} is not supported for ASan")
    endif()
  else()
    message(WARNING "ENABLE_ASAN requested, but platform is not supported for ASan")
  endif()
endif()

# Setup the include paths.
include_directories(
  ${guess_SOURCE_DIR}/framework
  ${guess_SOURCE_DIR}/libraries/gutil
  ${guess_SOURCE_DIR}/libraries/plib
  ${guess_SOURCE_DIR}/libraries/guessnc
  ${guess_SOURCE_DIR}/modules
  ${guess_SOURCE_DIR}/cru/guessio)

# The following directories contain source code and 
# additional CMakeLists.txt files
add_subdirectory(framework)
add_subdirectory(modules)
add_subdirectory(cru)
add_subdirectory(libraries)
if (WIN32)
  add_subdirectory(windows_version)
endif()

add_subdirectory(tests)

# Add the command line program's target
if (WIN32)
  # Let the exe be called guesscmd so it doesn't collide with the dll target
  set(guess_command_name "guesscmd")
else()
  # On Unix we don't have the dll so the command line binary can be called guess
  set(guess_command_name "guess")
endif()

# Build the common model code once, then include all of its objects in the
# executables below. This ensures that static input/output module registration
# initializers are not discarded by the static-library linker.
add_library(guess_core OBJECT ${guess_sources})

# Specify the executable to build, and which sources to build it from
add_executable(${guess_command_name}
  command_line_version/main.cpp
  $<TARGET_OBJECTS:guess_core>)
target_link_libraries(${guess_command_name} ${LIBS})

# Rule for building the unit test binary
add_executable(runtests ${test_sources} $<TARGET_OBJECTS:guess_core>)
target_link_libraries(runtests ${LIBS})

if (WIN32)
  # Create guess.dll (used with the graphical Windows shell)
  add_library(guess SHARED
    windows_version/dllmain.cpp
    $<TARGET_OBJECTS:guess_core>)

  # Specify libraries to link to the dll
  target_link_libraries(guess ${LIBS})
endif (WIN32)

# The custom build rule for generating the submit script from a template.
# The submit script is generated each time the command line binary is built.
if (UNIX)
  # Choose template file depending on system
  if (SYSTEM STREQUAL "")
    set(SYSTEM "simba")
  endif(SYSTEM STREQUAL "")  
  set(submit_template ${guess_SOURCE_DIR}/parallel_version/${SYSTEM}.tmpl)

  # A sed command which exchanges each occurrence of BINARY and DATE in
  # the template with appropriate path to the compiled binary and the
  # current date.
  add_custom_command(TARGET ${guess_command_name}
                     POST_BUILD
             COMMAND sed 's|BINARY|'${CMAKE_CURRENT_BINARY_DIR}/${guess_command_name}'|g\;s|DATE|'\"`date`\"'|g' ${submit_template} > ${CMAKE_BINARY_DIR}/submit.sh
			 COMMAND chmod +x ${CMAKE_BINARY_DIR}/submit.sh
             COMMENT "Generating submit script")
endif (UNIX)

# Rule for running unit tests automatically
if (UNIT_TESTS)
  add_custom_command(TARGET runtests 
    POST_BUILD 
    COMMAND runtests 
    COMMENT "Running tests")
endif()

if (UNIX)
  # Set default build type to Release on Unix
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
  endif(NOT CMAKE_BUILD_TYPE)
endif(UNIX)

if(APPLE)
  add_custom_command(TARGET guess POST_BUILD
    COMMAND codesign --force --sign - 
            --entitlements ${CMAKE_SOURCE_DIR}/debug.entitlements
            $<TARGET_FILE:guess>
    COMMENT "Ad-hoc signing LPJ-GUESS with debug entitlement"
  )
endif()
