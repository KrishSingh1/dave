
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/generated/config.hpp"
  @ONLY
)

set(headers 
  catch.hpp
  temporary_directory.h
  simulation_helper.h
  )

set(source  
  main.cpp
  date_test.cpp
  climate_test.cpp
  math_test.cpp
  kdtree_test.cpp
  indata_test.cpp
  ncompete_test.cpp
  cftime_test.cpp
  string_test.cpp
  guesscontainer_test.cpp
  management_test.cpp
  demo_integration_test.cpp
  demo_state_test.cpp
  fluxes_test.cpp
  temporary_directory.cpp
  simulation_helper.cpp
  )

include(add_test_sources)
