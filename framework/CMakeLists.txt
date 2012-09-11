set(headers 
  config.h 
  guess.h 
  guessmath.h
  outputchannel.h 
  archive.h 
  framework.h 
  shell.h
  partitionedmapserializer.h
  guessserializer.h
  parallel.h
)

set(source  
  guess.cpp 
  outputchannel.cpp 
  archive.cpp 
  framework.cpp 
  shell.cpp
  partitionedmapserializer.cpp
  guessserializer.cpp
  parallel.cpp
)

include(add_guess_sources)