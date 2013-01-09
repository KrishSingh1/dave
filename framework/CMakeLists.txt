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
  commandlinearguments.h
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
  commandlinearguments.cpp
)

include(add_guess_sources)