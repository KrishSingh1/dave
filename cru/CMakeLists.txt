file(GLOB headers "guessio/*.h")
file(GLOB source  "guessio/*.cpp")

set(guess_sources ${guess_sources} ${headers} ${source} PARENT_SCOPE)
