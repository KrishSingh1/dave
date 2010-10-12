file(GLOB headers "*.h")
file(GLOB source  "*.cpp")

set(guess_sources ${guess_sources} ${headers} ${source} PARENT_SCOPE)