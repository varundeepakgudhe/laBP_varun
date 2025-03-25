cmake_minimum_required(VERSION "3.26.5")

project("ARGSIM")

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set optimization flag
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")


# Include paths
include_directories(
    /Library/Developer/CommandLineTools/SDKs/MacOSX15.2.sdk/usr/include/c++/v1/
    /opt/homebrew/Cellar/boost/1.87.0/include
)

# Link directories
link_directories(
    /opt/homebrew/Cellar/boost/1.87.0/lib
)

set(SOURCES
    src/main.cpp 
    src/argnode.cpp 
    src/chromosome.cpp 
    src/chromrecomb.cpp 
    src/parameters.cpp 
    src/poisevents.cpp 
    src/simulate.cpp 
    src/sitenode.cpp 
    src/world.cpp
    src/migprob.cpp
)

# Build executable
add_executable(labp_v31 ${SOURCES})

# Set the output directory for the executable
set_target_properties(labp_v31 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/../executables")