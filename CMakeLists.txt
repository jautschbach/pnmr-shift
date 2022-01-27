cmake_minimum_required(VERSION 2.8)
project(PNMRShift CXX)
enable_language(Fortran)

include_directories(${CMAKE_SOURCE_DIR})

#set(Boost_DEBUG TRUE)

add_definitions("-Wall")

IF(MINGW)
  SET(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++")
ENDIF(MINGW)

#find_package(Boost)
#if(Boost_FOUND) 
#  include_directories( ${Boost_INCLUDE_DIRS} )
#endif(Boost_FOUND)

if (CMAKE_BUILD_TYPE STREQUAL "debug")
    add_definitions(-O0 -g -DTHROWSTREAM_EXCEPTIONSOURCE)
endif (CMAKE_BUILD_TYPE STREQUAL "debug")

if (CMAKE_BUILD_TYPE STREQUAL "release")
    add_definitions(-O2 -funroll-loops)
endif (CMAKE_BUILD_TYPE STREQUAL "release")

add_subdirectory(lapack)

set(PNMRSHIFTFILES ConvertSelectionString.cpp Output.cpp PNMRShift.cpp
ReadFile.cpp XYZFile.cpp)
 
add_executable(PNMRShift ${PNMRSHIFTFILES})
target_link_libraries(PNMRShift lapack_zheev)
                    
