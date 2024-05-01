cmake_minimum_required(VERSION 3.8)
project(gfq_ldpc CXX)
#C++ standard option 
set(CMAKE_CXX_STANDARD 23)

# set current directory for include files
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# workaround for DOWNLOAD_EXTERNAL_TIMESTAMP option
if(POLICY CMP0135)
	cmake_policy(SET CMP0135 NEW)
	set(CMAKE_POLICY_DEFAULT_CMP0135 NEW)
endif()

### Configure dependency install path, include path and library path
set(DEPEND_LIB_ROOT ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_BINARY_DIR})
include_directories(${PROJECT_BINARY_DIR}/include)
link_directories(${DEPEND_LIB_ROOT})

set(CMAKE_CXX_FLAGS "-O3 -Wall -g -fPIC -fmessage-length=0 -DDEBUG -I/opt/homebrew/include -I${PROJECT_BINARY_DIR}/include -L/opt/homebrew/lib")

set(galois_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/galois)
set(galois_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(galois_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${galois_INSTALL_DIR})

message("galois_PREFIX = '${galois_PREFIX}'")
message("galois_INSTALL_DIR = '${galois_INSTALL_DIR}'")

set(fftw3_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fftw3)
set(fftw3_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(fftw3_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/fftw3/build)
set(fftw3_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${fftw3_INSTALL_DIR})
set(fftw3_VERSION 3.3.10)


set(Boost_USE_STATIC_LIBS    ON)
set(Boost_USE_MULTITHREADED  ON)
set(Boost_INCLUDE_DIR "/opt/homebrew/include/boost")
set(Boost_LIBRARY_DIR "/opt/homebrew/lib")
 
find_package(Boost 1.85.0 REQUIRED)
message(STATUS "  Boost_INCLUDE_DIR=${Boost_INCLUDE_DIR}")
message(STATUS "  Boost_LIBRARY_DIR=${Boost_LIBRARY_DIR}")
message(STATUS "  Boost_LIBRARIES=${Boost_LIBRARIES}")

include(ExternalProject)
ExternalProject_Add(galois
### Set prefix
    PREFIX ${CMAKE_PREFIX_PATH}/galois
### Path to Git repository
    GIT_REPOSITORY https://github.com/wkjarosz/galois.git
### Specify Tag(dummy option)
#    GIT_TAG 3.1.2
    INSTALL_DIR ${galois_INSTALL_DIR}
    CMAKE_ARGS ${galois_CMAKE_ARGS}
### BUILD_COMMAND
    BUILD_COMMAND make
### install command
    INSTALL_COMMAND cp ${galois_PREFIX}/src/galois-build/libgalois++.a ${galois_INSTALL_DIR}/ && 
                    cp -r ${galois_PREFIX}/src/galois/include ${PROJECT_BINARY_DIR}/
)

ExternalProject_Add(fftw3
### Set prefix
    PREFIX ${fftw3_PREFIX}
### Set URL
    URL https://www.fftw.org/fftw-${fftw3_VERSION}.tar.gz
### Set MD5
    URL_MD5 8ccbf6a5ea78a16dbc3e1306e234cc5c
### Set source directory
    SOURCE_DIR ${fftw3_PREFIX}/src
### prohibit extraction
    DOWNLOAD_NO_EXTRACT true
### Set download directory
    TMP_DIR ${fftw3_PREFIX}/tmp
    DOWNLOAD_DIR ${fftw3_PREFIX}/tarballs
    BINARY_DIR ${fftw3_BINARY_DIR}
    INSTALL_DIR ${fftw3_INSTALL_DIR}
    CMAKE_ARGS ${galois_CMAKE_ARGS}
### CONFIGURE_COMMAND
    CONFIGURE_COMMAND cd ${fftw3_PREFIX}/src && tar xvfz ${fftw3_PREFIX}/tarballs/fftw-${fftw3_VERSION}.tar.gz && cd ${fftw3_PREFIX}/src/fftw-${fftw3_VERSION} && ./configure --enable-threads && make
### BUILD_COMMAND
    BUILD_COMMAND ""
### UPDATE COMMAND
    UPDATE_COMMAND ""
### install command
    INSTALL_COMMAND cp ${fftw3_PREFIX}/src/fftw-${fftw3_VERSION}/.libs/libfftw3.a ${fftw3_INSTALL_DIR}/ && 
                    cp -r ${fftw3_PREFIX}/src/fftw-${fftw3_VERSION}/api/fftw3.h ${PROJECT_BINARY_DIR}/include
)

include_directories("/usr/include" "/opt/homebrew/include" "${Boost_INCLUDE_DIRS}" "${PROJECT_BINARY_DIR}")
link_directories("/usr/lib" "/opt/homebrew/lib" "${Boost_LIBRARY_DIR}")

add_executable(make-gen make-gen.cpp get_generating_matrix.cpp generating_matrix.hpp gfq_alist.hpp split.cpp split.hpp)
add_dependencies(make-gen galois)
target_link_libraries(make-gen ${DEPEND_LIB_ROOT}/libgalois++.a m)

add_executable(gfq_simulator gfq_simulator.cpp decode.cpp encode.cpp transmit.cpp get_generating_matrix.cpp generating_matrix.hpp gfq_alist.hpp edge.hpp variable_node.hpp factor_node.hpp transmit.cpp factor_graph.cpp split.cpp split.hpp)
add_dependencies(gfq_simulator galois fftw3)
target_link_libraries(gfq_simulator ${DEPEND_LIB_ROOT}/libgalois++.a ${DEPEND_LIB_ROOT}/libfftw3.a m)