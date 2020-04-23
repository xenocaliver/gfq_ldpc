cmake_minimum_required(VERSION 3.8)
project(gfq_ldpc CXX)
# 最適化・警告等のオプション
set(CMAKE_CXX_STANDARD 17)

# curdirをinclude file探索対象にする
set(CMAKE_INCLUDE_CURRENT_DIR ON)

### Configure dependency install path, include path and library path
set(DEPEND_LIB_ROOT ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_BINARY_DIR})
include_directories(${PROJECT_BINARY_DIR}/include)
link_directories(${DEPEND_LIB_ROOT})

set(CMAKE_CXX_FLAGS "-O3 -Wall -g -fPIC -fmessage-length=0 -I/usr/local/include -I${PROJECT_BINARY_DIR}/include -L/usr/local/lib")

set(galois_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/galois)
set(galois_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(galois_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${galois_INSTALL_DIR})

message("galois_PREFIX = '${galois_PREFIX}'")
message("galois_INSTALL_DIR = '${galois_INSTALL_DIR}'")

set(Boost_USE_STATIC_LIBS    ON)
set(Boost_USE_MULTITHREADED  ON)
set(Boost_INCLUDE_DIR "/usr/local/include/boost")
set(Boost_LIBRARY_DIR "/usr/local/lib")
 
find_package(Boost 1.72.0 REQUIRED)
message(STATUS "  Boost_INCLUDE_DIR=${Boost_INCLUDE_DIR}")
message(STATUS "  Boost_LIBRARY_DIR=${Boost_LIBRARY_DIR}")
message(STATUS "  Boost_LIBRARIES  =${Boost_LIBRARIES}")

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

include_directories("/usr/include" "/usr/local/include" "${Boost_INCLUDE_DIRS}")
link_directories("/usr/lib" "/usr/local/lib" "${Boost_LIBRARY_DIR}")

add_executable(make-gen make-gen.cpp get_generating_matrix.cpp generating_matrix.hpp gfq_alist.hpp)
add_dependencies(make-gen galois)
target_link_libraries(make-gen ${DEPEND_LIB_ROOT}/libgalois++.a Boost::boost)