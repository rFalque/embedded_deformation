cmake_minimum_required(VERSION 3.1)

project(embedded_deformation)

set(CMAKE_CXX_STANDARD 14)          
set(CMAKE_CXX_STANDARD_REQUIRED ON)  
set(CMAKE_CXX_EXTENSIONS OFF) 

SET(INCLUDE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/include")
SET(SRC_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/src")

find_package(Eigen3 REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(Ceres REQUIRED)
find_package(OpenCV REQUIRED)
find_package(CUDA REQUIRED)

include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/polyscope/include)

file(GLOB containers
    ${INCLUDE_ROOT}/containers/*.hpp
    ${INCLUDE_ROOT}/containers/*.cuh
    ${INCLUDE_ROOT}/containers/*.h
)
cuda_add_library(containers ${containers})

file(GLOB SRCS
    # embedded_deformation
    ${SRC_ROOT}/embedded_deformation/*.cpp
    ${SRC_ROOT}/imageProcessor/*.cu
    ${SRC_ROOT}/imageProcessor/*.cpp
    ${SRC_ROOT}/math/*.cpp
    ${SRC_ROOT}/rigidSolver/*.cpp
    ${SRC_ROOT}/rigidSolver/*.cu
)

file(GLOB HEADERS
    # embedded_deformation
    ${INCLUDE_ROOT}/embedded_deformation/*.hpp
    ${INCLUDE_ROOT}/imageProcessor/*.h
    ${INCLUDE_ROOT}/math/*.h
    ${INCLUDE_ROOT}/math/*.hpp
    ${INCLUDE_ROOT}/containers/*.h
    ${INCLUDE_ROOT}/containers/*.hpp
    ${INCLUDE_ROOT}/containers/*.cpp
    ${INCLUDE_ROOT}/containers/*.cuh
    ${INCLUDE_ROOT}/rigidSolver/*.h
    
)

cuda_add_library(embedded_deformation ${SRCS} ${HEADERS})

target_include_directories(embedded_deformation PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_include_directories(embedded_deformation PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include/embedded_deformation")
target_include_directories(embedded_deformation PRIVATE ${EIGEN3_INCLUDE_DIR})
target_include_directories(embedded_deformation PRIVATE ${YAML_CPP_INCLUDE_DIR})
target_include_directories(embedded_deformation PRIVATE ${CERES_INCLUDE_DIRS})
target_include_directories(embedded_deformation PRIVATE ${LIB_GRAPH_CPP_INCLUDE_DIR})

target_link_libraries(embedded_deformation libGraphCpp)
target_link_libraries(embedded_deformation ${CERES_LIBRARIES} )
target_link_libraries(embedded_deformation ${YAML_CPP_LIBRARIES} ${OpenCV_LIBRARIES} ${CUDA_LIBRARIES} containers)
