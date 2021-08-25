cmake_minimum_required(VERSION 3.1)

project(embedded_deformation)

SET(INCLUDE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/include/")
SET(SRC_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/src/")

find_package(Eigen3 REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(Ceres REQUIRED)

set(SRCS
    # embedded_deformation
    ${SRC_ROOT}/embedded_deformation/embedDeform.cpp
    ${SRC_ROOT}/embedded_deformation/greedySearch.cpp
)

set(HEADERS
    # embedded_deformation
    ${INCLUDE_ROOT}/embedded_deformation/costFunction.hpp
    ${INCLUDE_ROOT}/embedded_deformation/downsampling.hpp
    ${INCLUDE_ROOT}/embedded_deformation/embedDeform.hpp
    ${INCLUDE_ROOT}/embedded_deformation/farther_sampling.hpp
    ${INCLUDE_ROOT}/embedded_deformation/getMinMax.hpp
    ${INCLUDE_ROOT}/embedded_deformation/greedySearch.hpp
    ${INCLUDE_ROOT}/embedded_deformation/nanoflann.hpp
    ${INCLUDE_ROOT}/embedded_deformation/nanoflannWrapper.hpp
    ${INCLUDE_ROOT}/embedded_deformation/options.hpp
)

add_library(embedded_deformation ${SRCS} ${HEADERS})

target_include_directories(embedded_deformation PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_include_directories(embedded_deformation PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include/embedded_deformation")
target_include_directories(embedded_deformation PRIVATE ${EIGEN3_INCLUDE_DIR})
target_include_directories(embedded_deformation PRIVATE ${YAML_CPP_INCLUDE_DIR})
target_include_directories(embedded_deformation PRIVATE ${CERES_INCLUDE_DIRS})
target_include_directories(embedded_deformation PRIVATE ${LIB_GRAPH_CPP_INCLUDE_DIR})

target_link_libraries(embedded_deformation libGraphCpp)
target_link_libraries(embedded_deformation ${CERES_LIBRARIES} )
target_link_libraries(embedded_deformation ${YAML_CPP_LIBRARIES} )
