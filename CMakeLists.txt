cmake_minimum_required(VERSION 3.1)

project(embedded_deformation)

include(FetchContent)
FetchContent_Declare(
    polyscope

    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
    GIT_TAG        56b82cd
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/polyscope"
)
FetchContent_Populate(polyscope)

include(FetchContent)
FetchContent_Declare(
    utils

    GIT_REPOSITORY https://github.com/rFalque/utils.git
    GIT_TAG        master
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/utils"
)
FetchContent_Populate(utils)

FetchContent_Declare(
    libGraphCpp

    GIT_REPOSITORY https://github.com/rFalque/libGraphCpp.git
    GIT_TAG        library
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/libGraphCpp"
)
FetchContent_Populate(libGraphCpp)

add_subdirectory(utils)
add_subdirectory(polyscope)
add_subdirectory(libGraphCpp)
add_subdirectory(embedded_deformation)



set(CMAKE_BUILD_TYPE Release)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# specific version of eigen needed
find_package(Eigen3 REQUIRED)

# yaml-cpp if not found, run sudo apt-get install libyaml-cpp-dev
find_package(yaml-cpp REQUIRED)

find_package(Ceres REQUIRED)


file(GLOB_RECURSE my_c_list RELATIVE ${CMAKE_SOURCE_DIR} "app/*")

foreach(file_path ${my_c_list})
    string( REPLACE ".cpp" "" new_name ${file_path} )
    get_filename_component(filename ${new_name} NAME)
    add_executable( ${filename} ${file_path})
    set_target_properties(${filename} PROPERTIES CXX_STANDARD 11 CXX_STANDARD_REQUIRED YES)
    include_directories(${filename}
                        ${EIGEN3_INCLUDE_DIR}
                        ${YAML_CPP_INCLUDE_DIR}
                        ${CERES_INCLUDE_DIRS}
                        ${PYTHON_INCLUDE_DIRS}
                        )

    target_link_libraries(${filename}
                        utils
                        polyscope
                        libGraphCpp
                        embedded_deformation
                        ${CERES_LIBRARIES}
                        ${YAML_CPP_LIBRARIES}
                        ${PYTHON_LIBRARIES}
                        )

endforeach()
