cmake_minimum_required(VERSION 3.1)

project(embedded_deformation)

### automatically download submodules
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external_libraries/libGraphCpp/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

add_subdirectory(external_libraries)
include_directories(utils)
include_directories(include)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# specific version of eigen needed
find_package(Eigen3 REQUIRED)

# yaml-cpp if not found, run sudo apt-get install libyaml-cpp-dev
find_package(yaml-cpp REQUIRED)

find_package(Ceres REQUIRED)

set(SOURCE_FILES include/embedded_deformation/greedySearch.cpp
		 include/embedded_deformation/embedDeform.cpp)

file(GLOB_RECURSE my_c_list RELATIVE ${CMAKE_SOURCE_DIR} "app/*")

foreach(file_path ${my_c_list})
    string( REPLACE ".cpp" "" new_name ${file_path} )
    get_filename_component(filename ${new_name} NAME)
    add_executable( ${filename} ${file_path} ${SOURCE_FILES})
    set_target_properties(${filename} PROPERTIES CXX_STANDARD 11 CXX_STANDARD_REQUIRED YES)
    include_directories(${filename}
                        ${PROJECT_SOURCE_DIR}/external_libraries/
                        ${YAML_CPP_INCLUDE_DIR}
			${CERES_INCLUDE_DIRS}
                        ${EIGEN3_INCLUDE_DIR}
                        ${PYTHON_INCLUDE_DIRS}
			${LIB_GRAPH_CPP_INCLUDE_DIR}
			)

    target_link_libraries(${filename}
                        polyscope
                        ${PYTHON_LIBRARIES}
			${CERES_LIBRARIES}
 			${YAML_CPP_LIBRARIES}
			)

endforeach()
