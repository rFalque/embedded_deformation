cmake_minimum_required(VERSION 2.8)

project(embedded_deformation)

include_directories(utils)
include_directories(include)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# specific version of eigen needed
find_package(Eigen3 REQUIRED)

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
	            		${CERES_INCLUDE_DIRS}
                        ${EIGEN3_INCLUDE_DIR}
                        ${PYTHON_INCLUDE_DIRS}
                        )

    target_link_libraries(${filename}
                        ${PYTHON_LIBRARIES}
                        ${CERES_LIBRARIES}
                        )

endforeach()
