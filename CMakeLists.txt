# MIT License

# Copyright (c) 2014 Jay McMahon

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

# @file   CMakeLists.txt
# @Author Benjamin Bercovici (bebe0705@colorado.edu)
# @date   March, 2018
# @brief  CMake listing enabling compilation and installation of the YORPLib library

################################################################################
#
#
# 		The following should normally not require any modification
# 				Unless new files are added to the build tree
#
#
################################################################################

cmake_minimum_required(VERSION 3.0)


get_filename_component(dirName ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(LIB_NAME "YORPLib" CACHE STRING "Name of library to be created.")

project(${LIB_NAME})

# Specify the version used
if (${CMAKE_MAJOR_VERSION} LESS 3)
	message(FATAL_ERROR " You are running an outdated version of CMake")
endif()

cmake_minimum_required(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.0)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/source/cmake)

# Compiler flags
add_definitions(-Wall -O2 )

# Enable C++17 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

# Include directories
include_directories(include)

# Find OpenMP
if (NOT NO_OMP)
	find_package(OpenMP)
endif()


# Add source files in root directory
add_library(${LIB_NAME} 
	SHARED
	source/Body.cpp
	source/Facet.cpp
	source/FCoeffs.cpp
	source/Shadow.cpp
	source/SRPModel.cpp
	source/VoxelGrid.cpp
	)


if (OpenMP_FOUND AND NOT NO_OMP)
	target_link_libraries(${LIB_NAME} ${library_dependencies} OpenMP::OpenMP_CXX)
else()
	target_link_libraries(${LIB_NAME} ${library_dependencies})
endif()


# If ${BREW} is true, installing from homebrew formula
if(NOT BREW)
	install (TARGETS ${LIB_NAME} DESTINATION /usr/local/lib/)
	install (DIRECTORY ${PROJECT_SOURCE_DIR}/share/${LIB_NAME} DESTINATION /usr/local/share)
	install (DIRECTORY ${PROJECT_SOURCE_DIR}/include/YORPLib DESTINATION /usr/local/include/ FILES_MATCHING PATTERN "*.h")
endif()








