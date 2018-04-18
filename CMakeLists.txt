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
# 								User-defined paths
#						Should be checked for consistency
#						Before running 'cmake ..' in build dir
#
################################################################################

# OMP-friendly GCC Major version (will only be used on Mac OS systems)
set (OMP_FRIENDLY_GCC_MAJOR_VERSION 7)

# Path to OMP-Friendly GCC Compiler (will only be used on Mac OS systems)
set (OMP_FRIENDLY_GCC_PATH /usr/local/Cellar/gcc/7.3.0_1/bin/)


################################################################################
#
#
# 		The following should normally not require any modification
# 				Unless new files are added to the build tree
#
#
################################################################################

# Building procedure
get_filename_component(dirName ${CMAKE_CURRENT_SOURCE_DIR} NAME)
set(LIB_NAME YORPLib)

if(APPLE)
	if(EXISTS ${OMP_FRIENDLY_GCC_PATH}) 
		# Else running on a MAC. Will attempt to switch compiler to get the OMP-friendly GCC 6.3.0_1 from Homebrew
		message("Switching to OMP-friendly GCC ")
		set(CMAKE_C_COMPILER ${OMP_FRIENDLY_GCC_PATH}gcc-${OMP_FRIENDLY_GCC_MAJOR_VERSION} CACHE STRING "C Compiler" FORCE)
		set(CMAKE_CXX_COMPILER ${OMP_FRIENDLY_GCC_PATH}g++-${OMP_FRIENDLY_GCC_MAJOR_VERSION} CACHE STRING "C++ Compiler" FORCE)

	else()
		# no OMP-compliant compiler was found on this mac. 
		message( "No OMP-compliant compiler was found on this Mac.")

	endif()
else() 
	# Running on Linux. Will switch back to compiler in /usr/local/bin
	message("Switching to /usr/local/gcc ")
	set(CMAKE_C_COMPILER "/usr/local/bin/gcc" CACHE STRING "C Compiler" FORCE)
	set(CMAKE_CXX_COMPILER "/usr/local/bin/g++" CACHE STRING "C++ Compiler" FORCE)
endif()

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
# find_package(OpenMP REQUIRED)
find_package(OpenMP)
if(OpenMP_C_FLAGS)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
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


# If ${BREW} is true, installing from homebrew formula
if(NOT BREW)
	install (TARGETS ${LIB_NAME} DESTINATION /usr/local/lib/)
	install (FILES ${PROJECT_SOURCE_DIR}/${LIB_NAME} DESTINATION /usr/local/share)
	install (DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION /usr/local/include/YORPLib FILES_MATCHING PATTERN "*.h")
endif()







