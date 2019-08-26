# Copyright (C) 2015
#
# Copying and distribution of this file, with or without modification, are
# permitted in any medium without royalty provided the copyright notice and this
# notice are preserved. This file is offered as-is, without any warranty.

cmake_minimum_required(VERSION 2.8)

if(CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++1y -Wall -Wextra")
endif()

#find_package(GSL REQUIRED)
#find_package(Boost REQUIRED COMPONENTS math_tr1)

#include_directories(${CMAKE_SOURCE_DIR} ${GSL_INCLUDE_DIRS})
#include_directories(${CMAKE_SOURCE_DIR} ${Boost_INCLUDE_DIRS})
add_library(libproj STATIC
	proj.cpp proj.hpp
        disk-square.cpp disk-square.hpp
	sphere-lambert-azimuthal-equal-area.cpp sphere-lambert-azimuthal-equal-area.hpp
	sphere-breusing-harmonic.cpp sphere-breusing-harmonic.hpp
	sphere-equidistant-azimuthal.cpp sphere-equidistant-azimuthal.hpp
	hemisphere-stereographic.cpp hemisphere-stereographic.hpp
	hemisphere-mix-lambert-stereographic.cpp hemisphere-mix-lambert-stereographic.hpp
        hemisphere-perspective.cpp hemisphere-perspective.hpp)
#target_link_libraries(libproj ${GSL_LIBRARIES})
#target_link_libraries(libproj ${Boost_LIBRARIES})
#target_link_libraries(libproj -lquadmath)
set_target_properties(libproj PROPERTIES OUTPUT_NAME proj)
install(TARGETS libproj RUNTIME DESTINATION bin LIBRARY DESTINATION "lib" ARCHIVE DESTINATION "lib")
