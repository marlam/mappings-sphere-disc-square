# Copyright (C) 2015
# Computer Graphics Group, University of Siegen
# Written by Martin Lambers <martin.lambers@uni-siegen.de>
#
# Copying and distribution of this file, with or without modification, are
# permitted in any medium without royalty provided the copyright notice and this
# notice are preserved. This file is offered as-is, without any warranty.

cmake_minimum_required(VERSION 2.8)
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

if(CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++1y -Wall -Wextra")
endif()

# Required libraries
find_package(Qt5Gui REQUIRED)

# The projection library
add_subdirectory(libproj)

# The programs
include_directories(${CMAKE_SOURCE_DIR}/libproj)

add_executable(disk-square-precision-test disk-square-precision-test.cpp disk-square-common.hpp)
target_link_libraries(disk-square-precision-test libproj)
install(TARGETS disk-square-precision-test RUNTIME DESTINATION bin)

add_executable(disk-square-analyzer disk-square-analyzer.cpp disk-square-common.hpp pngvis.hpp pngvis.cpp)
target_link_libraries(disk-square-analyzer libproj Qt5::Gui)
install(TARGETS disk-square-analyzer RUNTIME DESTINATION bin)

add_executable(disk-square-plotter disk-square-plotter.cpp disk-square-common.hpp pdfplot.hpp pdfplot.cpp)
target_link_libraries(disk-square-plotter libproj Qt5::Gui)
install(TARGETS disk-square-plotter RUNTIME DESTINATION bin)

add_executable(sphere-map-precision-test sphere-map-precision-test.cpp sphere-map-common.hpp)
target_link_libraries(sphere-map-precision-test libproj)
install(TARGETS sphere-map-precision-test RUNTIME DESTINATION bin)

add_executable(sphere-map-analyzer sphere-map-analyzer.cpp sphere-map-common.hpp pngvis.hpp pngvis.cpp)
target_link_libraries(sphere-map-analyzer libproj Qt5::Gui)
install(TARGETS sphere-map-analyzer RUNTIME DESTINATION bin)

qt5_add_resources(SPHERE_MAP_PLOTTER_RESOURCES sphere-map-plotter.qrc)
add_executable(sphere-map-plotter sphere-map-plotter.cpp sphere-map-common.hpp
	pdfplot.hpp pdfplot.cpp ${SPHERE_MAP_PLOTTER_RESOURCES})
target_link_libraries(sphere-map-plotter libproj Qt5::Gui)
install(TARGETS sphere-map-plotter RUNTIME DESTINATION bin)

# The data
#configure_file(land-polygons-z5.csv land-polygons-z5.csv COPYONLY)
