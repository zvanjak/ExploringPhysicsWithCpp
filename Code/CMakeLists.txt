cmake_minimum_required(VERSION 3.20)

set (CMAKE_CXX_STANDARD 20)

project(Demo_Chapter_01 VERSION 1.0)
set(DEMO_CHAPTER_01_SOURCES  Chapter_01_basic_objects/main_chapter_01.cpp
                            Chapter_01_basic_objects/demo_baseutils.cpp
                            Chapter_01_basic_objects/demo_geometry.cpp
                            Chapter_01_basic_objects/demo_geometry_2d.cpp
                            Chapter_01_basic_objects/demo_geometry_3d.cpp

                      )
add_executable(Demo_Chapter_01 ${DEMO_CHAPTER_01_SOURCES})
target_include_directories(Demo_Chapter_01 PRIVATE Chapter_01_basic_objects/mml)

project(Demo_Chapter_02 VERSION 1.0)
set(DEMO_CHAPTER_02_SOURCES  Chapter_02_basic_algorithms/main_chapter_02.cpp 
                      )
add_executable(Demo_Chapter_02 ${DEMO_CHAPTER_02_SOURCES})
target_include_directories(Demo_Chapter_02 PRIVATE Chapter_02_basic_algorithms/mml)