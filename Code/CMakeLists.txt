cmake_minimum_required(VERSION 3.20)

set (CMAKE_CXX_STANDARD 20)

project(Chapter_01_Basic_objects VERSION 1.0)
set(DEMO_CHAPTER_01_SOURCES  Chapter_01_basic_objects/main_chapter_01.cpp
                             Chapter_01_basic_objects/demo_baseutils.cpp
                             Chapter_01_basic_objects/demo_geometry.cpp
                             Chapter_01_basic_objects/demo_geometry_2d.cpp
                             Chapter_01_basic_objects/demo_geometry_3d.cpp
                             Chapter_01_basic_objects/demo_vector.cpp
                             Chapter_01_basic_objects/demo_matrix.cpp
                             Chapter_01_basic_objects/demo_matrixnm.cpp
                             Chapter_01_basic_objects/demo_functions.cpp
                             Chapter_01_basic_objects/demo_real_function.cpp
                      )
add_executable(Chapter_01_Basic_objects ${DEMO_CHAPTER_01_SOURCES} )
target_include_directories(Chapter_01_Basic_objects PRIVATE mml)


project(Chapter_02_Vizualization VERSION 1.0)
set(DEMO_CHAPTER_02_SOURCES "Chapter_02_vizualization/main_chapter_02.cpp" 
                            Chapter_02_vizualization/example2_wpf_visualizators.cpp
                      )
add_executable(Chapter_02_Vizualization ${DEMO_CHAPTER_02_SOURCES} )
target_include_directories(Chapter_02_Vizualization PRIVATE mml)


project(Chapter_03_Basic_alg VERSION 1.0)
set(DEMO_CHAPTER_03_SOURCES "Chapter_03_basic_algorithms/main_chapter_03.cpp" 
                            Chapter_03_basic_algorithms/demo_derivation.cpp
                            Chapter_03_basic_algorithms/demo_integration.cpp
                            Chapter_03_basic_algorithms/demo_integration_2d.cpp
                            Chapter_03_basic_algorithms/demo_integration_3d.cpp
                            Chapter_03_basic_algorithms/demo_interpolation.cpp
                            Chapter_03_basic_algorithms/demo_lin_alg_sys_solvers.cpp
                            Chapter_03_basic_algorithms/demo_root_finding.cpp
                      )
add_executable(Chapter_03_Basic_alg ${DEMO_CHAPTER_03_SOURCES} )
target_include_directories(Chapter_03_Basic_alg PRIVATE mml)

project(Chapter_04_Collision_calc VERSION 1.0)
set(DEMO_CHAPTER_04_SOURCES Chapter_04_collision_simulator/main_chapter_04.cpp 
                            Chapter_04_collision_simulator/collision_simulator_2d.cpp
                            Chapter_04_collision_simulator/collision_simulator_3d.cpp
                      )
add_executable(Chapter_04_Collision_calc ${DEMO_CHAPTER_04_SOURCES} )
target_include_directories(Chapter_04_Collision_calc PRIVATE mml)
target_include_directories(Chapter_04_Collision_calc PRIVATE mpl)


project(Chapter_05_Throwing_things VERSION 1.0)
set(DEMO_CHAPTER_05_SOURCES Chapter_05_throwing_thing_in_the_air/main_chapter_05.cpp )

add_executable(Chapter_05_Throwing_things ${DEMO_CHAPTER_05_SOURCES} )
target_include_directories(Chapter_05_Throwing_things PRIVATE mml)
target_include_directories(Chapter_05_Throwing_things PRIVATE mpl)

project(Chapter_06_Pendulum VERSION 1.0)
set(DEMO_CHAPTER_06_SOURCES Chapter_06_pendulum/main_chapter_06.cpp 
                      )
add_executable(Chapter_06_Pendulum ${DEMO_CHAPTER_06_SOURCES} )
target_include_directories(Chapter_06_Pendulum PRIVATE mml)
target_include_directories(Chapter_06_Pendulum PRIVATE mpl)


project(Chapter_07_Double_spherical_pendulum VERSION 1.0)
set(DEMO_CHAPTER_07_SOURCES Chapter_07_double_and_spherical_pendulum/main_chapter_07.cpp 
                      )
add_executable(Chapter_07_Double_spherical_pendulum ${DEMO_CHAPTER_07_SOURCES} )
target_include_directories(Chapter_07_Double_spherical_pendulum PRIVATE mml)
target_include_directories(Chapter_07_Double_spherical_pendulum PRIVATE mpl)


project(Chapter_08_Gravity VERSION 1.0)
set(DEMO_CHAPTER_08_SOURCES Chapter_08_gravity/main_chapter_08.cpp 
                      )
add_executable(Chapter_08_Gravity  ${DEMO_CHAPTER_08_SOURCES} )
target_include_directories(Chapter_08_Gravity PRIVATE mml)
target_include_directories(Chapter_08_Gravity PRIVATE mpl)


project(Chapter_09_NBody_problem VERSION 1.0)
set(DEMO_CHAPTER_09_SOURCES Chapter_09_gravity_N_body_problem/main_chapter_09.cpp 
                      )
add_executable(Chapter_09_NBody_problem ${DEMO_CHAPTER_09_SOURCES} )
target_include_directories(Chapter_09_NBody_problem PRIVATE mml)
target_include_directories(Chapter_09_NBody_problem PRIVATE mpl)


project(Chapter_13_Rigid_body VERSION 1.0)
set(DEMO_CHAPTER_13_SOURCES Chapter_13_rigid_body/main_chapter_13.cpp 
                      )
add_executable(Chapter_13_Rigid_body ${DEMO_CHAPTER_13_SOURCES} )
target_include_directories(Chapter_13_Rigid_body PRIVATE mml)
target_include_directories(Chapter_13_Rigid_body PRIVATE mpl)

project(Chapter_18_Static_magnetic_fields VERSION 1.0)
set(DEMO_CHAPTER_18_SOURCES Chapter_18_static_magnetic_fields/main_chapter_18.cpp 
                      )
add_executable(Chapter_18_Static_magnetic_fields ${DEMO_CHAPTER_18_SOURCES} )
target_include_directories(Chapter_18_Static_magnetic_fields PRIVATE mml)