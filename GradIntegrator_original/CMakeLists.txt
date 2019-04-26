CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

add_definitions(-std=c++11)
if(POLICY CMP0012)
  cmake_policy(SET CMP0012 NEW)
endif()

PROJECT(GInt)

ADD_SUBDIRECTORY(src)

option(VTK_SUPPORT_ENABLED "Enable VTK" FALSE)
if(VTK_SUPPORT_ENABLED)
  ADD_DEFINITIONS(-DVTK_ENABLED)
endif()

if(CMAKE_GENERATOR STREQUAL Xcode)
  set(APPLE_CLANG_COMPILER /usr/local/opt/llvm/bin/clang)
  set(CMAKE_CXX_FLAGS "-fopenmp -std=c++11")
  link_directories(/usr/local/Cellar/llvm/5.0.1/lib)
endif()

#option(FIND_OPENMP_ENABLED "Enable FindOpenMP, not necessary on Apple" OFF)

option(INTEGRATION_ENABLED "Build integration" FALSE)
if(INTEGRATION_ENABLED)
  ADD_SUBDIRECTORY(integrator)
endif()

option(ONLY_INTEGRATOR_ENABLED "Build only integration" FALSE)
if(ONLY_INTEGRATOR_ENABLED)
  ADD_SUBDIRECTORY(onlyintegration)
endif()

option(ON_DEMAND_ACCURATE_ENABLED "Build ondemand accurate" TRUE)
if(ON_DEMAND_ACCURATE_ENABLED)
  ADD_SUBDIRECTORY(ondemandaccurate)
  ADD_SUBDIRECTORY(ondemandaccurate2d)
endif()

option(STREMLINES_2D_ENABLED "Build 2d streamline integrator" TRUE)
if(STREMLINES_2D_ENABLED)
  ADD_SUBDIRECTORY(streamlines2d)
endif()

option(SPARSEBLOCKEDCOMPUTE_ENABLED "Build sparse blocked compute integrator" TRUE)
if(SPARSEBLOCKEDCOMPUTE_ENABLED)
  ADD_SUBDIRECTORY(sparseblockedcompute)
endif()

option(ALL_INTEGRATION_MSC_ENABLED "Build all integration MSC" FALSE)
if(ALL_INTEGRATION_MSC_ENABLED)
  ADD_SUBDIRECTORY(allintegrationmorsecomplex)
endif()

option(CONVERGENT_ENABLED "Build convergent" FALSE)
if(CONVERGENT_ENABLED)
  ADD_SUBDIRECTORY(convergent)
endif()

option(EXTRACT_GRAPH_ENABLED "Build extract graph" FALSE)
if(EXTRACT_GRAPH_ENABLED)
  ADD_SUBDIRECTORY(extractgraph)
endif()

option(EXTRACT_MSC_ENABLED "Build extract msc" true)
if(EXTRACT_MSC_ENABLED)
  ADD_SUBDIRECTORY(extractmsc)
endif()

option(EXTRACT_CROSS_SECTION_ENABLED "Build extract cross section" FALSE)
if(EXTRACT_CROSS_SECTION_ENABLED)
  ADD_SUBDIRECTORY(extractcrosssection)
endif()

option(EXTRACT_CURVATURE_ENABLED "Build extract curvature" FALSE)
if(EXTRACT_CURVATURE_ENABLED)
  ADD_SUBDIRECTORY(extractcurvature)
endif()

option(LABEL_CONN_COMP_ENABLED "Label connected components" FALSE)
if(LABEL_CONN_COMP_ENABLED)
  ADD_SUBDIRECTORY(labelconnectedcomponents)
endif()

option(STEEPEST_ENABLED "Build steepest" FALSE)
if(STEEPEST_ENABLED)
  ADD_SUBDIRECTORY(steepest)
endif()

option(STEEPEST_LSTAR_ENABLED "Build steepest lstar" FALSE)
if(STEEPEST_LSTAR_ENABLED)
  ADD_SUBDIRECTORY(steepest_lstar)
endif()

option(MAKE_CONST_SMALL_GRID_ENABLED "Build make const small grid" FALSE)
if(MAKE_CONST_SMALL_GRID_ENABLED)
  ADD_SUBDIRECTORY(make_const_small_grid)
endif()

option(MORSE_HIERARCHY_2D_ENABLED "Build Morse hierarchy 2D" FALSE)
if(MORSE_HIERARCHY_2D_ENABLED)
  ADD_SUBDIRECTORY(morsehierarchy2d)
endif()

option(PATH_COMPRESS_INTEGRATOR_ENABLED "Build compress integrator" FALSE)
if(PATH_COMPRESS_INTEGRATOR_ENABLED)
  ADD_SUBDIRECTORY(pathcompressintegrator)
endif()

option(MIN_NETWORK_ENABLED "Build min network" FALSE)
if(MIN_NETWORK_ENABLED)
  ADD_SUBDIRECTORY(buildminnetwork)
endif()

option(TRIANGLE_TEST_ENABLED "Build triangle test" FALSE)
if(TRIANGLE_TEST_ENABLED)
  ADD_SUBDIRECTORY(triangletest)
endif()

option(GENERATE_ISOSURFACE_ENABLED "Generate isosurface util" FALSE)
if(GENERATE_ISOSURFACE_ENABLED)
  ADD_SUBDIRECTORY(generateisosurface)
endif()

option(EXPERIMENTAL_ENABLED "Build experimental" FALSE)
if(EXPERIMENTAL_ENABLED)
  ADD_SUBDIRECTORY(experimental)
endif()  
