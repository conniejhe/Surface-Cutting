cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)

set(CMAKE_MACOSX_RPATH 1)

if (NOT ParaView_BINARY_DIR)
  find_package(ParaView REQUIRED)
  include(${PARAVIEW_USE_FILE})
endif()

include(ParaViewPlugins)

set(BLA_VENDOR Intel10_64lp)
find_package(LAPACK)

# create a paraview plugin containing server manager xml and the server
# manager classes to build
# this plugin can be loaded on the server side

ADD_PARAVIEW_PLUGIN(SurfaceTrackerText "1.0"
  SERVER_MANAGER_XML vtkDijkstraGraphGeodesicPath1.xml
  SERVER_MANAGER_SOURCES vtkDijkstraGraphGeodesicPath1.cxx
  ../Curvature/vtkCurvatures1.cxx
  ../Curvature/vtkCurvatures1.h)

  target_link_libraries(SurfaceTrackerText PRIVATE
        /usr/local/opt/lapack/lib/libblas.dylib
        /usr/local/opt/lapack/lib//liblapack.dylib
        -O2 BLAS LAPACK
  )
