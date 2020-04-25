cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)

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

ADD_PARAVIEW_PLUGIN(SurfaceTrackerManual "1.0"
  SERVER_MANAGER_XML vtkMyDGGP.xml
  SERVER_MANAGER_SOURCES vtkMyDGGP.cxx
  ../SurfaceTracker-TextEntry/vtkDijkstraGraphGeodesicPath1.cxx
  ../Curvature/vtkCurvatures1.cxx)

  TARGET_LINK_LIBRARIES(SurfaceTrackerManual
    ${DGGP1_LIBRARY}
    /usr/lib/liblapack.so
    /usr/lib/libblas.so
  )
