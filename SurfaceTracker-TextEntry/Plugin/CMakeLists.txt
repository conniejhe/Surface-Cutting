paraview_add_plugin(SurfaceTrackerText
  REQUIRED_ON_SERVER
  VERSION "1.0"
  MODULES SurfaceTrackerTextFilter
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Filter/vtk.module")

  target_link_libraries(SurfaceTrackerText PRIVATE
      /usr/local/opt/lapack/lib/libblas.dylib
      /usr/local/opt/lapack/lib//liblapack.dylib
      -O2 BLAS LAPACK
  )
