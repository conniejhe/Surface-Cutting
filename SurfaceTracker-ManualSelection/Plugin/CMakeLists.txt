paraview_add_plugin(SurfaceTrackerManual
  REQUIRED_ON_SERVER
  VERSION "1.0"
  MODULES SurfaceTrackerManualFilter
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Filter/vtk.module")

  target_link_libraries(SurfaceTrackerManual PRIVATE
      /usr/local/opt/lapack/lib/libblas.dylib
      /usr/local/opt/lapack/lib//liblapack.dylib
      -O2 BLAS LAPACK
  )
