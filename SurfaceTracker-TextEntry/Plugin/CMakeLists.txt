paraview_add_plugin(SurfaceTrackerText
  REQUIRED_ON_SERVER
  VERSION "1.0"
  MODULES SurfaceTrackerTextFilter
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Filter/vtk.module")
