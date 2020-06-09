paraview_add_plugin(SurfaceTrackerCut
  REQUIRED_ON_SERVER
  VERSION "1.0"
  MODULES SurfaceTrackerCutPlugin
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Filter/vtk.module")
