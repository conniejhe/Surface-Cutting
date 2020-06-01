paraview_add_plugin(myCurvature
  REQUIRED_ON_SERVER
  VERSION "1.0"
  MODULES myCurvatureFilter
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Filter/vtk.module")
