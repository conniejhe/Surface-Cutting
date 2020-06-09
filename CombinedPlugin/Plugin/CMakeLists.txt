paraview_add_plugin(CombinedFilter
  REQUIRED_ON_SERVER
  VERSION "1.0"
  MODULES EverythingFilter
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Filter/vtk.module")
