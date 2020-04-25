paraview_add_plugin(surfaceCut
  REQUIRED_ON_SERVER
  VERSION "1.0"
  MODULES surfaceCut::surfaceCutFilter
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Filter/vtk.module")
