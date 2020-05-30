paraview_add_plugin(SurfaceCut
  REQUIRED_ON_SERVER
  VERSION "1.0"
  MODULES SurfaceCutFilter
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Filter/vtk.module")
