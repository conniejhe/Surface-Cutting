paraview_add_plugin(combinedPlugin
    REQUIRED_ON_SERVER
    VERSION "1.0"
    MODULES myCurvatureFilter SurfaceCutFilter SurfaceTrackerManualFilter SurfaceTrackerTextFilter
    MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Plugin/myCurvatureFilter/Plugin/Filter/vtk.module"
    "${CMAKE_CURRENT_SOURCE_DIR}/Plugin/SurfaceCutFilter/Plugin/Filter/vtk.module"
    "${CMAKE_CURRENT_SOURCE_DIR}/Plugin/SurfaceTrackerManualFilter/Plugin/Filter/vtk.module"
    "${CMAKE_CURRENT_SOURCE_DIR}/Plugin/SurfaceTrackerTextFilter/Plugin/Filter/vtk.module"
)
