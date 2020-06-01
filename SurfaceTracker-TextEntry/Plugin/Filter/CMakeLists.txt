set(classes
  vtkDijkstraGraphGeodesicPath1
  vtkCurvatures1)

vtk_module_add_module(SurfaceTrackerTextFilter
  CLASSES ${classes})

paraview_add_server_manager_xmls(
  XMLS  vtkDijkstraGraphGeodesicPath1.xml)

target_link_libraries(SurfaceTrackerTextFilter PRIVATE
    /usr/local/opt/lapack/lib/libblas.dylib
    /usr/local/opt/lapack/lib//liblapack.dylib
)
