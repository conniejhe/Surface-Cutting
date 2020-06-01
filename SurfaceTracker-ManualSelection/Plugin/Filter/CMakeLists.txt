set(classes
  vtkMyDGGP vtkCurvatures1 vtkDijkstraGraphGeodesicPath1)

vtk_module_add_module(SurfaceTrackerManualFilter
  CLASSES ${classes})

paraview_add_server_manager_xmls(
  XMLS  vtkMyDGGP.xml)

# Change to correct paths of blas and lapack libraries
# .so or .dylib file extensions supported
target_link_libraries(SurfaceTrackerManualFilter PRIVATE
    /usr/local/opt/lapack/lib/libblas.dylib
    /usr/local/opt/lapack/lib//liblapack.dylib)
