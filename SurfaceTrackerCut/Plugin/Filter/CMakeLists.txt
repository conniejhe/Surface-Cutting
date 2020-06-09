set(classes
  vtkCurvatures1 vtkDijkstraGraphGeodesicPath1 vtkMyDGGP surfaceCut)

vtk_module_add_module(SurfaceTrackerCutPlugin
  CLASSES ${classes})

paraview_add_server_manager_xmls(
  XMLS  vtkCurvatures1.xml vtkDijkstraGraphGeodesicPath1.xml vtkMyDGGP.xml surfaceCut.xml)

# Change to correct paths of blas and lapack libraries
# .so or .dylib file extensions supported

target_link_libraries(SurfaceTrackerCutPlugin PRIVATE
  /usr/local/opt/lapack/lib/libblas.dylib
  /usr/local/opt/lapack/lib//liblapack.dylib
)
