set(classes
  vtkCurvatures1 vtkDijkstraGraphGeodesicPath1 vtkMyDGGP surfaceCut)

vtk_module_add_module(SurfaceTrackerCutPlugin
  CLASSES ${classes})

paraview_add_server_manager_xmls(
  XMLS  vtkCurvatures1.xml vtkDijkstraGraphGeodesicPath1.xml vtkMyDGGP.xml surfaceCut.xml)

find_package(lapack REQUIRED)

target_link_libraries(SurfaceTrackerCutPlugin PRIVATE lapack)
