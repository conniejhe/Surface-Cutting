set(classes
  vtkCurvatures1 vtkDijkstraGraphGeodesicPath1 vtkMyDGGP surfaceCut)

vtk_module_add_module(SurfaceTrackerCutPlugin
  CLASSES ${classes})

paraview_add_server_manager_xmls(
  XMLS  vtkCurvatures1.xml vtkDijkstraGraphGeodesicPath1.xml vtkMyDGGP.xml surfaceCut.xml)

find_package(LAPACK REQUIRED)
vtk_module_link(SurfaceTrackerCutPlugin
  PRIVATE ${LAPACK_LIBRARIES})
if (LAPACK_LINKER_FLAGS)
  vtk_module_link_options(SurfaceTrackerCutPlugin
    PRIVATE ${LAPACK_LINKER_FLAGS})
endif()
