set(classes
  vtkDijkstraGraphGeodesicPath1)

vtk_module_add_module(SurfaceTrackerText::SurfaceTrackerTextFilter
  CLASSES ${classes})

paraview_add_server_manager_xmls(
  XMLS  vtkDijkstraGraphGeodesicPath1.xml)
