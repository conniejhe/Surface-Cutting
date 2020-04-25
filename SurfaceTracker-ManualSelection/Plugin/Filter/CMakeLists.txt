set(classes
  vtkMyDGGP)

vtk_module_add_module(SurfaceTrackerManual::SurfaceTrackerManualFilter
  CLASSES ${classes})

paraview_add_server_manager_xmls(
  XMLS  vtkMyDGGP.xml)
