set(classes
  vtkCurvatures1)

vtk_module_add_module(myCurvature::myCurvatureFilter
  CLASSES ${classes})

paraview_add_server_manager_xmls(
  XMLS  vtkCurvatures1.xml)
