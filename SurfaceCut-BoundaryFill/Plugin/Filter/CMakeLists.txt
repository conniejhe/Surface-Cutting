set(classes
  surfaceCut)

vtk_module_add_module(SurfaceCutFilter
  CLASSES ${classes})

paraview_add_server_manager_xmls(
  XMLS  surfaceCut.xml)
