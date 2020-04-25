set(classes
  surfaceCut)

vtk_module_add_module(surfaceCut::surfaceCutFilter
  CLASSES ${classes})

paraview_add_server_manager_xmls(
  XMLS  surfaceCut.xml)
