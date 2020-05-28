set(classes
  vtkCurvatures1)

vtk_module_add_module(myCurvature::myCurvatureFilter
  CLASSES ${classes})

paraview_add_server_manager_xmls(
  XMLS  vtkCurvatures1.xml)

target_link_libraries(myCurvatureFilter PRIVATE
  /usr/local/opt/lapack/lib/libblas.dylib
  /usr/local/opt/lapack/lib//liblapack.dylib
  -O2 BLAS LAPACK
)
