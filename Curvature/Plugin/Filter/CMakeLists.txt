set(classes
  vtkCurvatures1)

vtk_module_add_module(myCurvatureFilter
  CLASSES ${classes})

paraview_add_server_manager_xmls(
  XMLS  vtkCurvatures1.xml)

# Change to correct paths of blas and lapack libraries
# .so or .dylib file extensions supported

target_link_libraries(myCurvatureFilter PRIVATE
  /usr/local/opt/lapack/lib/libblas.dylib
  /usr/local/opt/lapack/lib//liblapack.dylib
)
