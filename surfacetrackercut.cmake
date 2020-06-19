superbuild_add_project(surfacetrackercut
    DEPENDS paraview lapack
    DEPENDS_OPTIONAL qt5
    CMAKE_ARGS
        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS})
