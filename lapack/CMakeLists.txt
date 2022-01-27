file(GLOB LAPACK_SRC "*.f")
add_library(lapack_zheev STATIC ${LAPACK_SRC})
