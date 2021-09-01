enable_language(Fortran)
add_custom_command(OUTPUT elasticity.f90
    COMMAND "maxima" ARGS  --userdir=${CMAKE_CURRENT_SOURCE_DIR} -b elasticity.mac 
    DEPENDS elasticity.mac
)

add_library(autodiff elasticity.f90)