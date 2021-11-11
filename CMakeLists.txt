cmake_minimum_required(VERSION 3.10.0)
project(physics)


# set(CMAKE_BUILD_TYPE RELEASE)

if (POLICY CMP0072)
  set(OpenGL_GL_PREFERENCE LEGACY)
endif()

# find the required packages
find_package(Eigen3 REQUIRED)
if(Eigen3_FOUND)
    message(STATUS "EIGEN3 included at ${EIGEN3_INCLUDE_DIR}")
    include_directories(${EIGEN3_INCLUDE_DIR})
endif()

include_directories(./)
add_subdirectory(core)

# add_executable(rigidbody rigidbody.h rigidbody.cpp jsoncpp.cpp)
# add_executable(explicitfem explicitfem.h explicitfem.cpp jsoncpp.cpp)
# add_executable(implicitfem implicitfem.h implicitfem.cpp jsoncpp.cpp)
# add_executable(fluid fluid.h fluid.cpp jsoncpp.cpp)
# add_executable(springmass springmass.h springmass.cpp springmass.cpp jsoncpp.cpp)
# add_executable(partviewer partviewer.cpp)
# target_link_libraries(partviewer glut GL GLU) 


# add_executable(delta "delta/test.cpp")
# add_executable(svd_dd "delta/svd_dd.cpp")
# target_link_libraries(svd_dd Eigen3::Eigen)

add_executable(main "main.cpp")
target_link_libraries(main Eigen3::Eigen)

add_executable(test_lapack "test_lapack.cpp")
target_link_libraries(test_lapack lapack blas)


add_executable(cons_law "constitutive_law.cpp")
target_link_libraries(cons_law autodiff Eigen3::Eigen)

add_executable(fem simplefem.cpp jsoncpp.cpp)
target_link_libraries(fem magic)


add_executable(fluid "fluid.cpp")
target_link_libraries(fluid Eigen3::Eigen)


add_executable(mpm "mpm.cpp")
target_link_libraries(mpm Eigen3::Eigen)

add_executable(clothsim "clothsim.cpp")
target_link_libraries(clothsim magic)

add_executable(clothbw clothbw.cpp jsoncpp.cpp)
target_link_libraries(clothbw magic)

list(APPEND CMAKE_MODULE_PATH external/Bartels/cmake)
include_directories(bartels_INCLUDE_DIR)