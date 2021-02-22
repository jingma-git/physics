cmake_minimum_required(VERSION 3.10.0)
project(physics)

set(CMAKE_BUILD_TYPE RELEASE)
if (POLICY CMP0072)
  set(OpenGL_GL_PREFERENCE LEGACY)
endif()

# find the required packages
find_package(Eigen3 REQUIRED)
if(Eigen3_FOUND)
    message(STATUS "EIGEN3 included at ${EIGEN3_INCLUDE_DIR}")
    include_directories(${EIGEN3_INCLUDE_DIR})
endif()

# add_executable(rigidbody rigidbody.h rigidbody.cpp jsoncpp.cpp)
# add_executable(explicitfem explicitfem.h explicitfem.cpp jsoncpp.cpp)
add_executable(implicitfem implicitfem.h implicitfem.cpp jsoncpp.cpp)
# add_executable(fluid fluid.h fluid.cpp jsoncpp.cpp)
# add_executable(springmass springmass.h springmass.cpp springmass.cpp jsoncpp.cpp)
# add_executable(partviewer partviewer.cpp)
# target_link_libraries(partviewer glut GL GLU) 


