
list(APPEND CMAKE_CXX_FLAGS "-fopenmp -std=c++11 -fpermissive")

# find the required packages
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    message(STATUS "Found OpenMP")
endif()

find_package(Eigen3 REQUIRED)
if(Eigen3_FOUND)
    message(STATUS "EIGEN3 included at ${EIGEN3_INCLUDE_DIR}")
    include_directories(${EIGEN3_INCLUDE_DIR})
endif()

find_package(Boost COMPONENTS filesystem REQUIRED)
if(Boost_FOUND)
    message("-- Boost_INCLUDE ${Boost_INCLUDE_DIRS}")
    message("-- Boost_LIB ${Boost_LIBRARIES}")
endif()

add_subdirectory(autodiff)

file(GLOB src "*.cpp")
add_library(magic ${src})
target_link_libraries(magic autodiff Eigen3::Eigen OpenMP::OpenMP_CXX ${Boost_LIBRARIES})