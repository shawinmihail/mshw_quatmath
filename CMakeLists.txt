cmake_minimum_required(VERSION 3.1)

project(mshw_quatmath)

# https://eigen.tuxfamily.org/dox/TopicCMakeGuide.html
find_package (Eigen3 REQUIRED NO_MODULE)

include_directories(include ${EIGEN3_INCLUDE_DIRS})

add_executable(mshw_quatmath src/mshw_quatmath_test.cpp)
