cmake_minimum_required(VERSION 3.0.0)
project(mylbm VERSION 0.1.0)

include(CTest)
enable_testing()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# add the MathFunctions library
add_subdirectory(src)

# add the executable
add_executable(mylbm main.cpp)

target_link_libraries(mylbm PUBLIC grid problem)

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
target_include_directories(mylbm PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/src"
                          )
