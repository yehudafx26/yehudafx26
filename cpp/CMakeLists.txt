cmake_minimum_required(VERSION 3.25)

project(yehucpp)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# Get EIGEN
include(FetchContent)
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(eigen)

# Set compliation flags
if(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /arch:AVX2 /GL")
else()
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -funroll-all-loops")
endif()

add_subdirectory(src)
add_subdirectory(mex)
