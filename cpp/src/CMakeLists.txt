set(YEHU_SRC
    yehu.h )

add_library(yehu SHARED ${YEHU_SRC})

target_include_directories(yehu
    PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR})

if(APPLE)
    target_link_libraries(yehu PUBLIC Eigen3::Eigen "-framework Accelerate")

elseif(UNIX)
    find_package(BLAS QUIET)
    find_package(LAPACK QUIET)
    if(BLAS_FOUND)
        message(STATUS "BLAS found — enabling Eigen BLAS backend")
        target_compile_definitions(yehu PUBLIC EIGEN_USE_BLAS)
        target_link_libraries(yehu PUBLIC Eigen3::Eigen ${BLAS_LIBRARIES})
        if(LAPACK_FOUND)
            target_compile_definitions(yehu PUBLIC EIGEN_USE_LAPACKE)
            target_link_libraries(yehu PUBLIC ${LAPACK_LIBRARIES})
        endif()
    else()
        message(STATUS "BLAS not found — Eigen will use its built-in implementations")
        target_link_libraries(yehu PUBLIC Eigen3::Eigen)
    endif()

elseif(WIN32)
    find_package(BLAS QUIET)
    find_package(LAPACK QUIET)
    if(BLAS_FOUND)
        message(STATUS "BLAS found (${BLA_VENDOR}) — enabling Eigen BLAS backend")
        target_compile_definitions(yehu PUBLIC EIGEN_USE_BLAS)
        target_link_libraries(yehu PUBLIC Eigen3::Eigen ${BLAS_LIBRARIES})
        if(LAPACK_FOUND)
            target_compile_definitions(yehu PUBLIC EIGEN_USE_LAPACKE)
            target_link_libraries(yehu PUBLIC ${LAPACK_LIBRARIES})
        endif()
    else()
        message(STATUS "No BLAS found on Windows — using Eigen built-in implementations")
        target_link_libraries(yehu PUBLIC Eigen3::Eigen)
    endif()
endif()

set_target_properties(yehu PROPERTIES LINKER_LANGUAGE CXX)
