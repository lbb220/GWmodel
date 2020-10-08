cmake_minimum_required(VERSION 3.9)
project(GWmodelCUDA LANGUAGES CXX CUDA)

if(WIN32 AND NOT MSVC)
    message(FATAL_ERROR "CUDA only support MSVC complier on Windows!")
endif(WIN32 AND NOT MSVC)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(cuBLAS REQUIRED)
include_directories(${CUBLAS_INCLUDE_DIRS})

add_definitions(-DCREATDLL_EXPORTS)

set(SOURCES_CXX
    CGWmodelCUDA.cpp
    IGWmodelCUDA.cpp
)

set(SOURCES_CUDA
    GWmodelKernel.cu
)

set(HEADERS
    CGWmodelCUDA.h
    GWmodelKernel.h
    helper.h
    IGWmodelCUDA.h
)


add_library(GWmodelCUDA64 SHARED ${SOURCES_CUDA} ${SOURCES_CXX} ${HEADERS})
set_property(TARGET GWmodelCUDA64 PROPERTY CUDA_ARCHITECTURES 35 37 52 60 61 70 75)
target_link_libraries(GWmodelCUDA64 PRIVATE 
    ${CUDA_LIBRARIES}
    ${CUBLAS_LIBRARIES}
    ${ARMADILLO_LIBRARIES}
)
target_compile_features(GWmodelCUDA64 PUBLIC cxx_std_11)
add_custom_command(TARGET GWmodelCUDA64 POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:GWmodelCUDA64> "${CMAKE_CURRENT_SOURCE_DIR}/.."
)
add_custom_command(TARGET GWmodelCUDA64 POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_LINKER_FILE:GWmodelCUDA64> "${CMAKE_CURRENT_SOURCE_DIR}/.."
)
add_custom_command(TARGET GWmodelCUDA64 POST_BUILD 
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_LINKER_FILE:GWmodelCUDA64> "${CMAKE_CURRENT_SOURCE_DIR}"
)