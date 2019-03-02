cmake_minimum_required(VERSION 3.5)

project (fvm)

set (CMAKE_CXX_STANDARD 14)

file (GLOB sources src/*.cpp src/*.h)

add_executable (fvm ${sources})

########
# Libs #
########

# tecio

add_library (tecio INTERFACE IMPORTED GLOBAL)

set_property (
    TARGET tecio PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES
    ${CMAKE_SOURCE_DIR}/lib/tecio/include/
)

set_property (
    TARGET tecio PROPERTY
    INTERFACE_LINK_LIBRARIES
    ${CMAKE_SOURCE_DIR}/lib/tecio/lib/libtecio.a
)

target_link_libraries(
    fvm
        tecio
)

