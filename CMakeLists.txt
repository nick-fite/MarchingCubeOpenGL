cmake_minimum_required(VERSION 3.10)

project(ShopliftingBoy)

set(CMAKE_CXX_STANDARD 17)
if(MSVC)
  add_compile_options(/EHsc /MP)
else()
  add_compile_options(-fexceptions)
endif()
# Ensure vcpkg integration
set(CMAKE_TOOLCHAIN_FILE "C:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")




find_package(glfw3 REQUIRED)
find_package(GLEW REQUIRED)
find_package(Freetype REQUIRED)
find_package(assimp REQUIRED)
find_package(glm REQUIRED)
find_package(freeimage REQUIRED)
find_package(OpenGL REQUIRED)

set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/src")

file(GLOB_RECURSE SRC_FILES 
    "${SRC}/*.cpp"
    "${SRC}/*.h"
    "${SRC}/*.hpp"
    "${SRC}/*.c"
)

add_executable(ShopliftingBoy main.cpp ${SRC_FILES})

target_link_libraries(ShopliftingBoy PRIVATE
    glfw
    GLEW::GLEW
    Freetype::Freetype
    assimp::assimp
    glm::glm
    freeimage::FreeImage
    freeimage::FreeImagePlus
    OpenGL::GL
)