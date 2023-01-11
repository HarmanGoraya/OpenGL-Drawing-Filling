# welcome to ECS 175
cmake_minimum_required(VERSION 3.0)
project(ecs175-template)

# enable rpath on macos (advanced configuration)
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_MACOSX_RPATH 1)
endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
option(OFFSCREEN_RENDERING "Enable offline rendering (only for TA)" OFF)

# --------------------------------------
# load/build required libraries
#
# ---> OpenGL (advanced configuration)
# reference https://cmake.org/cmake/help/v3.18/module/FindOpenGL.html
set(OpenGL_GL_PREFERENCE GLVND)
find_package(OpenGL REQUIRED)
# correctly select OpenGL targets (modern CMake)
if(TARGET OpenGL::OpenGL)
  list(APPEND gl_libs OpenGL::OpenGL)
else()
  list(APPEND gl_libs OpenGL::GL)
endif()
if(TARGET OpenGL::GLU)
  list(APPEND gl_libs OpenGL::GLU)
endif()

# ---> build glfw from source
if (NOT OFFSCREEN_RENDERING)
  # use GLFW for interactive GL rendering 
  # as it is cross platform compatable
  if(TARGET OpenGL::GLX)
    list(APPEND gl_libs OpenGL::GLX)
  endif()
  set(GLFW_USE_OSMESA OFF)
else()
  # enable software rendering via OSMesa
  set(GLFW_USE_OSMESA ON)
endif()
# configure other options
set(GLFW_BUILD_DOCS OFF)
set(GLFW_BUILD_EXAMPLES OFF)
set(GLFW_BUILD_TESTS OFF)
set(GLFW_INSTALL OFF)
add_subdirectory(${PROJECT_SOURCE_DIR}/glfw-3.3.2)
# hide these options by default
mark_as_advanced(GLFW_INSTALL 
  GLFW_BUILD_DOCS GLFW_BUILD_TESTS GLFW_BUILD_EXAMPLES
  GLFW_USE_OSMESA GLFW_USE_WAYLAND GLFW_VULKAN_STATIC)
list(APPEND gl_libs glfw)

# ---> build glad from source
add_subdirectory(glad-compatible)
add_subdirectory(glad-core-3.3)
add_subdirectory(glad-core-4.1)
set(ALL_GL_VERSIONS compatible core-3.3 core-4.1)

# ---> build "Dear ImGui"
add_subdirectory(imgui)

# ---> build the project utility library
add_subdirectory(util)

# build sub-projects
add_subdirectory(demo)
add_subdirectory(projects/p1)
