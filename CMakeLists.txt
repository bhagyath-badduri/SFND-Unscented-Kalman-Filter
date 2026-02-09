cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(playback)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# Compiler warnings

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")


# Eigen fixes (IMPORTANT)
# Prevent runtime alignment crash

add_definitions(
  -DEIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
  -DEIGEN_DONT_VECTORIZE
)


# PCL

find_package(PCL REQUIRED)

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

# Some PCL builds include vtkproj4 — remove if present
list(REMOVE_ITEM PCL_LIBRARIES "vtkproj4")

# Build target

add_executable(ukf_highway
  src/main.cpp
  src/ukf.cpp
  src/tools.cpp
  src/render/render.cpp
)

target_link_libraries(ukf_highway ${PCL_LIBRARIES})
