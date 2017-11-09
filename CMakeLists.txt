cmake_minimum_required(VERSION 3.1)

project("dune-spgrid" C CXX)

find_package(dune-common)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules")

include(DuneMacros)
dune_project()

add_subdirectory(cmake/modules)
add_subdirectory(doc)
add_subdirectory(dune)

if(dune-python_FOUND)
  add_subdirectory(python)
  dune_python_install_package(PATH "python")
endif()

finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
