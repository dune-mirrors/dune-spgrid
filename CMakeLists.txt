cmake_minimum_required(VERSION 3.16)

project("dune-spgrid" C CXX)

find_package(dune-common)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/modules")

include(DuneMacros)

# deactivate global include-directories
dune_policy(SET DP_DEFAULT_INCLUDE_DIRS NEW)

# deactivate global calls to add_dune_all_flags in tests
dune_policy(SET DP_TEST_ADD_ALL_FLAGS NEW)

# do not make suggested dependencies required automatically
dune_policy(SET DP_SUGGESTED_MODULE_DEPENDENCIES_REQUIRED_DOWNSTREAM NEW)

# start a dune project with information from dune.module
dune_project()

dune_add_library(dunespgrid INTERFACE
  EXPORT_NAME SPGrid
  LINK_LIBRARIES Dune::Grid)

dune_default_include_directories(dunespgrid INTERFACE)
link_libraries(Dune::SPGrid)

#target_link_libraries(dunespgrid INTERFACE Dune::Grid)


add_subdirectory(cmake/modules)
add_subdirectory(doc)
add_subdirectory(dune)

# if Python bindings are enabled, include necessary sub directories.
if( DUNE_ENABLE_PYTHONBINDINGS )
  add_subdirectory(python)
endif()

finalize_dune_project()
