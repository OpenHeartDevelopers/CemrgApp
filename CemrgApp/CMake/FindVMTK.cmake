# FindVMTK
# ---------
#
# Find VMTK libraries and applications
#
# The module defines the following variables::
#
#  VMTK_INCLUDE_DIRS  - Directories to include to use VMTK
#  VMTK_LIBRARIES     - Files to link against to use VMTK
#  VMTK_FOUND         - If false, don't try to use VMTK
#  VMTK_BUILD_DIR           - (optional) Source directory for VMTK
set(_vmtk_dir_description "The directory of VMTK build or install tree")
set(VMTK_BUILD_DIR "${MITK_EXTERNAL_PROJECT_PREFIX}/src/VMTK-build/VMTK-Build/"
CACHE PATH ${_vmtk_dir_description} FORCE)

mark_as_advanced(VMTK_BUILD_DIR)

set(_SAVED_VMTK_DIR ${VMTK_BUILD_DIR})
message("VMTK_BUILD_DIR: ${VMTK_BUILD_DIR}")

# Step1: Attempt to find a version of VMTK providing a VMTKConfig.cmake file.
if(NOT VMTK_FIND_QUIETLY)
  message(STATUS "Trying to find VMTK expecting VMTKConfig.cmake")
endif()
find_package(VMTK QUIET CONFIG
  PATHS ${VMTK_BUILD_DIR}
)

if(VMTK_FOUND)
  #if(NOT VMTK_FIND_QUIETLY)
  message(STATUS "Trying to find VMTK expecting VMTKConfig.cmake - ok")

  INCLUDE_DIRECTORIES(${VMTK_INCLUDE_DIRS})
  LINK_DIRECTORIES(${VMTK_LIBRARY_DIRS})
  LINK_LIBRARIES(vtkvmtkCommon vtkvmtkComputationalGeometry
    vtkvmtkContrib vtkvmtkDifferentialGeometry vtkvmtkIO
    vtkvmtkITK vtkvmtkMisc vtkvmtkSegmentation nl tet)
  #endif()
  
  return()
else()
  if(NOT VMTK_FIND_QUIETLY)
    message(STATUS "Trying to find VMTK expecting VMTKConfig.cmake - failed")
    message(STATUS "Trying to find VMTK expecting FindVMTK.cmake.")

    set(VMTK_BUILD_DIR ${_SAVED_VMTK_DIR} CACHE PATH ${_VMTK_dir_description} FORCE)

  endif()
endif()

#

# # Custom find script for VMTK
# include(FindPackageHandleStandardArgs)
# include(SelectLibraryConfigurations)
#
# set(VMTK_BUILD_DIR "${MITK_EXTERNAL_PROJECT_PREFIX}/src/VMTK-build/VMTK-Build/")
#
# find_path(VMTK_INCLUDE_DIR
#   NAMES vtkvmtkConfigure.h
#   PATHS ${VMTK_BUILD_DIR}
#   # PATH_SUFFIXES include/vmtk
#   PATH_SUFFIXES vtkVmtk
# )
#
# set(VMTK_LIBRARY_DIR "${VMTK_BUILD_DIR}/bin")
#
# foreach(lib Common ComputationalGeometry DifferentialGeometry IO ITK Misc Segmentation)
#   string(TOUPPER ${lib} LIB)
#   if(NOT VMTK_${LIB}_LIBRARIES)
#     find_library(VMTK_${LIB}_LIBRARY_RELEASE vtkvmtk${lib} ${VMTK_LIBRARY_DIR})
#     find_library(VMTK_${LIB}_LIBRARY_DEBUG vtkvmtk${lib}d ${VMTK_LIBRARY_DIR})
#     select_library_configurations(VMTK_${LIB})
#   endif()
# endforeach()
#
# find_package_handle_standard_args(VMTK DEFAULT_MSG
#   VMTK_INCLUDE_DIR
#   VMTK_COMMON_LIBRARIES
#   VMTK_COMPUTATIONALGEOMETRY_LIBRARIES
#   VMTK_DIFFERENTIALGEOMETRY_LIBRARIES
#   VMTK_IO_LIBRARIES
#   VMTK_ITK_LIBRARIES
#   VMTK_MISC_LIBRARIES
#   VMTK_SEGMENTATION_LIBRARIES
# )
#
# if(VMTK_FOUND)
#   set(VMTK_INCLUDE_DIRS ${VMTK_INCLUDE_DIR})
#   set(VMTK_LIBRARIES
#     ${VMTK_COMMON_LIBRARIES}
#     ${VMTK_COMPUTATIONALGEOMETRY_LIBRARIES}
#     ${VMTK_DIFFERENTIALGEOMETRY_LIBRARIES}
#     ${VMTK_IO_LIBRARIES}
#     ${VMTK_ITK_LIBRARIES}
#     ${VMTK_MISC_LIBRARIES}
#     ${VMTK_SEGMENTATION_LIBRARIES}
#   )
# endif()
#
# mark_as_advanced(
#   VMTK_INCLUDE_DIR
# )
