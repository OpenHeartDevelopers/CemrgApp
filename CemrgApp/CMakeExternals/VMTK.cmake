# ----------
# VMTK
# ----------
set(proj VMTK)
set(proj_DEPENDENCIES ITK VTK)
set(${proj}_DEPENDS ${proj})

if(MITK_USE_VMTK)
  # 'Sanity' checks
  if(DEFINED ${proj}_DIR AND NOT EXISTS ${${proj}_DIR})
    message(FATAL_ERROR "${proj}_DIR variable is defined but corresponds to non-existing directory")
  endif()

  if(NOT DEFINED ${proj}_DIR)
    set(additional_cmake_args)

    if(CTEST_USE_LAUNCHERS)
      list(APPEND additional_cmake_args
      -DCMAKE_PROJECT_VTK_VMTK_INCLUDE:FILEPATH=${CMAKE_ROOT}/Modules/CTestUseLaunchers.cmake
    )
    endif()

    if(MITK_USE_Python)
      list(APPEND additional_cmake_args
        -DVTK_VMTK_WRAP_PYTHON:BOOL=ON
        -DPYTHON_EXECUTABLE:FILEPATH=${PYTHON_EXECUTABLE}
        -DPYTHON_INCLUDE_DIR:PATH=${PYTHON_INCLUDE_DIR}
        -DPYTHON_LIBRARY:FILEPATH=${PYTHON_LIBRARY}
        -DVTK_PYTHON_VERSION:STRING=3
      )
    else()
      list(APPEND additional_cmake_args
        -DVTK_VMTK_WRAP_PYTHON:BOOL=OFF
      )
    endif()

    set(VMTK_GIT_REPOSITORY "https://github.com/vmtk/vmtk.git" CACHE STRING "The git repository for cloning VMTK")
    set(VMTK_GIT_TAG "v1.4.0" CACHE STRING "The git tag/hash to be used when cloning from VMTK_GIT_REPOSITORY")
    mark_as_advanced(VMTK_GIT_REPOSITORY VMTK_GIT_TAG)

    if (UNIX)
      ExternalProject_Add(${proj}
        PREFIX ${ep_prefix}
        GIT_REPOSITORY ${VMTK_GIT_REPOSITORY}
        GIT_TAG ${VMTK_GIT_TAG}
        GIT_SHALLOW TRUE
        CMAKE_ARGS
        ${ep_common_args}
        -DSUPERBUILD_INSTALL_PREFIX:PATH=${ep_prefix}
        -DUSE_SYSTEM_ITK:BOOL=ON
        -DITK_DIR:PATH=${ep_prefix}/src/ITK-build
        -DUSE_SYSTEM_VTK:BOOL=ON
        -DVTK_DIR:PATH=${ep_prefix}/src/VTK-build
        ${additional_cmake_args}
        CMAKE_CACHE_ARGS ${ep_common_cache_args}
        -DGSL_CXX_STANDARD:STRING=${MITK_CXX_STANDARD}
        -DGSL_TEST:BOOL=OFF
        CMAKE_CACHE_DEFAULT_ARGS ${ep_common_cache_default_args}
        BUILD_COMMAND "make"
        INSTALL_COMMAND "make"
        DEPENDS ${proj_DEPENDENCIES}
      )
    else(UNIX)
      ExternalProject_Add(${proj}
        PREFIX ${ep_prefix}
        GIT_REPOSITORY ${VMTK_GIT_REPOSITORY}
        GIT_TAG ${VMTK_GIT_TAG}
        GIT_SHALLOW TRUE
        CMAKE_ARGS
        ${ep_common_args}
        -DSUPERBUILD_INSTALL_PREFIX:PATH=${ep_prefix}
        -DUSE_SYSTEM_ITK:BOOL=ON
        -DITK_DIR:PATH=${ep_prefix}/src/ITK-build
        -DUSE_SYSTEM_VTK:BOOL=ON
        -DVTK_DIR:PATH=${ep_prefix}/src/VTK-build
        ${additional_cmake_args}
        CMAKE_CACHE_ARGS ${ep_common_cache_args}
        -DGSL_CXX_STANDARD:STRING=${MITK_CXX_STANDARD}
        -DGSL_TEST:BOOL=OFF
        CMAKE_CACHE_DEFAULT_ARGS ${ep_common_cache_default_args}
        DEPENDS ${proj_DEPENDENCIES}
      )
    endif (UNIX)


    set(${proj}_DIR ${ep_prefix})
  else()
    mitkMacroEmptyExternalProject(${proj} "${proj_DEPENDENCIES}")
  endif()
endif()
