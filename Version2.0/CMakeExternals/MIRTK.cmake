# ----------
# MIRTK. NOT supported at the moment
# ----------
# Please note. This external project might be included in a following version 
set(proj MIRTK)
set(proj_DEPENDENCIES VTK)
set(${proj}_DEPENDS ${proj})

if(MITK_USE_MIRTK)
  # 'Sanity' checks
  if(DEFINED ${proj}_DIR AND NOT EXISTS ${${proj}_DIR})
    message(FATAL_ERROR "${proj}_DIR variable is defined but corresponds to non-existing directory")
  endif()

  if(NOT DEFINED ${proj}_DIR)
    set(additional_cmake_args

    # if(CTEST_USE_LAUNCHERS)
    #   list(APPEND additional_cmake_args
    #   -DCMAKE_PROJECT_VTK_VMTK_INCLUDE:FILEPATH=${CMAKE_ROOT}/Modules/CTestUseLaunchers.cmake
    # )
    # endif()

    if(MITK_USE_Python)
       list(APPEND additional_cmake_args
        -DDEPENDS_Python_DIR:FILEPATH=${PYTHON_EXECUTABLE}
       )
    # else()
    #   list(APPEND additional_cmake_args
    #     -DVTK_VMTK_WRAP_PYTHON:BOOL=OFF
    #   )
     endif()
    # git clone --depth 1 -- https://github.com/BioMedIA/MIRTK.git
    set(MIRTK_GIT_REPOSITORY "https://github.com/BioMedIA/MIRTK.git" CACHE STRING "The git repository for cloning VMTK")

    ExternalProject_Add(${proj}
      SOURCE_DIR ${CMAKE_BINARY_DIR}/${ep_prefix}
      GIT_REPOSITORY ${MIRTK_GIT_REPOSITORY}
      GIT_SHALLOW ON
      CMAKE_ARGS
          ${ep_common_args}
          -DDEPENDS_Boost_DIR:PATH=${ep_prefix}/src/Boost
          -DDEPENDS_Eigen3_DIR:PATH=${ep_prefix}/src/Eigen
          -DWITH_VTK:BOOL=ON
          -DDEPENDS_VTK:PATH=${ep_prefix}/src/VTK-build
          ${additional_cmake_args}
      CMAKE_CACHE_ARGS ${ep_common_cache_args}
        -DGSL_CXX_STANDARD:STRING=${MITK_CXX_STANDARD}
        -DGSL_TEST:BOOL=OFF
      CMAKE_CACHE_DEFAULT_ARGS ${ep_common_cache_default_args}
      BUILD_COMMAND "make"
      INSTALL_COMMAND "make install"
    )

    set(${proj}_DIR ${ep_prefix})
  else()
    mitkMacroEmptyExternalProject(${proj} "${proj_DEPENDENCIES}")
  endif()
endif()
