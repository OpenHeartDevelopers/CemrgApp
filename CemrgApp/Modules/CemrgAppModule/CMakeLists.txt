mitk_create_module(CemrgAppModule
  DEPENDS PUBLIC MitkSegmentation MitkQtWidgetsExt MitkQtWidgets MitkAlgorithmsExt MitkCore
  PACKAGE_DEPENDS PRIVATE VMTK ITK|ITKLabelMap+ITKMathematicalMorphology VTK|IOImage
)

# zlib is used by CemrgCommonUtils to inspect and repair NIfTI headers inside gzipped
# (.nii.gz) files, which cannot be patched in place.
find_package(ZLIB REQUIRED)
target_link_libraries(MitkCemrgAppModule PRIVATE ZLIB::ZLIB)

# VMTK is located by CMake/FindVMTK.cmake, which publishes its libraries through the
# directory-scoped link_libraries() command. That only affects targets created afterwards, so
# it does not reliably reach this module on a re-configure. Link the imported targets directly.
target_link_libraries(MitkCemrgAppModule PRIVATE
  vtkvmtkCommon
  vtkvmtkComputationalGeometry
  vtkvmtkContrib
  vtkvmtkDifferentialGeometry
  vtkvmtkIO
  vtkvmtkITK
  vtkvmtkMisc
  vtkvmtkSegmentation
)

add_subdirectory(cmdapps)

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
