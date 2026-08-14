mitk_create_module(CemrgAppModule
  DEPENDS PUBLIC MitkSegmentation MitkQtWidgetsExt MitkQtWidgets MitkAlgorithmsExt MitkCore
  PACKAGE_DEPENDS PRIVATE VMTK ITK|ITKLabelMap+ITKMathematicalMorphology VTK|IOImage
)

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
