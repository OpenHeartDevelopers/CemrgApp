mitk_create_module(CemrgAppModule
  DEPENDS PUBLIC MitkSegmentation MitkQtWidgetsExt MitkQtWidgets MitkAlgorithmsExt MitkCore
  PACKAGE_DEPENDS PRIVATE VMTK ITK|ITKLabelMap+ITKMathematicalMorphology VTK|IOImage
)

# NIfTI headers inside .nii.gz cannot be patched in place
find_package(ZLIB REQUIRED)
target_link_libraries(MitkCemrgAppModule PRIVATE ZLIB::ZLIB)

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
