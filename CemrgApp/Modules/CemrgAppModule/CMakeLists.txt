mitk_create_module(CemrgAppModule
  DEPENDS PUBLIC MitkSegmentation MitkQtWidgetsExt MitkQtWidgets MitkAlgorithmsExt MitkCore
  PACKAGE_DEPENDS PRIVATE VMTK ITK
)

add_subdirectory(cmdapps)
add_subdirectory(test)
