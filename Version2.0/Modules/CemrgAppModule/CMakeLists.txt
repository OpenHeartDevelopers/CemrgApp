mitk_create_module(CemrgAppModule
  DEPENDS PUBLIC MitkQtWidgetsExt MitkQtWidgets MitkAlgorithmsExt MitkCore
  PACKAGE_DEPENDS PRIVATE VMTK ITK
)

add_subdirectory(cmdapps)
