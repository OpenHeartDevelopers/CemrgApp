project(my_cemrgproject_wathca)

mitk_create_plugin(
  EXPORT_DIRECTIVE WATHCA_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt MyCemrgLib
  PACKAGE_DEPENDS ITK
)
