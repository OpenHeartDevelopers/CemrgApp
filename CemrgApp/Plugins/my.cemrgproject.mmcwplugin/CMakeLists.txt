project(my_cemrgproject_mmcwplugin)

mitk_create_plugin(
  EXPORT_DIRECTIVE AWESOME_EXPORTS
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt MyCemrgLib
  PACKAGE_DEPENDS ITK
)
