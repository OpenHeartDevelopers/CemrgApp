project(kcl_cemrgapp_mmcwplugin)

mitk_create_plugin(
  EXPORT_DIRECTIVE MMCWPLUGIN_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt MitkCemrgAppModule
  PACKAGE_DEPENDS PRIVATE ITK
)
