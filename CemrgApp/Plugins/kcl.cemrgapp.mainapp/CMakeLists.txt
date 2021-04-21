project(kcl_cemrgapp_mainapp)

mitk_create_plugin(
  EXPORT_DIRECTIVE MAINAPP_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt MitkCemrgAppModule
)
