project(kcl_cemrgapp_powertrans)

mitk_create_plugin(
  EXPORT_DIRECTIVE POWERTRANS_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt MitkCemrgAppModule
)
