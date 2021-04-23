project(kcl_cemrgapp_mmeasurement)

mitk_create_plugin(
  EXPORT_DIRECTIVE MMEASUREMENT_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt MitkCemrgAppModule
)
