PROJECT(my_cemrgproject_cemrgapp)

mitk_create_plugin(
  EXPORT_DIRECTIVE CEMRGAPP_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt
)
