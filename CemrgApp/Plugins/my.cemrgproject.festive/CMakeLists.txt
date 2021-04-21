project(my_cemrgproject_festive)

mitk_create_plugin(
  EXPORT_DIRECTIVE FESTIVE_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt MyCemrgLib
)
