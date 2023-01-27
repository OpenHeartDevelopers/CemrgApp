project(kcl_cemrgapp_atrialfibres)

mitk_create_plugin(
  EXPORT_DIRECTIVE ATRIALFIBRES_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkSegmentation MitkQtWidgetsExt MitkCemrgAppModule
  PACKAGE_DEPENDS ITK
)
