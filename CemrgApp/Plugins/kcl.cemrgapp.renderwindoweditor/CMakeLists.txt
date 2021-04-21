project(kcl_cemrgapp_renderwindoweditor)

mitk_create_plugin(
  EXPORT_DIRECTIVE RENDERWINDOWEDITOR_EXPORT
  EXPORTED_INCLUDE_SUFFIXES src
  MODULE_DEPENDS MitkQtWidgetsExt
  SUBPROJECTS MITK-CoreUI
)
