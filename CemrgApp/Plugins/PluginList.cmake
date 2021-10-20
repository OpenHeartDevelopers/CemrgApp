list(APPEND MITK_PLUGIN_REGEX_LIST "^kcl_cemrgapp_[a-zA-Z0-9_]+$")

set(MITK_PLUGINS
  # plugins from this module
  kcl.cemrgapp.mainapp:ON
  kcl.cemrgapp.easi:ON
  kcl.cemrgapp.powertrans:ON
  kcl.cemrgapp.atrialfibres:ON
  kcl.cemrgapp.mmcwplugin:ON
  kcl.cemrgapp.mmeasurement:ON
  kcl.cemrgapp.scar:ON
  kcl.cemrgapp.wathca:ON
)
