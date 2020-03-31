set(SRC_CPP_FILES

)

set(INTERNAL_CPP_FILES
  kcl_cemrgapp_scar_Activator.cpp
  AtrialScarClipperView.cpp
  AtrialScarView.cpp
  YZSegView.cpp
  ScarCalculationsView.cpp
)

set(UI_FILES
  src/internal/AtrialScarClipperViewControls.ui
  src/internal/AtrialScarClipperViewLabels.ui
  src/internal/AtrialScarViewControls.ui
  src/internal/AtrialScarViewUIMeshing.ui
  src/internal/AtrialScarViewUIScar.ui
  src/internal/AtrialScarViewUISQuant.ui
  src/internal/AtrialScarViewUIcemrgnet.ui
  src/internal/YZSegViewControls.ui
  src/internal/ScarCalculationsViewControls.ui
  src/internal/ScarCalculationsViewUICorridor.ui
)

set(MOC_H_FILES
  src/internal/kcl_cemrgapp_scar_Activator.h
  src/internal/AtrialScarClipperView.h
  src/internal/AtrialScarView.h
  src/internal/YZSegView.h
  src/internal/ScarCalculationsView.h
)

# list of resource files which can be used by the plug-in
# system without loading the plug-ins shared library,
# for example the icon used in the menu and tabs for the
# plug-in views in the workbench
set(CACHED_RESOURCE_FILES
  resources/icon.xpm
  plugin.xml
)

# list of Qt .qrc files which contain additional resources
# specific to this plugin
set(QRC_FILES

)

set(CPP_FILES )

foreach(file ${SRC_CPP_FILES})
  set(CPP_FILES ${CPP_FILES} src/${file})
endforeach(file ${SRC_CPP_FILES})

foreach(file ${INTERNAL_CPP_FILES})
  set(CPP_FILES ${CPP_FILES} src/internal/${file})
endforeach(file ${INTERNAL_CPP_FILES})
