set(SRC_CPP_FILES

)

include_directories(../kcl.cemrgapp.atrialfibres/src/internal)

set(INTERNAL_CPP_FILES
        kcl_cemrgapp_atrialstrainmotion_Activator.cpp
        AtrialStrainMotionView.cpp
)

set(UI_FILES
        src/internal/AtrialStrainMotionViewControls.ui
        ../kcl.cemrgapp.atrialfibres/src/internal/AtrialFibresViewUIAnalysisSelector.ui
        ../kcl.cemrgapp.atrialfibres/src/internal/AtrialFibresViewUIMeshing.ui
        ../kcl.cemrgapp.atrialfibres/src/internal/AtrialFibresClipperViewControls.ui
        ../kcl.cemrgapp.atrialfibres/src/internal/AtrialFibresClipperViewLabels.ui
        ../kcl.cemrgapp.atrialfibres/src/internal/AtrialFibresClipperViewUICorridor.ui
        ../kcl.cemrgapp.atrialfibres/src/internal/AtrialFibresClipperViewUIRadius.ui
        ../kcl.cemrgapp.atrialfibres/src/internal/AtrialFibresViewControls.ui
        ../kcl.cemrgapp.atrialfibres/src/internal/AtrialFibresViewUIEditLabels.ui
        ../kcl.cemrgapp.atrialfibres/src/internal/AtrialFibresViewUIUacSelector.ui
        ../kcl.cemrgapp.atrialfibres/src/internal/AtrialFibresViewUIRemesh.ui
        ../kcl.cemrgapp.atrialfibres/src/internal/AtrialFibresViewUIConvert.ui
        ../kcl.cemrgapp.atrialfibres/src/internal/AtrialFibresLandmarksViewRough.ui
        ../kcl.cemrgapp.atrialfibres/src/internal/AtrialFibresLandmarksViewRough.ui
        ../kcl.cemrgapp.atrialfibres/src/internal/AtrialFibresVisualiseViewControls.ui
        ../kcl.cemrgapp.scar/src/internal/AtrialScarViewUIcemrgnet.ui
)

set(MOC_H_FILES
        src/internal/kcl_cemrgapp_atrialstrainmotion_Activator.h
        src/internal/AtrialStrainMotionView.h
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
