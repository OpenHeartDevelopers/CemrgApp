SET(SRC_CPP_FILES
)

SET(INTERNAL_CPP_FILES
  my_cemrgproject_renderwindoweditor_Activator.cpp
  QmitkCemrgRenderWindowEditor.cpp
)

SET(MOC_H_FILES
  src/internal/my_cemrgproject_renderwindoweditor_Activator.h
  src/internal/QmitkCemrgRenderWindowEditor.h
)

SET(UI_FILES
)

SET(CACHED_RESOURCE_FILES
  plugin.xml
)

SET(QRC_FILES
)

SET(CPP_FILES )

foreach(file ${SRC_CPP_FILES})
  SET(CPP_FILES ${CPP_FILES} src/${file})
endforeach(file ${SRC_CPP_FILES})

foreach(file ${INTERNAL_CPP_FILES})
  SET(CPP_FILES ${CPP_FILES} src/internal/${file})
endforeach(file ${INTERNAL_CPP_FILES})
