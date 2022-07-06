#cmake_minimum_required(VERSION 3.13)
if(MITK_BUILD_APP_CemrgApp)
  project(MainApp)

  set(MY_APP_NAME CemrgApp)

  # Create a cache entry for the provisioning file which is used to export
  # the file name in the MITKConfig.cmake file. This will keep external projects
  # which rely on this file happy.
  set(MITK_CEMRGAPP_PROVISIONING_FILE
    "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/CemrgApp.provisioning"
    CACHE INTERNAL "CemrgApp provisioning file" FORCE)

  # The _plugins variable can be used to specifically load certain plugins. This
  # might be useful to load them without having to build all the plugins, however
  # if all the MITK plugins are being built, these are not needed.
  set(_plugins
    # Require external plug-ins
    org.blueberry.ui.qt

    org.mitk.gui.qt.datamanager
    org.mitk.gui.qt.basicimageprocessing
    org.mitk.gui.qt.dicombrowser
    #org.mitk.gui.qt.diffusionimaging
    #org.mitk.gui.qt.diffusionimaging.registration
    org.mitk.gui.qt.pointsetinteraction
    org.mitk.gui.qt.segmentation
    org.mitk.gui.qt.remeshing

    # Enable plug-ins from this project
    kcl.cemrgapp.mainapp
    kcl.cemrgapp.easi
    kcl.cemrgapp.powertrans
    kcl.cemrgapp.atrialfibres
    kcl.cemrgapp.mmcwplugin
    kcl.cemrgapp.mmeasurement
    kcl.cemrgapp.scar
    kcl.cemrgapp.wathca
  )

  # Plug-ins listed below will not be
  #  - added as a build-time dependency to the executable
  #  - listed in the provisioning file for the executable
  #  - installed if they are external plug-ins

  set(_exclude_plugins
    org.blueberry.test
    org.blueberry.uitest
    org.mitk.gui.qt.coreapplication
    org.mitk.gui.qt.diffusionimagingapp
    org.mitk.example.gui.customviewer
    org.mitk.example.gui.customviewer.views
    org.mitk.example.gui.selectionservicemitk
    org.mitk.example.gui.selectionservicemitk.views
    org.mitk.example.gui.selectionserviceqt
    org.mitk.example.gui.extensionpointcontribution
    org.mitk.example.gui.extensionpointdefinition
    org.mitk.example.gui.minimalapplication
    org.mitk.example.gui.multipleperspectives
  )
  set(_src_files
    CemrgApp.cpp
  )

  qt5_add_resources(_src_files splashscreen.qrc)

  mitkFunctionCreateBlueBerryApplication(
    NAME ${MY_APP_NAME}
    DESCRIPTION "MITK - ${MY_APP_NAME} Application"
    #PLUGINS ${_plugins}
    EXCLUDE_PLUGINS ${_exclude_plugins}
    SOURCES ${_src_files}
  )

  mitk_use_modules(TARGET ${MY_APP_NAME}
    MODULES MitkAppUtil
    PACKAGES Qt5|Widgets
  )

  # Add meta dependencies (e.g. on auto-load modules from depending modules)
  if(TARGET ${CMAKE_PROJECT_NAME}-autoload)
    add_dependencies(${MY_APP_NAME} ${CMAKE_PROJECT_NAME}-autoload)
  endif()

  #Setting application icon for macOS systems
  set_target_properties(${MY_APP_NAME} PROPERTIES MACOSX_BUNDLE_ICON_FILE "icon.icns")

  if(APPLE)
    install(FILES "icons/icon.icns" DESTINATION "${MY_APP_NAME}.app/Contents/Resources")
    set(MACOSX_BUNDLE_NAMES ${MACOSX_BUNDLE_NAMES} ${MY_APP_NAME})
  endif(APPLE)

  # Add a build time dependency to legacy BlueBerry bundles.
  if(MITK_MODULES_ENABLED_PLUGINS)
    add_dependencies(${MY_APP_NAME} ${MITK_MODULES_ENABLED_PLUGINS})
  endif()

endif()
