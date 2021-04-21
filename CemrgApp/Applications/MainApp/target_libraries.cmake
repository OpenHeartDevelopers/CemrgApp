# A list of plug-in targets which should be automatically enabled
# (or be available in external projects) for this application.

set(target_libraries
  # Require external plug-ins
  org_blueberry_ui_qt

  org_mitk_gui_qt_datamanager
  org_mitk_gui_qt_basicimageprocessing
  org_mitk_gui_qt_dicom
  #org_mitk_gui_qt_diffusionimaging
  #org_mitk_gui_qt_diffusionimaging_registration
  org_mitk_gui_qt_pointsetinteraction
  org_mitk_gui_qt_segmentation

  # Enable plug-ins from this project
  kcl_cemrgapp_mainapp
)
