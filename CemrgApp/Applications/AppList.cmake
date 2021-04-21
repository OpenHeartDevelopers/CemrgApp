option(MITK_BUILD_APP_CemrgApp "Build the MITK - CemrgApp" ON)
# option(MITK_BUILD_APP_Workbench "Switch the MITK Workbench executable off." OFF)

set(MITK_APPS
  MainApp^^MITK_BUILD_APP_CemrgApp^^CemrgApp
)
