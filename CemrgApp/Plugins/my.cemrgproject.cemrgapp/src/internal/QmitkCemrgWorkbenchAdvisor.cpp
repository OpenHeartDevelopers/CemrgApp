/*=========================================================================

 Program:   Medical Imaging & Interaction Toolkit
 Language:  C++
 Date:      $Date$
 Version:   $Revision$

 Copyright (c) German Cancer Research Center, Division of Medical and
 Biological Informatics. All rights reserved.
 See MITKCopyright.txt or http://www.mitk.org/copyright.html for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/
/*=========================================================================
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#include "QmitkCemrgWorkbenchAdvisor.h"
#include "internal/mitkCemrgAppPluginActivator.h"
#include <QmitkExtWorkbenchWindowAdvisor.h>
#include <mitkWorkbenchUtil.h>


const QString QmitkCemrgWorkbenchAdvisor::DEFAULT_PERSPECTIVE_ID = "my.cemrgproject.CemrgJBPerspective";

void QmitkCemrgWorkbenchAdvisor::Initialize(berry::IWorkbenchConfigurer::Pointer configurer) {

    berry::QtWorkbenchAdvisor::Initialize(configurer);

    configurer->SetSaveAndRestore(true);

    ctkPluginContext* context = mitk::CemrgAppPluginActivator::GetDefault()->GetPluginContext();
    mitk::WorkbenchUtil::SetDepartmentLogoPreference(":/CemrgApp/MyLogo.png", context);
}

berry::WorkbenchWindowAdvisor*
QmitkCemrgWorkbenchAdvisor::CreateWorkbenchWindowAdvisor(berry::IWorkbenchWindowConfigurer::Pointer configurer) {

    // -------------------------------------------------------------------
    // Here you could pass your custom Workbench window advisor
    // -------------------------------------------------------------------
    QmitkExtWorkbenchWindowAdvisor* advisor = new QmitkExtWorkbenchWindowAdvisor(this, configurer);

    advisor->ShowViewToolbar(false);
    advisor->SetWindowIcon(":/CemrgApp/icon_research.xpm");

    return advisor;
}

QString QmitkCemrgWorkbenchAdvisor::GetInitialWindowPerspectiveId() {

    return DEFAULT_PERSPECTIVE_ID;
}
