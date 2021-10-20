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
 * CemrgApp Main App
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#include "QmitkCemrgApplication.h"
#include <berryPlatformUI.h>
#include "QmitkCemrgWorkbenchAdvisor.h"
#include <QDebug>

QmitkCemrgApplication::QmitkCemrgApplication() {
}

QmitkCemrgApplication::QmitkCemrgApplication(const QmitkCemrgApplication& other) : QObject(other.parent()) {
    Q_UNUSED(other)
    throw std::runtime_error("Copy constructor not implemented");
}

QVariant QmitkCemrgApplication::Start(berry::IApplicationContext* /*context*/) {

    berry::Display* display = berry::PlatformUI::CreateDisplay();
    int code = berry::PlatformUI::CreateAndRunWorkbench(display, new QmitkCemrgWorkbenchAdvisor());

    //exit the application with an appropriate return code
    return code == berry::PlatformUI::RETURN_RESTART ? EXIT_RESTART : EXIT_OK;
}

void QmitkCemrgApplication::Stop() {
}
