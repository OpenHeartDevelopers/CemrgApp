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

#ifndef QMITKCEMRGWORKBENCHADVISOR_H_
#define QMITKCEMRGWORKBENCHADVISOR_H_
#ifdef __MINGW32__
// We need to inlclude winbase.h here in order to declare
// atomic intrinsics like InterlockedIncrement correctly.
// Otherwhise, they would be declared wrong within qatomic_windows.h .
#include <windows.h>
#endif

#include <berryQtWorkbenchAdvisor.h>

class QmitkCemrgWorkbenchAdvisor: public berry::QtWorkbenchAdvisor {

public:

    static const QString DEFAULT_PERSPECTIVE_ID;

    void Initialize(berry::IWorkbenchConfigurer::Pointer configurer);
    berry::WorkbenchWindowAdvisor* CreateWorkbenchWindowAdvisor(berry::IWorkbenchWindowConfigurer::Pointer configurer) override;
    QString GetInitialWindowPerspectiveId() override;
};

#endif /*QMITKCEMRGWORKBENCHADVISOR_H_*/
