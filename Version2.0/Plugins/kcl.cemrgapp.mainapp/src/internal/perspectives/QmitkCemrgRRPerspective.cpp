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

#include "QmitkCemrgRRPerspective.h"
#include "berryIViewLayout.h"

QmitkCemrgRRPerspective::QmitkCemrgRRPerspective() {
}

QmitkCemrgRRPerspective::QmitkCemrgRRPerspective(const QmitkCemrgRRPerspective& other) : QObject() {

    Q_UNUSED(other)
    throw std::runtime_error("Copy constructor not implemented");
}

void QmitkCemrgRRPerspective::CreateInitialLayout(berry::IPageLayout::Pointer layout) {

    QString editorArea = layout->GetEditorArea();
    layout->AddView("org.mitk.views.motionmeasurement", berry::IPageLayout::LEFT, 0.35f, editorArea);
    berry::IFolderLayout::Pointer folder = layout->CreateFolder(
                "folder", berry::IPageLayout::BOTTOM, 0.8f, "org.mitk.views.motionmeasurement");
    folder->AddView("org.mitk.views.datamanager");
    berry::IViewLayout::Pointer lo = layout->GetViewLayout("org.mitk.views.motionmeasurement");
    lo->SetCloseable(false);
}
