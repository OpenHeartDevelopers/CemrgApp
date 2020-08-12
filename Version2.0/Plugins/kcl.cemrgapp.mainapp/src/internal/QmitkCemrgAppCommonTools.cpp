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
 * CemrgApp Common Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbench.h>
#include <berryIWorkbenchWindow.h>
#include <berryIWorkbenchPage.h>

// Qmitk
#include <mitkCoreObjectFactory.h>
#include <mitkProgressBar.h>
#include <QmitkIOUtil.h>

// CemrgAppModule
#include <CemrgCommonUtils.h>
#include "QmitkCemrgAppCommonTools.h"

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QSignalMapper>

const std::string QmitkCemrgAppCommonTools::VIEW_ID = "org.mitk.views.cemrgappcommontools";

void QmitkCemrgAppCommonTools::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::LoadMesh);
    connect(m_Controls.button_2, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::ConvertToCarto);
}

void QmitkCemrgAppCommonTools::SetFocus() {
}

void QmitkCemrgAppCommonTools::OnSelectionChanged(
        berry::IWorkbenchPart::Pointer /*source*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

void QmitkCemrgAppCommonTools::LoadMesh() {

    QString path = "";
    path = QFileDialog::getOpenFileName(
                NULL, "Open Mesh Data File", QmitkIOUtil::GetFileOpenFilterString());
    CemrgCommonUtils::AddToStorage(
                CemrgCommonUtils::LoadVTKMesh(path.toStdString()), "Mesh", this->GetDataStorage());
}

void QmitkCemrgAppCommonTools::ConvertToCarto() {

    QString path = "";
    path = QFileDialog::getOpenFileName(
                NULL, "Open Mesh Data File", QmitkIOUtil::GetFileOpenFilterString());

    if (path.isEmpty() || !path.endsWith(".vtk")) {
        QMessageBox::warning(NULL, "Attention", "Select Correct Input File!");
        return;
    }

    CemrgCommonUtils::ConvertToCarto(path.toStdString());
    QMessageBox::information(NULL, "Attention", "Conversion Completed!");
}
