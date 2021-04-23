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
 * Eikonal Activation Simulation (EASI) Plugin for MITK
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/


// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>

// Qmitk
#include "QmitkCemrgRenderWindowEditor.h"

// Qt
#include <QMessageBox>

// mitk image
#include <mitkImage.h>

//const std::string QmitkCemrgRenderWindowEditor::VIEW_ID = "org.mitk.views.renderwindoweditor";
const std::string QmitkCemrgRenderWindowEditor::EDITOR_ID = "org.mitk.editors.stdmultiwidget";

QmitkCemrgRenderWindowEditor::QmitkCemrgRenderWindowEditor()
    : m_standrdWindow(0)
{
}

QmitkRenderWindow *QmitkCemrgRenderWindowEditor::GetActiveQmitkRenderWindow() const
{
    if (m_standrdWindow) return m_standrdWindow->GetRenderWindow1();
    return 0;
}

//*** NEEDS DOUBLE CHECKING ***//
//*****************************//
QHash<QString, QmitkRenderWindow*> wnds;
QHash<QString, QmitkRenderWindow *> QmitkCemrgRenderWindowEditor::GetQmitkRenderWindows() const
{
    wnds.insert("axial", m_standrdWindow->GetRenderWindow1());
    wnds.insert("sagittal", m_standrdWindow->GetRenderWindow2());
    wnds.insert("coronal", m_standrdWindow->GetRenderWindow3());
    wnds.insert("3d", m_standrdWindow->GetRenderWindow4());
    return wnds;
}

QmitkRenderWindow *QmitkCemrgRenderWindowEditor::GetQmitkRenderWindow(const QString &id) const
{
    if (wnds.contains(id))
      return wnds[id];
    return 0;
}

mitk::Point3D QmitkCemrgRenderWindowEditor::GetSelectedPosition(const QString & /*id*/) const
{
    return m_standrdWindow->GetCrossPosition();
}

void QmitkCemrgRenderWindowEditor::SetSelectedPosition(const mitk::Point3D &pos, const QString &/*id*/)
{
    m_standrdWindow->MoveCrossToPosition(pos);
}

void QmitkCemrgRenderWindowEditor::EnableDecorations(bool enable, const QStringList &decorations)
{
    if (decorations.isEmpty() || decorations.contains(DECORATION_BORDER))
    {
      enable ? m_standrdWindow->EnableColoredRectangles()
             : m_standrdWindow->DisableColoredRectangles();
    }
    if (decorations.isEmpty() || decorations.contains(DECORATION_LOGO))
    {
      enable ? m_standrdWindow->EnableDepartmentLogo()
             : m_standrdWindow->DisableDepartmentLogo();
    }
    if (decorations.isEmpty() || decorations.contains(DECORATION_MENU))
    {
      m_standrdWindow->ActivateMenuWidget(enable);
    }
    if (decorations.isEmpty() || decorations.contains(DECORATION_BACKGROUND))
    {
      enable ? m_standrdWindow->EnableGradientBackground()
             : m_standrdWindow->DisableGradientBackground();
    }
}

bool QmitkCemrgRenderWindowEditor::IsDecorationEnabled(const QString &decoration) const
{
    if (decoration == DECORATION_BORDER)
    {
      return m_standrdWindow->IsColoredRectanglesEnabled();
    }
    else if (decoration == DECORATION_LOGO)
    {
      return m_standrdWindow->IsColoredRectanglesEnabled();
    }
    else if (decoration == DECORATION_MENU)
    {
      return m_standrdWindow->IsMenuWidgetEnabled();
    }
    else if (decoration == DECORATION_BACKGROUND)
    {
      return m_standrdWindow->GetGradientBackgroundFlag();
    }
    return false;
}

QStringList QmitkCemrgRenderWindowEditor::GetDecorations() const
{
    QStringList decorations;
    decorations << DECORATION_BORDER << DECORATION_LOGO << DECORATION_MENU << DECORATION_BACKGROUND;
    return decorations;
}

void QmitkCemrgRenderWindowEditor::SetFocus()
{
    if (m_standrdWindow != 0)
        m_standrdWindow->setFocus();
}

void QmitkCemrgRenderWindowEditor::CreateQtPartControl(QWidget* parent)
{
    QGridLayout* layout = new QGridLayout(parent);
    layout->setContentsMargins(0,0,0,0);

    m_standrdWindow = new QmitkStdMultiWidget(parent);
    m_RenderWindow0 = new QmitkRenderWindow(parent);
    m_RenderWindow1 = new QmitkRenderWindow(parent);

    layout->addWidget(m_standrdWindow, 0, 0, -1, 1);
    layout->addWidget(m_RenderWindow0, 0, 1,  1, 1);
    layout->addWidget(m_RenderWindow1, 1, 1,  1, 1);

    mitk::DataStorage::Pointer ds = this->GetDataStorage();
    m_standrdWindow->SetDataStorage(ds);

    //Initialize the QmitkRenderWindow as transversal to all data objects in DataStorage
    //mitk::TimeSlicedGeometry::Pointer geo = ds->ComputeBoundingGeometry3D(ds->GetAll());
    //m_RenderWindow->
    //mitk::RenderingManager::GetInstance()->InitializeViews(geo);

    this->RequestUpdate();
}


// void QmitkCemrgRenderWindowEditor::SetFocus()
// {
//   m_Controls.buttonPerformImageProcessing->setFocus();
// }
//
// void QmitkCemrgRenderWindowEditor::CreateQtPartControl(QWidget *parent)
// {
//   // create GUI widgets from the Qt Designer's .ui file
//   m_Controls.setupUi(parent);
//   connect(m_Controls.buttonPerformImageProcessing, &QPushButton::clicked, this, &QmitkCemrgRenderWindowEditor::DoImageProcessing);
// }
//
// void QmitkCemrgRenderWindowEditor::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*source*/,
//                                                 const QList<mitk::DataNode::Pointer>DataNode::Pointer>& /*nodes*/)
// {
//   // iterate all selected objects, adjust warning visibility
//   foreach (mitk::DataNode::Pointer node, nodes)
//   {
//     if (node.IsNotNull() && dynamic_cast<mitk::Image *>(node->GetData()))
//     {
//       m_Controls.labelWarning->setVisible(false);
//       m_Controls.buttonPerformImageProcessing->setEnabled(true);
//       return;
//     }
//   }
//
//   m_Controls.labelWarning->setVisible(true);
//   m_Controls.buttonPerformImageProcessing->setEnabled(false);
// }
//
