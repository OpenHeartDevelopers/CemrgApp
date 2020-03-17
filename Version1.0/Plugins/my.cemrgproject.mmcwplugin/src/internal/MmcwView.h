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
 * Motion Measurement of Cardiac Wall (MMCW) Plugin for MITK
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef MmcwView_h
#define MmcwView_h

#include <berryISelectionListener.h>
#include <QmitkAbstractView.h>
#include <mitkSurface.h>
#include "ui_MmcwViewControls.h"
#include "ui_MmcwViewUIMeshing.h"
#include "ui_MmcwViewUITracking.h"
#include "ui_MmcwViewUIApplying.h"

/*!
  \brief MmcwView

  \warning  This application module is not yet documented. Use "svn blame/praise/annotate" and ask the author to provide basic documentation.

  \sa QmitkFunctionality
  \ingroup Functionalities
*/
class MmcwView : public QmitkAbstractView
{  
    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)
    Q_OBJECT

public:

    static const std::string VIEW_ID;
    virtual void CreateQtPartControl(QWidget *parent);

protected slots:

    /// \brief Called when the user clicks the GUI button
    void LoadDICOM();
    void ProcessIMGS();
    void ConvertNII();
    void CropinIMGS();
    void ResampIMGS();
    void SegmentIMGS();
    void CreateSurf();
    void TrackingIMGS();
    void BrowseT(const QString& buttDir);
    void BrowseA(const QString& buttDir);
    void Tracking();
    void Applying();
    void Demoings();
    void LandmarkSelection();
    void Plotting();
    void Reset();

protected:

    virtual void SetFocus();
    /// \brief called by QmitkAbstractView when DataManager's selection has changed
    virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer source, const QList<mitk::DataNode::Pointer>& nodes);

    Ui::MmcwViewControls m_Controls;
    Ui::MmcwViewUIMeshing m_UIMeshing;
    Ui::MmcwViewUITracking m_UITracking;
    Ui::MmcwViewUIApplying m_UIApplying;

private:

    mitk::Surface::Pointer ReadVTKMesh(std::string meshPath);
    void AddToStorage(std::string nodeName, mitk::BaseData *data);
    QString directory;
    int timePoints;
};

#endif // MmcwView_h
