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
 * Motion Quantification
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
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

/**
  \brief MmcwView
  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.
  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
class MmcwView: public QmitkAbstractView {

    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)
    Q_OBJECT

public:

    static const std::string VIEW_ID;

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

    virtual void CreateQtPartControl(QWidget *parent) override;
    virtual void SetFocus() override;
    /// \brief called by QmitkFunctionality when DataManager's selection has changed
    virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer source, const QList<mitk::DataNode::Pointer>& nodes) override;

    Ui::MmcwViewControls m_Controls;
    Ui::MmcwViewUIMeshing m_UIMeshing;
    Ui::MmcwViewUITracking m_UITracking;
    Ui::MmcwViewUIApplying m_UIApplying;

private:

    int timePoints;
    QString directory;
};

#endif // MmcwView_h
