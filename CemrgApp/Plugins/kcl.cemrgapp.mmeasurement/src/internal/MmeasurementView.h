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
 * Anatomical Measurements
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef MmeasurementView_h
#define MmeasurementView_h

#include <berryISelectionListener.h>
#include <QmitkAbstractView.h>
#include <QmitkPlotWidget.h>
#include "ui_MmeasurementViewControls.h"
#include "ui_MmeasurementViewUIApplying.h"
#include "ui_MmeasurementViewUITracking.h"

/**
  \brief MmeasurementView

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
class MmeasurementView : public QmitkAbstractView {

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
    void TrackingButton();
    void BrowseT(const QString& buttDir);
    void BrowseA(const QString& buttDir);
    void Tracking();
    void Applying();
    void SelectLandmark();
    void WriteFileButton();
    void CalcDistButton();
    void CalcPeriButton();
    void CalcAreaButton();
    void FindCentre();
    void Reset();

protected:

    virtual void CreateQtPartControl(QWidget *parent) override;
    virtual void SetFocus() override;
    /// \brief called by QmitkFunctionality when DataManager's selection has changed
    virtual void OnSelectionChanged(
            berry::IWorkbenchPart::Pointer source, const QList<mitk::DataNode::Pointer>& nodes) override;

    Ui::MmeasurementViewControls m_Controls;
    Ui::MmeasurementViewUITracking m_UITracking;
    Ui::MmeasurementViewUIApplying m_UIApplying;

private:

    int timePoints;
    int smoothness;
    QString directory;
    std::vector<double> plotValueVectors;
};

#endif // MmeasurementView_h
