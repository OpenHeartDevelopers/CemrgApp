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
 * Power Transmitter Calculations Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * angela.lee@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef powertransView_h
#define powertransView_h

#include <berryISelectionListener.h>
#include <QmitkAbstractView.h>
#include <mitkSurface.h>
#include "ui_powertransViewControls.h"
#include "ui_powertransViewUIRibSpacing.h"
#include "ui_powertransViewUIAhaInput.h"
#include "CemrgPower.h"

/**
  \brief powertransView

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
class powertransView: public QmitkAbstractView {

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
    void VolumeRendering();
    void MapPowerTop();
    void ResetRibSpacing();
    void LandmarkSelection();
    void MapPowerTransLM();
    //void ConfirmSite();
    void CalcPowerTop();
    void LoadMesh();
    void CalculatePower();
    void MapAHATop();
    void AHALandmarkSelection();
    void MapAHAfromInput();
    void MapAHA();
    void Reset();

protected:

    virtual void CreateQtPartControl(QWidget *parent) override;
    virtual void SetFocus() override;
    /// \brief called by QmitkFunctionality when DataManager's selection has changed
    virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer source, const QList<mitk::DataNode::Pointer>& nodes) override;

    Ui::powertransViewControls m_Controls;
    Ui::powertransViewUIRibSpacing m_UIRibSpacing;
    Ui::powertransViewUIAhaInput m_UIAhaInput;

private:

    QString directory;
    int ribSpacing = 0;
    std::unique_ptr<CemrgPower> power;
};

#endif // powertransView_h
