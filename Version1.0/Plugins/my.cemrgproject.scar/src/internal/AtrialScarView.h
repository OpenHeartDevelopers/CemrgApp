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
 * Atrial Scar (AS) Plugin for MITK
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef AtrialScarView_h
#define AtrialScarView_h

#include <berryISelectionListener.h>
#include <QmitkAbstractView.h>
#include <mitkSurface.h>
#include <CemrgScar3D.h>
#include "ui_AtrialScarViewUIScar.h"
#include "ui_AtrialScarViewUISQuant.h"
#include "ui_AtrialScarViewControls.h"
#include "ui_AtrialScarViewUIMeshing.h"


class AtrialScarView : public QmitkAbstractView {

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
    void SegmentIMGS();
    void Registration();
    void Register();
    void Transform();
    void ClipPVeins();
    void CreateSurf();
    void ClipperMV();
    void SelectLandmarks();
    void ClipMitralValve();
    void ScarMap();
    void Threshold();
    void Sphericity();
    void ResetMain();

protected:

    /// \brief called by QmitkFunctionality when DataManager's selection has changed
    virtual void SetFocus();
    virtual void NodeAdded(const mitk::DataNode* node);

    Ui::AtrialScarViewUIScar m_UIScar;
    Ui::AtrialScarViewUISQuant m_UISQuant;
    Ui::AtrialScarViewControls m_Controls;
    Ui::AtrialScarViewUIMeshing m_UIMeshing;

private:

    void Reset(bool allItems);
    mitk::DataNode::Pointer AddToStorage(std::string nodeName, mitk::BaseData* data);
    mitk::Surface::Pointer ReadVTKMesh(std::string meshPath);

    QString fileName;
    QString directory;
    std::unique_ptr<CemrgScar3D> scar;
};

#endif // AtrialScarView_h
