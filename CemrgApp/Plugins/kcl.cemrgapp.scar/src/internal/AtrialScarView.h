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
 * Scar Quantification
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
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
#include "ui_AtrialScarViewUIcemrgnet.h"

/**
  \brief AtrialScarView

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
class AtrialScarView: public QmitkAbstractView {

    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)
    Q_OBJECT

public:

    static const std::string VIEW_ID;
    AtrialScarView();

protected slots:

    /// \brief Called when the user clicks the GUI button
    void LoadDICOM();
    void ProcessIMGS();
    void ConvertNII();
    void AnalysisChoice();
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
    void ScarDebug();
    void Threshold();
    void Sphericity();
    void ExtraCalcs();
    void ResetMain();

protected:

    virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer source, const QList<mitk::DataNode::Pointer>& nodes) override;
    virtual void CreateQtPartControl(QWidget *parent) override;
    virtual void SetFocus() override;
    virtual void NodeAdded(const mitk::DataNode* node);

    Ui::AtrialScarViewUIScar m_UIScar;
    Ui::AtrialScarViewUISQuant m_UISQuant;
    Ui::AtrialScarViewControls m_Controls;
    Ui::AtrialScarViewUIMeshing m_UIMeshing;
    Ui::AtrialScarViewUIcemrgnet m_UIcemrgnet;

private:

    void AutomaticAnalysis();
    void Reset(bool allItems);

    // helper functions
    bool RequestProjectDirectoryFromUser();

    QString fileName;
    QString directory;
    QString debugSCARname;
    QString alternativeNiftiFolder;
    std::unique_ptr<CemrgScar3D> scar;
    bool _useDockerInPlugin = false; // change to FALSE to use MIRTK static libraries
};

#endif // AtrialScarView_h
