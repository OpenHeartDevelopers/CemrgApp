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
 * Morphological Quantification
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef WallThicknessCalculationsView_h
#define WallThicknessCalculationsView_h

#include <QmitkAbstractView.h>
#include <mitkSurface.h>
#include "ui_WallThicknessCalculationsViewControls.h"
#include "ui_WallThicknessCalculationsViewUIMeshing.h"
#include "ui_WallThicknessCalculationsViewUIThickness.h"
#include "QmitkRenderWindow.h"
#include "mitkCommon.h"
#include "mitkDataStorage.h"
#include "mitkDataNode.h"
#include "mitkSurface.h"
#include "vtkRenderer.h"
#include "vtkTextActor.h"

/**
  \brief WallThicknessCalculationsView

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
class WallThicknessCalculationsView: public QmitkAbstractView {

    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)
    Q_OBJECT

public:

    static const std::string VIEW_ID;
    WallThicknessCalculationsView();

protected slots:

    /// \brief Called when the user clicks the GUI button
    void LoadDICOM();
    void ProcessIMGS();
    void ConvertNII();
    void CropIMGS();
    void ResampIMGS();
    void ApplyFilter();
    void SegmentIMGS();
    void SelectROI();
    void SelectLandmarks();
    void CombineSegs();
    void ClipperPV();
    void MorphologyAnalysis();
    void ThicknessAnalysis();
    void ConvertNRRD();
    void Browse();
    void ThicknessCalculator();
    void Reset();

protected:

    virtual void CreateQtPartControl(QWidget *parent) override;
    virtual void SetFocus() override;
    /// \brief called by QmitkFunctionality when DataManager's selection has changed
    virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer source, const QList<mitk::DataNode::Pointer>& nodes) override;

    Ui::WallThicknessCalculationsViewControls m_Controls;
    Ui::WallThicknessCalculationsViewUIMeshing m_UIMeshing;
    Ui::WallThicknessCalculationsViewUIThickness m_Thickness;

private:

    QString fileName;
    QString directory;
    const int APPENDAGECUT = 19;
    const int APPENDAGEUNCUT = 20;
};

#endif // WallThicknessCalculationsView_h
