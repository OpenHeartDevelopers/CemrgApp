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
 * Wall Thickness Calculations (WATHCA) Plugin for MITK
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
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


class WallThicknessCalculationsView : public QmitkAbstractView {

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

    virtual void SetFocus();

    Ui::WallThicknessCalculationsViewControls m_Controls;
    Ui::WallThicknessCalculationsViewUIMeshing m_UIMeshing;
    Ui::WallThicknessCalculationsViewUIThickness m_Thickness;

private:

    mitk::Surface::Pointer ReadVTKMesh(std::string meshPath);
    void AddToStorage(std::string nodeName, mitk::BaseData* data);

    QString fileName;
    QString directory;

    const int APPENDAGECUT   = 19;
    const int APPENDAGEUNCUT = 20;
};

#endif // WallThicknessCalculationsView_h
