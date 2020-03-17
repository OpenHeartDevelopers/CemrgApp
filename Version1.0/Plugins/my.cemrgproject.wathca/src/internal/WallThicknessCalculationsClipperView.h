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

#ifndef WallThicknessCalculationsClipperView_h
#define WallThicknessCalculationsClipperView_h

#include <QmitkAbstractView.h>
#include <QMessageBox>
#include <vtkIdList.h>
#include <vtkActor.h>
#include <CemrgAtriaClipper.h>
#include "ui_WallThicknessCalculationsClipperViewControls.h"
#include "ui_WallThicknessCalculationsClipperViewLabels.h"
#include "ui_WallThicknessCalculationsViewUIMeshing.h"


class WallThicknessCalculationsClipperView : public QmitkAbstractView {

    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)
    Q_OBJECT

public:

    static const std::string VIEW_ID;
    virtual void CreateQtPartControl(QWidget *parent);
    static void SetDirectoryFile(const QString directory, const QString fileName);
    ~WallThicknessCalculationsClipperView();

protected slots:

    /// \brief Called when the user clicks the GUI button
    void CtrLines();
    void CtrPlanes();
    void ClipperImage();
    void CtrPlanesPlacer();
    void CtrLinesSelector(int);

protected:

    virtual void SetFocus();
    Ui::WallThicknessCalculationsClipperViewControls m_Controls;
    Ui::WallThicknessCalculationsClipperViewLabels m_Labels;
    Ui::WallThicknessCalculationsViewUIMeshing m_UIMeshing;

private:

    void iniPreSurf();
    void Visualiser();
    void PickCallBack();
    void ManualCutterCallBack();
    static void KeyCallBackFunc(vtkObject*, long unsigned int, void* ClientData, void*);

    mitk::Surface::Pointer surface;
    vtkSmartPointer<vtkActor> surfActor;
    std::vector<int> pickedSeedLabels;
    vtkSmartPointer<vtkIdList> pickedSeedIds;
    vtkSmartPointer<vtkPolyData> pickedLineSeeds;
    vtkSmartPointer<vtkPolyData> pickedCutterSeeds;
    std::unique_ptr<CemrgAtriaClipper> clipper;
    std::vector<vtkSmartPointer<vtkActor>> clipperActors;

    QDialog* inputs;
    static QString fileName;
    static QString directory;
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkCallbackCommand> callBack;
    vtkSmartPointer<vtkRenderWindowInteractor> interactor;
};

#endif // WallThicknessCalculationsClipperView_h
