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

#ifndef AtrialFibresClipperView_h
#define AtrialFibresClipperView_h

#include <berryISelectionListener.h>
#include <QmitkAbstractView.h>
#include <QMessageBox>
#include <vtkIdList.h>
#include <vtkActor.h>
#include <CemrgAtriaClipper.h>
#include "ui_AtrialFibresClipperViewControls.h"
#include "ui_AtrialFibresClipperViewLabels.h"
#include "ui_AtrialFibresViewUIMeshing.h"

/**
  \brief AtrialFibresClipperView

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
class AtrialFibresClipperView : public QmitkAbstractView {

    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)
    Q_OBJECT

public:

    static const std::string VIEW_ID;
    static void SetDirectoryFile(const QString directory, const QString fileName, const bool isAuto);
    ~AtrialFibresClipperView();

protected slots:

    /// \brief Called when the user clicks the GUI button
    void CtrLines();
    void CtrPlanes();
    void ClipperImage();
    void CtrPlanesPlacer();
    void CtrLinesSelector(int);

protected:

    virtual void CreateQtPartControl(QWidget *parent) override;
    virtual void SetFocus() override;
    /// \brief called by QmitkFunctionality when DataManager's selection has changed
    virtual void OnSelectionChanged(
            berry::IWorkbenchPart::Pointer source, const QList<mitk::DataNode::Pointer>& nodes) override;

    Ui::AtrialFibresClipperViewControls m_Controls;
    Ui::AtrialFibresClipperViewLabels m_Labels;
    Ui::AtrialFibresViewUIMeshing m_UIMeshing;

private:

    void iniPreSurf();
    void Visualiser(double opacity=1.0);
    void PickCallBack();
    void ManualCutterCallBack();
    static void KeyCallBackFunc(vtkObject*, long unsigned int, void* ClientData, void*);

    void InitialisePickerObjects();

    static QString fileName;
    static QString directory;
    static bool isAutomatic;

    mitk::Surface::Pointer surface;
    vtkSmartPointer<vtkActor> surfActor;
    std::vector<int> pickedSeedLabels;
    vtkSmartPointer<vtkIdList> pickedSeedIds;
    vtkSmartPointer<vtkPolyData> pickedLineSeeds;
    vtkSmartPointer<vtkPolyData> pickedCutterSeeds;
    std::unique_ptr<CemrgAtriaClipper> clipper;
    std::vector<vtkSmartPointer<vtkActor>> clipperActors;
    QDialog* inputs;
    double maxScalar, minScalar;
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkCallbackCommand> callBack;
    vtkSmartPointer<vtkRenderWindowInteractor> interactor;

};

#endif // AtrialFibresClipperView_h
