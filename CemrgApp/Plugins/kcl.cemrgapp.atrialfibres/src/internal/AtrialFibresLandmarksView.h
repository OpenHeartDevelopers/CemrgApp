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

#ifndef AtrialFibresLandmarksView_h
#define AtrialFibresLandmarksView_h

#include <berryISelectionListener.h>
#include <QmitkAbstractView.h>
#include <QMessageBox>
#include <vtkIdList.h>
#include <vtkActor.h>
#include <CemrgAtriaClipper.h>
#include <CemrgScarAdvanced.h>

#include "ui_AtrialFibresLandmarksViewControls.h"
#include "ui_AtrialFibresLandmarksViewRough.h"
#include "ui_AtrialFibresLandmarksViewRefined.h"

/**
  \brief AtrialFibresLandmarksView

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
typedef std::pair<vtkIdType, double> SeedRadiusPairType;
class AtrialFibresLandmarksView : public QmitkAbstractView {

    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)
    Q_OBJECT

public:

    static const std::string VIEW_ID;
    static void SetDirectoryFile(const QString directory, const QString fileName, const QString whichAtrium);
    ~AtrialFibresLandmarksView();

    // helper functions
    std::string GetShortcuts();
    std::string GetRoughPointsGuide();
    std::string GetRefinedPointsGiude();
    std::string GetStructureIdFromLabel(bool refinedLandmarks, int label);

    int GetIndex(std::vector<int> v, int value);

protected slots:

    /// \brief Called when the user clicks the GUI button
    void Help(bool firstTime=false);
    void SaveRoughPoints();
    void SaveRefinedPoints();

protected:

    virtual void CreateQtPartControl(QWidget *parent) override;
    virtual void SetFocus() override;
    /// \brief called by QmitkFunctionality when DataManager's selection has changed
    virtual void OnSelectionChanged(
            berry::IWorkbenchPart::Pointer source, const QList<mitk::DataNode::Pointer>& nodes) override;

    Ui::AtrialFibresLandmarksViewControls m_Controls;
    Ui::AtrialFibresLandmarksViewRough m_Rough;
    Ui::AtrialFibresLandmarksViewRefined m_Refined;

private:

    void iniPreSurf();
    void Visualiser(double opacity=1.0);
    void VisualiserAuto(double opacity);
    void VisualiserManual(double opacity);

    void SphereSourceVisualiser(vtkSmartPointer<vtkPolyData> pointSources, QString colour="1.0,0.0,0.0", double scaleFactor=0.01);
    void PickCallBack(bool refinedLandmarks=false);
    static void KeyCallBackFunc(vtkObject*, long unsigned int, void* ClientData, void*);

    void InitialisePickerObjects();

    void UserSelectPvLabel(bool refinedLandmarks=false);
    void UserSelectPvRoughLabel();
    void UserSelectPvRefinedLabel();

    void RoughUiEnableButtons();
    void RefinedUiEnableButtons();

    static QString fileName;
    static QString directory;
    static QString whichAtrium;

    mitk::Surface::Pointer surface;
    vtkSmartPointer<vtkActor> surfActor;
    bool isLeftAtrium;

    std::vector<int> roughSeedLabels;
    vtkSmartPointer<vtkIdList> roughSeedIds;
    vtkSmartPointer<vtkPolyData> roughLineSeeds;

    std::vector<int> refinedSeedLabels;
    vtkSmartPointer<vtkIdList> refinedSeedIds;
    vtkSmartPointer<vtkPolyData> refinedLineSeeds;

    QDialog* inputsRough;
    QDialog* inputsRefined;
    double maxScalar, minScalar;

    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkCallbackCommand> callBack;
    vtkSmartPointer<vtkRenderWindowInteractor> interactor;

};

#endif // AtrialFibresLandmarksView_h
