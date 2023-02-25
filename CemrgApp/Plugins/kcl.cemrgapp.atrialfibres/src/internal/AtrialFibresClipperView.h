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
#include <CemrgScarAdvanced.h>

#include "ui_AtrialFibresClipperViewControls.h"
#include "ui_AtrialFibresClipperViewLabels.h"
#include "ui_AtrialFibresViewUIMeshing.h"
#include "ui_AtrialFibresClipperViewUIRadius.h"
#include "ui_AtrialFibresClipperViewUICorridor.h"

/**
  \brief AtrialFibresClipperView

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
typedef std::pair<vtkIdType, double> SeedRadiusPairType;
class AtrialFibresClipperView : public QmitkAbstractView {

    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)
    Q_OBJECT

public:

    static const std::string VIEW_ID;
    static void SetDirectoryFile(const QString directory, const QString fileName, const bool isAuto);
    ~AtrialFibresClipperView();

    // helper functions
    void SetManualModeButtons(bool b);
    void SetAutomaticModeButtons(bool b);

    // helper functions
    std::string GetShortcuts();
    std::string GetHelp();
    bool IsPointSelectionControlsAvailable();
    bool IsClipperManualControlsAvailable();
    void UserSelectPvLabel();
    void LoadPickedSeedsFromFile();
    void CreateSphereClipperAndRadiiVectors(bool showOnRenderer);
    void SaveSphereClippers();
    void PrintCorridorIds();
    void UpdateClipperSeedIds(int newPickedId, int currentId);
    int GetPickedId();
    int GetUserFixMeshingLabel();

    inline void SetDebug(bool b){debugging = b;};
    inline void SetDebugOn(){SetDebug(true);};
    inline void SetDebugOff(){SetDebug(false);};

protected slots:

    /// \brief Called when the user clicks the GUI button
    // Manual Pipeline
    void CtrLines();
    void CtrPlanes();
    void ClipperImage();

    void CtrPlanesPlacer();
    void CtrLinesSelector(int);

    // Automatic Pipeline
    void SaveLabels();
    void ShowPvClippers();
    void InterPvSpacing();

    void PvClipperRadius();
    void PvClipperSelector(int);

    void ClipPVs();
    void Help();


protected:

    virtual void CreateQtPartControl(QWidget *parent) override;
    virtual void SetFocus() override;
    /// \brief called by QmitkFunctionality when DataManager's selection has changed
    virtual void OnSelectionChanged(
            berry::IWorkbenchPart::Pointer source, const QList<mitk::DataNode::Pointer>& nodes) override;

    Ui::AtrialFibresClipperViewControls m_Controls;
    Ui::AtrialFibresClipperViewLabels m_Labels;
    Ui::AtrialFibresViewUIMeshing m_UIMeshing;
    Ui::AtrialFibresClipperViewUIRadius m_UIRadius;
    Ui::AtrialFibresClipperViewUICorridor m_UICorridor;

private:

    void iniPreSurf();
    void Visualiser(double opacity=1.0);
    void VisualiserAuto(double opacity);
    void VisualiserManual(double opacity);
    void VisualisePolyData(vtkSmartPointer<vtkPolyData> pd);
    void VisualiseSphereAtPoint(int ptId, double radius);
    void SphereSourceVisualiser(vtkSmartPointer<vtkPolyData> pointSources, QString colour="1.0,0.0,0.0", double scaleFactor=0.01);
    void PickCallBack(bool pvCorridor=false);
    void ManualCutterCallBack();
    static void KeyCallBackFunc(vtkObject*, long unsigned int, void* ClientData, void*);

    void InitialisePickerObjects();
    void ResetCorridorObjects();

    static QString fileName;
    static QString directory;
    static bool isAutomatic;

    bool automaticPipeline, debugging;

    mitk::Surface::Pointer surface;
    vtkSmartPointer<vtkActor> surfActor;
    std::vector<int> pickedSeedLabels;
    vtkSmartPointer<vtkIdList> pickedSeedIds;
    vtkSmartPointer<vtkPolyData> pickedLineSeeds;
    vtkSmartPointer<vtkPolyData> pickedCutterSeeds;
    vtkSmartPointer<vtkIdList> corridorSeedIds;
    vtkSmartPointer<vtkPolyData> corridorLineSeeds;

    std::unique_ptr<CemrgAtriaClipper> clipper;
    std::unique_ptr<CemrgScarAdvanced> csadv;

    vtkSmartPointer<vtkIdList> pvClipperSeedIdx;
    std::vector<double> pvClipperRadii;

    std::vector<vtkSmartPointer<vtkActor>> clipperActors;

    QDialog* inputs;
    double maxScalar, minScalar, defaultClipperRadius, currentRadius;
    int corridorMax, corridorCount;
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkCallbackCommand> callBack;
    vtkSmartPointer<vtkRenderWindowInteractor> interactor;

};

#endif // AtrialFibresClipperView_h
