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
 * jose.solislemus@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef ScarCalculationsView_h
#define ScarCalculationsView_h

#include <berryISelectionListener.h>
#include <QmitkAbstractView.h>
#include <QMessageBox>
#include <vtkIdList.h>
#include <vtkActor.h>
#include <vtkDijkstraGraphGeodesicPath.h>
#include <string>
#include <sstream>
#include <CemrgAtriaClipper.h>
#include <CemrgScarAdvanced.h>
#include "ui_ScarCalculationsViewControls.h"
#include "ui_ScarCalculationsViewUICorridor.h"
#include "QmitkRenderWindow.h"
#include "mitkCommon.h"
#include "mitkDataStorage.h"
#include "mitkDataNode.h"
#include "mitkSurface.h"
#include "vtkRenderer.h"
#include "vtkTextActor.h"

/**
  \brief ScarCalculationsView

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
class ScarCalculationsView: public QmitkAbstractView {

    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)
    Q_OBJECT

public:

    static const std::string VIEW_ID;
    static void SetDirectoryFile(const QString directory, const QString fileName);
    static void SetCalculationsPaths(const QString directory);
    static void GetInputsFromFile();
    static bool CheckForRequiredFiles();
    static int SearchDirectory(QString searchDir);
    ScarCalculationsView();
    ~ScarCalculationsView();

protected slots:

    void DoImageProcessing(); // F&I T1
    void GapMeasurement(); // F&I T2
    void GapMeasurementVisualisation(const QString &text);
    void BeforeAndAfterComp(); // F7I T3
    void BeforeAndAfterCompVisualisation();
    void Sphericity();
    void CtrlPrePostSelection(const QString &text);
    void SetNewThreshold(const QString &text);
    void EditThreshold();
    void SaveNewThreshold();
    void CancelThresholdEdit();

protected:

    virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer source, const QList<mitk::DataNode::Pointer>& nodes) override;
    virtual void CreateQtPartControl(QWidget *parent) override;
    virtual void SetFocus() override;

    // Helper functions
    void TransformMeshesForComparison();
    void GetThresholdValuesFromFile(QString filepath);
    void SetThresholdValuesToFile(QString filepath);
    void SetShortcutLegend();
    void CopyScalarValues();

    Ui::ScarCalculationsViewControls m_Controls;
    Ui::ScarCalculationsViewUICorridor m_UICorridor;

private:

    void iniPreSurf();
    void Visualiser();
    void BinVisualiser();
    void PickCallBack();

    void InitialisePickerObjects();

    static void KeyCallBackFunc(vtkObject*, long unsigned int, void* ClientData, void*);
    static QStringList CheckForAdvancedDirectoryFiles();

    mitk::Surface::Pointer surface;
    vtkSmartPointer<vtkActor> surfActor;
    vtkSmartPointer<vtkIdList> pickedSeedIds;
    vtkSmartPointer<vtkPolyData> pickedLineSeeds;
    std::unique_ptr<CemrgScarAdvanced> csadv;
    std::vector<vtkSmartPointer<vtkActor> > dijkstraActors;

    QDialog* inputs;
    QString outprefix; // pre, post, tx_postl
    static QString fileName;
    static QString directory;
    static QString predir;
    static QString postdir;
    static QString advdir;
    static QString preScarFile;
    static QString postScarFile;
    int method;
    double value, mean, stdv, thres;
    double maxScalar, minScalar;
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkCallbackCommand> callBack;
    vtkSmartPointer<vtkRenderWindowInteractor> interactor;
};

#endif // ScarCalculationsView_h
