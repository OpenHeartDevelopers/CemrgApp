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

// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>
#include <berryIWorkbenchPage.h>

// Qmitk
#include <mitkImage.h>
#include <mitkIOUtil.h>
#include <mitkProgressBar.h>
#include <mitkNodePredicateProperty.h>
#include <mitkManualSegmentationToSurfaceFilter.h>

#include "FourChamberCommon.h"
#include "FourChamberView.h"
#include "FourChamberGuidePointsView.h"

// VTK
#include <vtkGlyph3D.h>
#include <vtkSphere.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkProperty.h>
#include <vtkCellPicker.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCell.h>
#include <vtkMath.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkFloatArray.h>
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>
#include <vtkCellDataToPointData.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkClipPolyData.h>

//ITK
#include <itkAddImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include "itkLabelObject.h"
#include "itkLabelMap.h"
#include "itkLabelImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelSelectionLabelMapFilter.h"

// Qt
#include <QMessageBox>
#include <QDesktopWidget>

// CemrgAppModule
#include <CemrgAtriaClipper.h>
#include <CemrgCommandLine.h>
#include <CemrgCommonUtils.h>
#include <CemrgFourChamberCmd.h>
#include <FourChamberCommon.h>

QString FourChamberGuidePointsView::fileName;
QString FourChamberGuidePointsView::directory;
QString FourChamberGuidePointsView::whichAtrium;

const std::string FourChamberGuidePointsView::VIEW_ID = "org.mitk.views.fourchamberguidepointsview";

void FourChamberGuidePointsView::CreateQtPartControl(QWidget *parent) {
    MITK_INFO << "FourChamberGuidePointsView";   

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_guide1, SIGNAL(clicked()), this, SLOT(Help()));
    connect(m_Controls.radio_load_la, SIGNAL(toggled(bool)), this, SLOT(LeftAtriumReactToToggle()));
    connect(m_Controls.radio_load_ra, SIGNAL(toggled(bool)), this, SLOT(RightAtriumReactToToggle()));
    connect(m_Controls.button_save, SIGNAL(clicked()), this, SLOT(Save()));
    connect(m_Controls.alpha_slider, SIGNAL(valueChanged(int)), this, SLOT(ChangeAlpha(int)));
    connect(m_Controls.button_opacity, SIGNAL(clicked()), this, SLOT(ConfirmAlpha()));

    inputsSelector = new QDialog(0, Qt::Window);
    m_Selector.setupUi(inputsSelector);
    connect(m_Selector.buttonBox, SIGNAL(accepted()), inputsSelector, SLOT(accept()));
    connect(m_Selector.buttonBox, SIGNAL(rejected()), inputsSelector, SLOT(reject()));

    //Setup renderer
    surfActor = vtkSmartPointer<vtkActor>::New();
    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(0.5,0.5,0.5);
    renderer->AutomaticLightCreationOn();
    renderer->LightFollowCameraOn();
    // renderer->TwoSidedLightingOn();
    // renderer->UpdateLightsGeometryToFollowCamera();
    vtkSmartPointer<vtkTextActor> txtActor = vtkSmartPointer<vtkTextActor>::New();
    std::string shortcuts = GetShortcuts();
    txtActor->SetInput(shortcuts.c_str());
    txtActor->GetTextProperty()->SetFontSize(14);
    txtActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
    renderer->AddActor2D(txtActor);

    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow =
            vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    m_Controls.widget_1->setRenderWindow(renderWindow);
    m_Controls.widget_1->renderWindow()->AddRenderer(renderer);

    //Setup keyboard interactor
    callBack = vtkSmartPointer<vtkCallbackCommand>::New();
    callBack->SetCallback(KeyCallBackFunc);
    callBack->SetClientData(this);
    interactor = m_Controls.widget_1->renderWindow()->GetInteractor();
    interactor->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New());
    interactor->GetInteractorStyle()->KeyPressActivationOff();
    interactor->GetInteractorStyle()->AddObserver(vtkCommand::KeyPressEvent, callBack);

    //Initialisation
    m_Controls.radio_load_la->setChecked(true);
    m_Selector.radioBtn_ATRIAL_SEPTUM->setEnabled(false); // Not yet implemented
    iniPreSurf();
    if (surface.IsNotNull()) {
        InitialisePickerObjects();
        Visualiser();
    }

    Help(true);
    alpha = 1.0;
    pluginLoaded = true;
}

void FourChamberGuidePointsView::SetFocus() {
    m_Controls.button_guide1->setFocus();
}

void FourChamberGuidePointsView::OnSelectionChanged(
        berry::IWorkbenchPart::Pointer /*src*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

FourChamberGuidePointsView::~FourChamberGuidePointsView() {
    inputsSelector->deleteLater();
}

// slots
void FourChamberGuidePointsView::Help(bool firstTime){
    std::string msg = "";
    if(firstTime){
        msg = "HELP\n";
        msg += "1. Check the boxes to load the left or right atrium.\n";
        msg += "2. Use the slider to adjust the opacity of the loaded atrium.\n";
        msg += "3. Hover on top of the surface and press SPACE to select points.\n";
    } else {
        msg = "HELP\n";
        msg += "1. Check the boxes to load the left or right atrium.\n";
        msg += "2. Use the slider to adjust the opacity of the loaded atrium.\n";
        msg += "3. Hover on top of the surface and press SPACE to select points.\n";
    }
    QMessageBox::information(NULL, "Help", msg.c_str());
}

void FourChamberGuidePointsView::SetDirectoryFile(const QString directory) {
    FourChamberGuidePointsView::directory = directory;
}

void FourChamberGuidePointsView::SetSubdirs(){
    QString dir = FourChamberGuidePointsView::directory;
    parfiles = dir + "/parfiles";
    apex_septum_files = parfiles + "/apex_septum_templates";
    surfaces_uvc_la = dir + "/surfaces_uvc_LA";
    surfaces_uvc_ra = dir + "/surfaces_uvc_RA";

    path_to_la = surfaces_uvc_la + "/la/la_vis.surfmesh.vtk";
    path_to_ra = surfaces_uvc_ra + "/ra/ra_vis.surfmesh.vtk";
}

bool FourChamberGuidePointsView::CreateVisualisationMesh(QString dir, QString visName, QString originalName) {
    MITK_INFO << "[INFO] Creating visualisation mesh";
    bool success = false;
    QString path = dir + "/" + visName + ".surfmesh.vtk";
    QFileInfo fi(path);
    if (!fi.exists()) {
        MITK_INFO << "[INFO] Visualisation mesh does not exist. Creating...";
        QString originalPath = dir + "/" + originalName;
        if (!QFile::exists(originalPath)) {
            MITK_WARN << "File" << originalPath.toStdString() << "does not exist!";
            return false;
        }

        std::unique_ptr<CemrgFourChamberCmd> fourch = std::make_unique<CemrgFourChamberCmd>();
        QStringList arguments;
        arguments << "-msh=" + originalName;
        arguments << "-surf=" + visName;
        arguments << "-ofmt=vtk_polydata";
        QString path = fourch->DockerMeshtoolGeneric(dir, "extract", "surface", arguments, visName);

        success = (QFile::exists(fi.absoluteFilePath()));
    } else {
        MITK_INFO << "[INFO] Visualisation mesh already exists!";
        success = true;
    }

    return success;
}

void FourChamberGuidePointsView::iniPreSurf() {
    MITK_INFO << "[INFO] Initialising pre-surface visualisation";
    SetSubdirs();
    bool laVisSuccess = CreateVisualisationMesh(surfaces_uvc_la, QString("la/la_vis"), QString("la/la.vtk"));
    bool raVisSuccess = CreateVisualisationMesh(surfaces_uvc_ra, QString("ra/ra_vis"), QString("ra/ra.vtk"));

    m_Controls.radio_load_la->setEnabled(laVisSuccess);
    m_Controls.radio_load_ra->setEnabled(raVisSuccess);

    if (!laVisSuccess && !raVisSuccess) {
        MITK_WARN << "No visualisation files found!";
        return;
    }

    // Find the selected node
    QString path = m_Controls.radio_load_la->isChecked() ? path_to_la : path_to_ra;
    mitk::Surface::Pointer shell = mitk::IOUtil::Load<mitk::Surface>(path.toStdString());
    if (shell.IsNull()) {
        MITK_WARN << "Shell is null!";
        return;
    }
    surface = shell;
}

void FourChamberGuidePointsView::Visualiser(double opacity){
    MITK_INFO << "[Visualiser]";
    double max_scalar = 4; // RA
    double min_scalar = 3; // LA

    renderer->RemoveAllViewProps();

    SphereSourceVisualiser(pickedPointsHandler->GetLineSeeds(), "0.4,0.1,0.0", 0.025);

    //Create a mapper and actor for surface
    vtkSmartPointer<vtkPolyDataMapper> surfMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    surfMapper->SetInputData(surface->GetVtkPolyData());
    surfMapper->SetScalarRange(min_scalar, max_scalar);
    surfMapper->SetScalarModeToUseCellData();
    surfMapper->ScalarVisibilityOn();

    vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
    lut->SetTableRange(min_scalar, max_scalar);
    lut->SetHueRange(0.5, 0.0);  // this is the way_neighbourhood_size you tell which colors you want to be displayed.
    lut->Build();     // this is important

    surfMapper->SetLookupTable(lut);

    vtkSmartPointer<vtkActor> surfActor = vtkSmartPointer<vtkActor>::New();
    surfActor->SetMapper(surfMapper);
    surfActor->GetProperty()->SetOpacity(opacity);
    renderer->AddActor(surfActor);
}

void FourChamberGuidePointsView::SphereSourceVisualiser(vtkSmartPointer<vtkPolyData> pointSources, QString colour, double scaleFactor){
    MITK_INFO << "[INFO] Sphere source visualiser"; 
    // e.g colour = "0.4,0.1,0.0" - values for r,g, and b separated by commas.
    double r, g, b;
    bool okr, okg, okb;
    r = colour.section(',',0,0).toDouble(&okr);
    g = colour.section(',',1,1).toDouble(&okg);
    b = colour.section(',',2,2).toDouble(&okb);

    if(!okr) r=1.0;
    if(!okg) g=0.0;
    if(!okb) b=0.0;

    vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
    vtkSmartPointer<vtkSphereSource> glyphSource = vtkSmartPointer<vtkSphereSource>::New();
    glyph3D->SetInputData(pointSources);
    glyph3D->SetSourceConnection(glyphSource->GetOutputPort());
    glyph3D->SetScaleModeToDataScalingOff();
    glyph3D->SetScaleFactor(surface->GetVtkPolyData()->GetLength()*scaleFactor);
    glyph3D->Update();

    //Create a mapper and actor for glyph
    vtkSmartPointer<vtkPolyDataMapper> glyphMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    glyphMapper->SetInputConnection(glyph3D->GetOutputPort());
    vtkSmartPointer<vtkActor> glyphActor = vtkSmartPointer<vtkActor>::New();
    glyphActor->SetMapper(glyphMapper);
    glyphActor->GetProperty()->SetColor(r,g,b);
    glyphActor->PickableOff();
    renderer->AddActor(glyphActor);
}

void FourChamberGuidePointsView::PickCallBack() {

    vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
    picker->SetTolerance(1E-4 * surface->GetVtkPolyData()->GetLength());
    int* eventPosition = interactor->GetEventPosition();
    int result = picker->Pick(float(eventPosition[0]), float(eventPosition[1]), 0.0, renderer);
    if (result == 0) return;
    double* pickPosition = picker->GetPickPosition();
    std::cout << "Pick position: " << pickPosition[0] << " " << pickPosition[1] << " " << pickPosition[2] << std::endl;
    vtkIdList* pickedCellPointIds = surface->GetVtkPolyData()->GetCell(picker->GetCellId())->GetPointIds();

    double distance;
    int pickedSeedId = -1;
    double minDistance = 1E10;
    for (int i=0; i<pickedCellPointIds->GetNumberOfIds(); i++) {
        distance = vtkMath::Distance2BetweenPoints(
                    pickPosition, surface->GetVtkPolyData()->GetPoint(pickedCellPointIds->GetId(i)));
        std::cout << "Point ID: " << pickedCellPointIds->GetId(i) << std::endl;
        std::cout << "Distance: " << distance << std::endl; 
        if (distance < minDistance) {
            minDistance = distance;
            pickedSeedId = pickedCellPointIds->GetId(i);
        }//_if
    }//_for
    std::cout << "Picked seed ID: " << pickedSeedId << std::endl;
    std::cout << "Picked seed distance: " << minDistance << std::endl;
    if (pickedSeedId == -1){
        std::cout << "No point was picked!" << std::endl;
        pickedSeedId = pickedCellPointIds->GetId(0);
    }
    
    // pushedLabel gets updated in UserSelectPvLabel function
    pickedPointsHandler->AddPointFromSurface(surface, pickedSeedId, pushedLabel);
    MITK_INFO << pickedPointsHandler->ToString();

    m_Controls.widget_1->renderWindow()->Render();
}

void FourChamberGuidePointsView::KeyCallBackFunc(vtkObject*, long unsigned int, void* ClientData, void*) {

    FourChamberGuidePointsView* self;
    self = reinterpret_cast<FourChamberGuidePointsView*>(ClientData);
    std::string key = self->interactor->GetKeySym();

    if (!self->m_Controls.radio_load_la->isChecked() && !self->m_Controls.radio_load_ra->isChecked()) {
        QMessageBox::information(NULL, "Info", "Please load a mesh first!");
        return;
    }

    if (self->m_Controls.radio_load_la->isChecked() && self->m_Controls.radio_load_ra->isChecked()) { 
        QMessageBox::information(NULL, "Info", "Please load only one mesh at a time!");
        return;
    }

    if (key == "space") {
        //Ask the labels
        self->UserSelectPvLabel();
        self->PickCallBack();

    } else if (key == "Delete") {

        int lastLabel = self->pickedPointsHandler->CleanupLastPoint();

        if (lastLabel != -1) {
            if (lastLabel == AtrialLandmarksType::LA_APEX)
                self->m_Selector.radioBtn_LA_APEX->setEnabled(true);
            else if (lastLabel == AtrialLandmarksType::LA_SEPTUM)
                self->m_Selector.radioBtn_LA_SEPTUM->setEnabled(true);
            else if (lastLabel == AtrialLandmarksType::RA_APEX)
                self->m_Selector.radioBtn_RA_APEX->setEnabled(true);
            else if (lastLabel == AtrialLandmarksType::RA_SEPTUM)
                self->m_Selector.radioBtn_RA_SEPTUM->setEnabled(true);
            else if (lastLabel == AtrialLandmarksType::RAA_APEX)
                self->m_Selector.radioBtn_RAA_APEX->setEnabled(true);
            else if (lastLabel == AtrialLandmarksType::ATRIAL_SEPTUM) 
                self->m_Selector.radioBtn_ATRIAL_SEPTUM->setEnabled(true);
            
        }//_if

        self->m_Controls.widget_1->renderWindow()->Render();
    } else if (key == "H" || key == "h"){
        self->Help();
    }
}

// helper functions
void FourChamberGuidePointsView::InitialisePickerObjects(){
    pickedPointsHandler = std::unique_ptr<CemrgMeshPointSelector>(new CemrgMeshPointSelector());
    QStringList names = QStringList();
    names << "la.lvapex.vtx" << "la.rvsept_pt.vtx" << "ra.lvapex.vtx" << "ra.rvsept_pt.vtx" << "raa_apex.txt";
    std::vector<int> labels = {11, 13, 15, 17, 19};

    pickedPointsHandler->SetAvailableLabels(names, labels);

}

std::string FourChamberGuidePointsView::GetShortcuts(){
    std::string res = "";
    
    res += " POINT SELECTION:\n\tSpace: select location\n\tDelete: remove location";
    res += "\nHELP:\n\tH/h: Guides";

    return res;
}

void FourChamberGuidePointsView::LeftAtriumReactToToggle() {
    if (pluginLoaded && m_Controls.radio_load_la->isChecked()) {
        QMessageBox::information(NULL, "Info", "Left Atrium selected");
        iniPreSurf();
        Visualiser();
    }
}

void FourChamberGuidePointsView::RightAtriumReactToToggle() {
    if (pluginLoaded && m_Controls.radio_load_ra->isChecked()) {
        QMessageBox::information(NULL, "Info", "Right Atrium selected");
        iniPreSurf();
        Visualiser();
    }
}

void FourChamberGuidePointsView::ChangeAlpha(int alphaChange) {
    alpha = alphaChange / 100.0;
}

void FourChamberGuidePointsView::ConfirmAlpha() {
    // LoadMeshes(2);
    Visualiser(alpha);
}

void FourChamberGuidePointsView::Save() {
    QMessageBox::information(NULL, "Info", "Saving points to files");
    QStringList saveTypes = {"vtx", "coord"};
    QStringList laPoints = {"la.lvapex.vtx", "la.rvsept_pt.vtx"}; // LA points
    QStringList raPoints = {"ra.lvapex.vtx", "ra.rvsept_pt.vtx", "raa_apex.txt"}; // RA points

    std::string la_path = surfaces_uvc_la.toStdString() + "/la/la.vtk";
    mitk::UnstructuredGrid::Pointer grid = mitk::IOUtil::Load<mitk::UnstructuredGrid>(la_path);
    if (grid.IsNull()) {
        MITK_WARN << "LA Grid is null!";
        return;
    }
    pickedPointsHandler->UpdateVtkIdFromUnstructuredGrid(grid, laPoints);

    std::string ra_path = surfaces_uvc_ra.toStdString() + "/ra/ra.vtk";
    grid = mitk::IOUtil::Load<mitk::UnstructuredGrid>(ra_path);
    if (grid.IsNull()) {
        MITK_WARN << "RA Grid is null!";
        return;
    }
    pickedPointsHandler->UpdateVtkIdFromUnstructuredGrid(grid, raPoints);

    QStringList points = laPoints + raPoints;

    // Save points
    for (auto& point : points) {
        bool isVtxType = pickedPointsHandler->GetPointData(point).pointName.contains("vtx");
        QString savetype = isVtxType ? "vtx" : "coord";
        pickedPointsHandler->SaveToFile(APEX_SEPTUM(point), (QStringList(point)), savetype);
        if (isVtxType) { 
            pickedPointsHandler->SaveToFile(APEX_SEPTUM(point), (QStringList(point.replace(".vtx", "_coord.txt"))), "coord");
        }
    }
    
    QMessageBox::information(NULL, "Info", "Points saved, please close this window.");
    m_Controls.button_save->setEnabled(false);
}

void FourChamberGuidePointsView::UserSelectPvLabel(){
    int dialogCode = inputsSelector->exec();
    QRect screenGeometry = QApplication::desktop()->screenGeometry();
    int x = (screenGeometry.width() - inputsSelector->width()) / 2;
    int y = (screenGeometry.height() - inputsSelector->height()) / 2;
    inputsSelector->move(x,y);

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {

        if (m_Selector.radioBtn_LA_APEX->isChecked()) {
            pushedLabel = AtrialLandmarksType::LA_APEX; // LA_APEX
            m_Selector.radioBtn_LA_APEX->setEnabled(false);
        } else if (m_Selector.radioBtn_LA_SEPTUM->isChecked()) {
            pushedLabel = AtrialLandmarksType::LA_SEPTUM; // LA_SEPTUM
            m_Selector.radioBtn_LA_SEPTUM->setEnabled(false);
        } else if (m_Selector.radioBtn_RA_APEX->isChecked()) {
            pushedLabel = AtrialLandmarksType::RA_APEX; // RA_APEX
            m_Selector.radioBtn_RA_APEX->setEnabled(false);
        } else if (m_Selector.radioBtn_RA_SEPTUM->isChecked()) {
            pushedLabel = AtrialLandmarksType::RA_SEPTUM; // RA_SEPTUM
            m_Selector.radioBtn_RA_SEPTUM->setEnabled(false);
        } else if (m_Selector.radioBtn_RAA_APEX->isChecked()) {
            pushedLabel = AtrialLandmarksType::RAA_APEX; // RAA_APEX
            m_Selector.radioBtn_RAA_APEX->setEnabled(false);
        } else if (m_Selector.radioBtn_ATRIAL_SEPTUM->isChecked()) {
            pushedLabel = AtrialLandmarksType::ATRIAL_SEPTUM; // ATRIAL_SEPTUM
            m_Selector.radioBtn_ATRIAL_SEPTUM->setEnabled(false);
        }
        MITK_INFO << "Label selected: " << pushedLabel;
        
    } else if (dialogCode == QDialog::Rejected) {
        inputsSelector->close();
    }//_if
}


/*
========================
 CemrgApp radiobtn codes
========================
=== Landmarks ===
SVC_POST - radioBtn_RA_SVC_POST - 29
IVC_POST - radioBtn_RA_IVC_POST - 31
RAA_VALVE_P - radioBtn_RAA_TCV -  33
CS_VALVE_P - radioBtn_RA_CS_TCV - 35
SVC_ANT  - radioBtn_RA_SVC_ANT -  37
IVC_ANT  - radioBtn_RA_IVC_ANT -  39

=== Region ===
IVC_ANT - radioBtn_RA_IVC_ANT -      29
CS_TOP - radioBtn_RA_CS -            31
IVC_SVC_ANT - radioBtn_RA_IvcSvc -   33
SVC_ANT - radioBtn_RA_SVC_ANT -      35
RAA_ANT - radioBtn_RAA_ANT -         37
RAA_CS_ANT - radioBtn_RAA_CS_ANT -   39
*/
