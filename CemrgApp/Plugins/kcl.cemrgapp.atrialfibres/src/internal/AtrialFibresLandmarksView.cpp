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
#include "AtrialFibresLandmarksView.h"
#include "AtrialFibresView.h"

// VTK
#include <vtkGlyph3D.h>
#include <vtkSphere.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
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

QString AtrialFibresLandmarksView::fileName;
QString AtrialFibresLandmarksView::directory;
QString AtrialFibresLandmarksView::whichAtrium;

const std::string AtrialFibresLandmarksView::VIEW_ID = "org.mitk.views.atrialfibreslandmarksview";

void AtrialFibresLandmarksView::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_guide1, SIGNAL(clicked()), this, SLOT(Help()));
    connect(m_Controls.button_save1_rough, SIGNAL(clicked()), this, SLOT(SaveRoughPoints()));
    connect(m_Controls.button_save2_refined, SIGNAL(clicked()), this, SLOT(SaveRefinedPoints()));

    isLeftAtrium = (AtrialFibresLandmarksView::whichAtrium.compare("LA", Qt::CaseInsensitive)==0);
    std::cout << (isLeftAtrium ? "Working on Left Atrium" : "Working on Right Atrium") << '\n';

    //Create GUI widgets
    inputsRough = new QDialog(0,0);
    m_Rough.setupUi(inputsRough);
    connect(m_Rough.buttonBox, SIGNAL(accepted()), inputsRough, SLOT(accept()));
    connect(m_Rough.buttonBox, SIGNAL(rejected()), inputsRough, SLOT(reject()));

    inputsRefined = new QDialog(0,0);
    m_Refined.setupUi(inputsRefined);
    connect(m_Refined.buttonBox, SIGNAL(accepted()), inputsRefined, SLOT(accept()));
    connect(m_Refined.buttonBox, SIGNAL(rejected()), inputsRefined, SLOT(reject()));

    RoughUiEnableButtons();
    RefinedUiEnableButtons();

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
    m_Controls.widget_1->SetRenderWindow(renderWindow);
    m_Controls.widget_1->GetRenderWindow()->AddRenderer(renderer);

    //Setup keyboard interactor
    callBack = vtkSmartPointer<vtkCallbackCommand>::New();
    callBack->SetCallback(KeyCallBackFunc);
    callBack->SetClientData(this);
    interactor = m_Controls.widget_1->GetRenderWindow()->GetInteractor();
    interactor->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New());
    interactor->GetInteractorStyle()->KeyPressActivationOff();
    interactor->GetInteractorStyle()->AddObserver(vtkCommand::KeyPressEvent, callBack);

    //Initialisation
    iniPreSurf();
    if (surface.IsNotNull()) {
        InitialisePickerObjects();
        Visualiser();
    }

    m_Controls.button_save2_refined->setEnabled(false);
    Help(true);
}

void AtrialFibresLandmarksView::SetFocus() {
    m_Controls.button_guide1->setFocus();
}

void AtrialFibresLandmarksView::OnSelectionChanged(
        berry::IWorkbenchPart::Pointer /*src*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

AtrialFibresLandmarksView::~AtrialFibresLandmarksView() {
    inputsRough->deleteLater();
    inputsRefined->deleteLater();
}

// slots
void AtrialFibresLandmarksView::Help(bool firstTime){
    std::string msg = "";
    if(firstTime){
        msg = "HELP\n Select Rough locations with the spacebar, ";
        msg += "then click the Save Rough Locations button to enable saving ";
        msg += "Refined Loactions.\n\n";
        msg += GetRoughPointsGuide();
    } else if(!m_Controls.button_save2_refined->isEnabled()) {
        msg = GetRoughPointsGuide();
    } else {
        msg = GetRefinedPointsGiude();
    }
    QMessageBox::information(NULL, "Help", msg.c_str());
}

void AtrialFibresLandmarksView::SaveRoughPoints(){
    if(roughSeedIds->GetNumberOfIds() == 0){
        QMessageBox::warning(NULL, "Rough landmarks not selected.", "Select the correct rough landmarks.");
        return;
    }

    MITK_INFO << "[SaveRoughPoints] Saving rough points to file.";
    QString prodPath = directory + "/";
    QString outname = (isLeftAtrium) ? "prodLaRoughLandmarks" : "prodRaLandmarks";
    ofstream fileRough, fileRoughLabels;

    MITK_INFO << "[SaveRoughPoints] Saving TXT file.";
    fileRough.open((prodPath + outname + ".txt").toStdString());
    fileRoughLabels.open((prodPath + outname + "-Labels.txt").toStdString());

	std::vector<int> roughPointsOrder;
	if (isLeftAtrium) {
		roughPointsOrder = { 15, 17, 13, 11, 21, 19 };
	} else {
		roughPointsOrder = { 29, 31, 33, 35, 37, 39 };
	}
    for (unsigned int ix = 0; ix<roughPointsOrder.size(); ix++) {
        int index = GetIndex(roughSeedLabels, roughPointsOrder.at(ix));
        if(index!=-1){
            std::cout << GetStructureIdFromLabel(false, roughSeedLabels.at(index)) << '\n';

            vtkIdType vId = roughSeedIds->GetId(index);
            double* point = surface->GetVtkPolyData()->GetPoint(vId);

            fileRough << std::setprecision(12) << point[0] << "," << point[1] << "," << point[2] << "\n";
            fileRoughLabels << GetStructureIdFromLabel(false, roughSeedLabels.at(index)) << "\n";
        } else{
            MITK_WARN << "[SaveRoughPoints] Value not found";
        }
    }
    fileRough.close();

    m_Controls.button_save1_rough->setEnabled(false);
    m_Controls.button_save2_refined->setEnabled(true);
    Help();


}

void AtrialFibresLandmarksView::SaveRefinedPoints(){
    if(refinedSeedIds->GetNumberOfIds() == 0){
        QMessageBox::warning(NULL, "Refined landmarks not selected.", "Select the correct refined landmarks.");
        return;
    }

    MITK_INFO << "[SaveRefinedPoints] Saving refined points to file.";
    QString prodPath = directory + "/";
    QString outname = (isLeftAtrium) ? "prodLaRefinedLandmarks" : "prodRaRegion";
    ofstream fileRefined;
    ofstream fileRough, fileRefinedLabels;

    MITK_INFO << "[SaveRefinedPoints] Saving TXT file.";
    fileRefined.open((prodPath + outname + ".txt").toStdString());
    fileRefinedLabels.open((prodPath + outname + "-Labels.txt").toStdString());

	std::vector<int> refinedPointsOrder;
	if (isLeftAtrium) {
		refinedPointsOrder = { 19, 22, 13, 17 };
	} else {
		refinedPointsOrder = { 29, 31, 33, 35, 37, 39 };
	}
    for (unsigned int ix = 0; ix<refinedPointsOrder.size(); ix++) {
        int index = GetIndex(refinedSeedLabels, refinedPointsOrder.at(ix));
        if(index!=-1){
            std::cout << GetStructureIdFromLabel(true, refinedSeedLabels.at(index)) << '\n';
            vtkIdType vId = refinedSeedIds->GetId(index);
            double* point = surface->GetVtkPolyData()->GetPoint(vId);

            fileRefined << std::setprecision(12) << point[0] << "," << point[1] << "," << point[2] << "\n";
            fileRefinedLabels << GetStructureIdFromLabel(true, refinedSeedLabels.at(index)) << "\n";

        } else{
            MITK_WARN << "[SaveRefinedPoints] Value not found";
        }
    }
    fileRefined.close();

    m_Controls.button_save2_refined->setEnabled(false);
    m_Controls.button_guide1->setEnabled(false);
}

void AtrialFibresLandmarksView::SetDirectoryFile(const QString directory, const QString fileName, const QString whichAtrium) {
    AtrialFibresLandmarksView::fileName = fileName;
    AtrialFibresLandmarksView::directory = directory;
    AtrialFibresLandmarksView::whichAtrium = whichAtrium;
}

void AtrialFibresLandmarksView::iniPreSurf() {
    //Find the selected node
    QString path = AtrialFibresLandmarksView::directory + "/" + AtrialFibresLandmarksView::fileName;
    // mitk::Surface::Pointer shell = mitk::IOUtil::Load<mitk::Surface>(path.toStdString());
    mitk::Surface::Pointer shell = CemrgCommonUtils::LoadVTKMesh(path.toStdString());
    CemrgCommonUtils::FlipXYPlane(shell, "", "");
    mitk::IOUtil::Save(shell, path.toStdString());

    surface = shell;
}

void AtrialFibresLandmarksView::Visualiser(double opacity){
    MITK_INFO << "[Visualiser]";
    double max_scalar = (isLeftAtrium) ? 19 : 7;
    double min_scalar = 1;
    vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
    // double max_scalar=-1, min_scalar=1e9,s;
    // vtkFloatArray *scalars = vtkFloatArray::New();
    // scalars = vtkFloatArray::SafeDownCast(surface->GetVtkPolyData()->GetCellData()->GetScalars());
    // for (vtkIdType i=0;i<surface->GetVtkPolyData()->GetNumberOfCells();i++) {
    //     s = scalars->GetTuple1(i);
    //     if (s > max_scalar)
    //         max_scalar = s;
    //     if (s < min_scalar)
    //         min_scalar = s;
    // }
    this->maxScalar = max_scalar;
    this->minScalar = min_scalar;

    SphereSourceVisualiser(roughLineSeeds);
    SphereSourceVisualiser(refinedLineSeeds, "0.0,0.0,1.0");

    //Create a mapper and actor for surface
    vtkSmartPointer<vtkPolyDataMapper> surfMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    surfMapper->SetInputData(surface->GetVtkPolyData());
    surfMapper->SetScalarRange(min_scalar, max_scalar);
    surfMapper->SetScalarModeToUseCellData();
    surfMapper->ScalarVisibilityOn();
    lut->SetTableRange(min_scalar, max_scalar);
    lut->SetHueRange(0.5, 0.0);  // this is the way_neighbourhood_size you tell which colors you want to be displayed.
    lut->Build();     // this is important

    vtkSmartPointer<vtkScalarBarActor> scalarBar =
            vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBar->SetLookupTable(surfMapper->GetLookupTable());
    scalarBar->SetWidth(0.075);
    scalarBar->SetHeight(0.3);
    scalarBar->SetTextPositionToPrecedeScalarBar();
    scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    scalarBar->GetPositionCoordinate()->SetValue( 0.9, 0.01 );
    scalarBar->SetNumberOfLabels(5);

    std::string scalarBarTitle = "Atrium Labels";
    scalarBar->SetTitle(scalarBarTitle.c_str());
    scalarBar->SetVerticalTitleSeparation(15);

    surfMapper->SetLookupTable(lut);
    scalarBar->SetLookupTable(lut);

    vtkSmartPointer<vtkActor> surfActor = vtkSmartPointer<vtkActor>::New();
    surfActor->SetMapper(surfMapper);
    surfActor->GetProperty()->SetOpacity(opacity);
    renderer->AddActor(surfActor);
    renderer->AddActor2D(scalarBar);
}

void AtrialFibresLandmarksView::SphereSourceVisualiser(vtkSmartPointer<vtkPolyData> pointSources, QString colour, double scaleFactor){
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

void AtrialFibresLandmarksView::PickCallBack(bool refinedLandmarks) {

    vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
    picker->SetTolerance(1E-4 * surface->GetVtkPolyData()->GetLength());
    int* eventPosition = interactor->GetEventPosition();
    int result = picker->Pick(float(eventPosition[0]), float(eventPosition[1]), 0.0, renderer);
    if (result == 0) return;
    double* pickPosition = picker->GetPickPosition();
    vtkIdList* pickedCellPointIds = surface->GetVtkPolyData()->GetCell(picker->GetCellId())->GetPointIds();

    double distance;
    int pickedSeedId = -1;
    double minDistance = 1E10;
    for (int i=0; i<pickedCellPointIds->GetNumberOfIds(); i++) {
        distance = vtkMath::Distance2BetweenPoints(
                    pickPosition, surface->GetVtkPolyData()->GetPoint(pickedCellPointIds->GetId(i)));
        if (distance < minDistance) {
            minDistance = distance;
            pickedSeedId = pickedCellPointIds->GetId(i);
        }//_if
    }//_for
    if (pickedSeedId == -1){
        pickedSeedId = pickedCellPointIds->GetId(0);
    }

    double* point = surface->GetVtkPolyData()->GetPoint(pickedSeedId);
    if(!refinedLandmarks){
        roughSeedIds->InsertNextId(pickedSeedId);
        roughLineSeeds->GetPoints()->InsertNextPoint(point);
        roughLineSeeds->Modified();
    } else{
        refinedSeedIds->InsertNextId(pickedSeedId);
        refinedLineSeeds->GetPoints()->InsertNextPoint(point);
        refinedLineSeeds->Modified();
    }

    m_Controls.widget_1->GetRenderWindow()->Render();
}

void AtrialFibresLandmarksView::KeyCallBackFunc(vtkObject*, long unsigned int, void* ClientData, void*) {

    AtrialFibresLandmarksView* self;
    self = reinterpret_cast<AtrialFibresLandmarksView*>(ClientData);
    std::string key = self->interactor->GetKeySym();

    if (key == "space") {
        //Ask the labels
        self->PickCallBack();
        self->UserSelectPvLabel();

    } else if (key == "Delete") {

        //Clean up last dropped seed point
        vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkPoints> points = self->roughLineSeeds->GetPoints();
        for (int i=0; i<points->GetNumberOfPoints()-1; i++){
            newPoints->InsertNextPoint(points->GetPoint(i));
        }
        self->roughLineSeeds->SetPoints(newPoints);
        vtkSmartPointer<vtkIdList> newSeedIds = vtkSmartPointer<vtkIdList>::New();
        newSeedIds->Initialize();
        vtkSmartPointer<vtkIdList> roughSeedIds = self->roughSeedIds;
        for (int i=0; i<roughSeedIds->GetNumberOfIds()-1; i++){
            newSeedIds->InsertNextId(roughSeedIds->GetId(i));
        }
        self->roughSeedIds = newSeedIds;

        if (self->roughSeedLabels.empty() == false) {
            int radioButtonNumber = self->roughSeedLabels.back() - 10;
            if (radioButtonNumber == 1)
            self->m_Rough.radioBtn_LA_LSPV->setEnabled(true);
            else if (radioButtonNumber == 3)
            self->m_Rough.radioBtn_LA_LIPV->setEnabled(true);
            else if (radioButtonNumber == 5)
            self->m_Rough.radioBtn_LA_RSPV->setEnabled(true);
            else if (radioButtonNumber == 7)
            self->m_Rough.radioBtn_LA_RIPV->setEnabled(true);
            else if (radioButtonNumber == 9)
            self->m_Rough.radioBtn_LAA_base->setEnabled(true);
            else if (radioButtonNumber == 11)
            self->m_Rough.radioBtn_LAA_tip->setEnabled(true);
            else if (radioButtonNumber == 19)
            self->m_Rough.radioBtn_RA_SVC_POST->setEnabled(true);
            else if (radioButtonNumber == 21)
            self->m_Rough.radioBtn_RA_IVC_POST->setEnabled(true);
            else if (radioButtonNumber == 23)
            self->m_Rough.radioBtn_RAA_TCV->setEnabled(true);
            else if (radioButtonNumber == 25)
            self->m_Rough.radioBtn_RA_CS_TCV->setEnabled(true);
            else if (radioButtonNumber == 27)
            self->m_Rough.radioBtn_RA_SVC_ANT->setEnabled(true);
            else if (radioButtonNumber == 29)
            self->m_Rough.radioBtn_RA_IVC_ANT->setEnabled(true);

            self->roughSeedLabels.pop_back();
        }//_if

        self->m_Controls.widget_1->GetRenderWindow()->Render();
    } else if (key == "X" || key == "x"){
        if(self->m_Controls.button_save2_refined->isEnabled()){
            bool refinedLandmarks = true;
            self->PickCallBack(refinedLandmarks);
            self->UserSelectPvLabel(refinedLandmarks);
        }
    } else if (key == "D" || key == "d"){
        if(self->m_Controls.button_save2_refined->isEnabled()){
            //Clean up last dropped seed point
            vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
            vtkSmartPointer<vtkPoints> points = self->refinedLineSeeds->GetPoints();
            for (int i=0; i<points->GetNumberOfPoints()-1; i++){
                newPoints->InsertNextPoint(points->GetPoint(i));
            }
            self->refinedLineSeeds->SetPoints(newPoints);
            vtkSmartPointer<vtkIdList> newSeedIds = vtkSmartPointer<vtkIdList>::New();
            newSeedIds->Initialize();
            vtkSmartPointer<vtkIdList> roughSeedIds = self->refinedSeedIds;
            for (int i=0; i<roughSeedIds->GetNumberOfIds()-1; i++){
                newSeedIds->InsertNextId(roughSeedIds->GetId(i));
            }
            self->refinedSeedIds = newSeedIds;

            if (self->refinedSeedLabels.empty() == false) {
                int radioButtonNumber = self->refinedSeedLabels.back() - 10;
                if (radioButtonNumber == 1)
                self->m_Refined.radioBtn_LA_LSPV->setEnabled(true);
                else if (radioButtonNumber == 3)
                self->m_Refined.radioBtn_LA_LspvBody->setEnabled(true);
                else if (radioButtonNumber == 5)
                self->m_Refined.radioBtn_LA_RSPV->setEnabled(true);
                else if (radioButtonNumber == 7)
                self->m_Refined.radioBtn_LA_RspvBody->setEnabled(true);
                else if (radioButtonNumber == 9)
                self->m_Refined.radioBtn_LA_LatWall->setEnabled(true);
                else if (radioButtonNumber == 12)
                self->m_Refined.radioBtn_LA_FO->setEnabled(true);

                else if (radioButtonNumber == 19)
                self->m_Refined.radioBtn_RA_IVC_ANT->setEnabled(true);
                else if (radioButtonNumber == 21)
                self->m_Refined.radioBtn_RA_CS->setEnabled(true);
                else if (radioButtonNumber == 23)
                self->m_Refined.radioBtn_RA_IvcSvc->setEnabled(true);
                else if (radioButtonNumber == 25)
                self->m_Refined.radioBtn_RA_SVC_ANT->setEnabled(true);
                else if (radioButtonNumber == 27)
                self->m_Refined.radioBtn_RAA_ANT->setEnabled(true);
                else if (radioButtonNumber == 29)
                self->m_Refined.radioBtn_RAA_CS_ANT->setEnabled(true);

                self->refinedSeedLabels.pop_back();
            }//_if

            self->m_Controls.widget_1->GetRenderWindow()->Render();
        }
    } else if (key == "H" || key == "h"){
        self->Help();
    }
}

// helper functions
void AtrialFibresLandmarksView::InitialisePickerObjects(){
    roughSeedIds = vtkSmartPointer<vtkIdList>::New();
    roughSeedIds->Initialize();
    roughLineSeeds = vtkSmartPointer<vtkPolyData>::New();
    roughLineSeeds->Initialize();
    roughLineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());


    refinedSeedIds = vtkSmartPointer<vtkIdList>::New();
    refinedSeedIds->Initialize();
    refinedLineSeeds = vtkSmartPointer<vtkPolyData>::New();
    refinedLineSeeds->Initialize();
    refinedLineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
}

std::string AtrialFibresLandmarksView::GetShortcuts(){
    std::string res = "";
    res += "ROUGH POINT SELECTION:\n\tSpace: select rough location\n\tDelete: remove rough location";
    res += "\n\nREFINED POINT SELECTION:\n\tX: Select refined landmark\n\tD: remove refined landmark";
    res += "\nHELP:\n\tH/h: Guides";

    return res;
}

std::string AtrialFibresLandmarksView::GetRoughPointsGuide(){
    std::string res = "ROUGH LANDMARKS GUIDE\n Select rough locations for:\n";
    if(isLeftAtrium){
        res += "LSPV, LIPV\n RSPV, RIPV\nLAA tip, and LAA base.\n";
    } else{
        res += "SVC posterior \n";
        res += "IVC posterior \n";
        res += "RAA/TCV posterior \n";
        res += "CS/TCV posterior \n";
        res += "SVC anterior \n";
        res += "IVC anterior \n";
    }

    return res;
}

std::string AtrialFibresLandmarksView::GetRefinedPointsGiude(){
    std::string res = "REFINED LANDMARKS GUIDE\n Select specific locations for: \n";
    if(isLeftAtrium){
        res += "Lateral wall (LAA) - Between LSPV and MV, away from LAA\n";
        res += "Septal wall (FO)\n";
        res += "Posterior segment of LSPV/LA junction\n";
        res += "Posterior segment of RSPV/LA junction\n";
    } else{
        res += "IVC anterior \n";
		res += "CS\n";
        res += "IVC/SVC anterior \n";
        res += "SVC anterior \n";
        res += "RAA anterior \n";
        res += "RAA/CS anterior \n";
    }

    return res;
}


void AtrialFibresLandmarksView::UserSelectPvLabel(bool refinedLandmarks){
    if(!refinedLandmarks){
        UserSelectPvRoughLabel();
    } else{
        UserSelectPvRefinedLabel();
    }
}

void AtrialFibresLandmarksView::UserSelectPvRoughLabel(){
    int dialogCode = inputsRough->exec();
    QRect screenGeometry = QApplication::desktop()->screenGeometry();
    int x = (screenGeometry.width() - inputsRough->width()) / 2;
    int y = (screenGeometry.height() - inputsRough->height()) / 2;
    inputsRough->move(x,y);

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {

        if(isLeftAtrium){
            if (m_Rough.radioBtn_LA_LSPV->isChecked()) {
                roughSeedLabels.push_back(11); // LSPV
                m_Rough.radioBtn_LA_LSPV->setEnabled(false);
            } else if (m_Rough.radioBtn_LA_LIPV->isChecked()) {
                roughSeedLabels.push_back(13); // LIPV
                m_Rough.radioBtn_LA_LIPV->setEnabled(false);
            } else if (m_Rough.radioBtn_LA_RSPV->isChecked()) {
                roughSeedLabels.push_back(15); // RSPV
                m_Rough.radioBtn_LA_RSPV->setEnabled(false);
            } else if (m_Rough.radioBtn_LA_RIPV->isChecked()) {
                roughSeedLabels.push_back(17); // RIPV
                m_Rough.radioBtn_LA_RIPV->setEnabled(false);
            } else if (m_Rough.radioBtn_LAA_base->isChecked()) {
                roughSeedLabels.push_back(19); // LAAP_1
                m_Rough.radioBtn_LAA_base->setEnabled(false);
            } else if(m_Rough.radioBtn_LAA_tip){
                roughSeedLabels.push_back(21);
                m_Rough.radioBtn_LAA_tip->setEnabled(false);
            }
        } else{
            if(m_Rough.radioBtn_RA_SVC_POST->isChecked()){
                roughSeedLabels.push_back(29);
                m_Rough.radioBtn_RA_SVC_POST->setEnabled(false);
            } else if(m_Rough.radioBtn_RA_IVC_POST->isChecked()){
                roughSeedLabels.push_back(31);
                m_Rough.radioBtn_RA_IVC_POST->setEnabled(false);
            } else if(m_Rough.radioBtn_RAA_TCV->isChecked()){
                roughSeedLabels.push_back(33);
                m_Rough.radioBtn_RAA_TCV->setEnabled(false);
            } else if(m_Rough.radioBtn_RA_CS_TCV->isChecked()){
                roughSeedLabels.push_back(35);
                m_Rough.radioBtn_RA_CS_TCV->setEnabled(false);
            } else if(m_Rough.radioBtn_RA_SVC_ANT->isChecked()){
                roughSeedLabels.push_back(37);
                m_Rough.radioBtn_RA_SVC_ANT->setEnabled(false);
            } else if(m_Rough.radioBtn_RA_IVC_ANT->isChecked()){
                roughSeedLabels.push_back(39);
                m_Rough.radioBtn_RA_IVC_ANT->setEnabled(false);
            }
        }

    } else if (dialogCode == QDialog::Rejected) {
        inputsRough->close();
    }//_if
}

void AtrialFibresLandmarksView::UserSelectPvRefinedLabel(){
    int dialogCode = inputsRefined->exec();
    QRect screenGeometry = QApplication::desktop()->screenGeometry();
    int x = (screenGeometry.width() - inputsRefined->width()) / 2;
    int y = (screenGeometry.height() - inputsRefined->height()) / 2;
    inputsRefined->move(x,y);

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {
        if(isLeftAtrium){
            if (m_Refined.radioBtn_LA_LSPV->isChecked()) {
                refinedSeedLabels.push_back(11); // LSPV
                m_Refined.radioBtn_LA_LSPV->setEnabled(false);
            } else if (m_Refined.radioBtn_LA_LspvBody->isChecked()) {
                refinedSeedLabels.push_back(13); // LspvBody
                m_Refined.radioBtn_LA_LspvBody->setEnabled(false);
            } else if (m_Refined.radioBtn_LA_RSPV->isChecked()) {
                refinedSeedLabels.push_back(15); // RSPV
                m_Refined.radioBtn_LA_RSPV->setEnabled(false);
            } else if (m_Refined.radioBtn_LA_RspvBody->isChecked()) {
                refinedSeedLabels.push_back(17); // RspvBody
                m_Refined.radioBtn_LA_RspvBody->setEnabled(false);
            } else if (m_Refined.radioBtn_LA_LatWall->isChecked()) {
                refinedSeedLabels.push_back(19); // LatWall
                m_Refined.radioBtn_LA_LatWall->setEnabled(false);
            } else if(m_Refined.radioBtn_LA_FO){ // LAAP_1
                refinedSeedLabels.push_back(22);
                m_Refined.radioBtn_LA_FO->setEnabled(false);
            }
        } else {
            if(m_Refined.radioBtn_RA_IVC_ANT->isChecked()){
                refinedSeedLabels.push_back(29);
                m_Refined.radioBtn_RA_IVC_ANT->setEnabled(false);
            } else if(m_Refined.radioBtn_RA_CS->isChecked()){
                refinedSeedLabels.push_back(31);
                m_Refined.radioBtn_RA_CS->setEnabled(false);
            } else if(m_Refined.radioBtn_RA_IvcSvc->isChecked()){
                refinedSeedLabels.push_back(33);
                m_Refined.radioBtn_RA_IvcSvc->setEnabled(false);
            } else if(m_Refined.radioBtn_RA_SVC_ANT->isChecked()){
                refinedSeedLabels.push_back(35);
                m_Refined.radioBtn_RA_SVC_ANT->setEnabled(false);
            } else if(m_Refined.radioBtn_RAA_ANT->isChecked()){
                refinedSeedLabels.push_back(37);
                m_Refined.radioBtn_RAA_ANT->setEnabled(false);
            } else if(m_Refined.radioBtn_RAA_CS_ANT->isChecked()){
                refinedSeedLabels.push_back(39);
                m_Refined.radioBtn_RAA_CS_ANT->setEnabled(false);
            }
        }

    } else if (dialogCode == QDialog::Rejected) {
        inputsRefined->close();
    }//_if
}

std::string AtrialFibresLandmarksView::GetStructureIdFromLabel(bool refinedLandmarks, int label){
    QString res;
    if(!refinedLandmarks){
        if(label==11){
            res = "LSPV";
        }else if(label==13){
            res = "LIPV";
        }else if(label==15){
            res = "RSPV";
        }else if(label==17){
            res = "RIPV";
        }else if(label==19){
            res = "LAA_BASE";
        }else if(label==21){
            res = "LAA_TIP";
        } else if(label==29){
            res = "SVC_POST";
        } else if(label==31){
            res = "IVC_POST";
        } else if(label==33){
            res = "RAA_VALVE_P";
        } else if(label==35){
            res = "CS_VALVE_P";
        } else if(label==37){
            res = "SVC_ANT";
        } else if(label==39){
            res = "IVC_ANT";
        }
    } else{
        if(label==11){
            res = "LSPV_ROOF";
        }else if(label==13){
            res = "LSPV_POST";
        }else if(label==15){
            res = "RSPV_ROOF";
        }else if(label==17){
            res = "RIPV_POST";
        }else if(label==19){
            res = "LAA";
        }else if(label==22){
            res = "FO";
        } else if(label==29){
            res = "IVC_ANT";
        } else if(label==31){
            res = "CS_TOP";
        } else if(label==33){
            res = "IVC_SVC_ANT";
        } else if(label==35){
            res = "SVC_ANT";
        } else if(label==37){
            res = "RAA_ANT";
        } else if(label==39){
            res = "RAA_CS_ANT";
        }
    }

    return res.toStdString();
}

int AtrialFibresLandmarksView::GetIndex(std::vector<int> v, int value){
    int index=-1;
    auto it = std::find(v.begin(), v.end(), value);
    if(it != v.end()){
        index = it - v.begin();
    }
    return index;
}

void AtrialFibresLandmarksView::RoughUiEnableButtons(){
    m_Rough.radioBtn_LAA_base->setVisible(isLeftAtrium);
    m_Rough.radioBtn_LAA_tip->setVisible(isLeftAtrium);
    m_Rough.radioBtn_LA_LSPV->setVisible(isLeftAtrium);
    m_Rough.radioBtn_LA_LIPV->setVisible(isLeftAtrium);
    m_Rough.radioBtn_LA_RSPV->setVisible(isLeftAtrium);
    m_Rough.radioBtn_LA_RIPV->setVisible(isLeftAtrium);

    m_Rough.radioBtn_RA_SVC_POST->setVisible(!isLeftAtrium);
    m_Rough.radioBtn_RA_IVC_POST->setVisible(!isLeftAtrium);
    m_Rough.radioBtn_RAA_TCV->setVisible(!isLeftAtrium);
    m_Rough.radioBtn_RA_CS_TCV->setVisible(!isLeftAtrium);
    m_Rough.radioBtn_RA_SVC_ANT->setVisible(!isLeftAtrium);
    m_Rough.radioBtn_RA_IVC_ANT->setVisible(!isLeftAtrium);
}

void AtrialFibresLandmarksView::RefinedUiEnableButtons(){

    m_Refined.radioBtn_LA_FO->setVisible(isLeftAtrium);
    m_Refined.radioBtn_LA_LSPV->setVisible(isLeftAtrium);
    m_Refined.radioBtn_LA_LatWall->setVisible(isLeftAtrium);
    m_Refined.radioBtn_LA_LspvBody->setVisible(isLeftAtrium);
    m_Refined.radioBtn_LA_RSPV->setVisible(isLeftAtrium);
    m_Refined.radioBtn_LA_RspvBody->setVisible(isLeftAtrium);

    m_Refined.radioBtn_RA_IVC_ANT->setVisible(!isLeftAtrium);
    m_Refined.radioBtn_RA_CS->setVisible(!isLeftAtrium);
    m_Refined.radioBtn_RA_IvcSvc->setVisible(!isLeftAtrium);
    m_Refined.radioBtn_RA_SVC_ANT->setVisible(!isLeftAtrium);
    m_Refined.radioBtn_RAA_ANT->setVisible(!isLeftAtrium);
    m_Refined.radioBtn_RAA_CS_ANT->setVisible(!isLeftAtrium);
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
