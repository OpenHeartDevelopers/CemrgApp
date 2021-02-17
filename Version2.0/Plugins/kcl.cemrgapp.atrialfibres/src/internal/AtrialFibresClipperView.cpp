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
#include "AtrialFibresClipperView.h"
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

QString AtrialFibresClipperView::fileName;
QString AtrialFibresClipperView::directory;
bool AtrialFibresClipperView::isAutomatic;
const std::string AtrialFibresClipperView::VIEW_ID = "org.mitk.views.atrialfibresclipperview";

void AtrialFibresClipperView::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_man1_ctrlines, SIGNAL(clicked()), this, SLOT(CtrLines()));
    connect(m_Controls.button_man2_clippers, SIGNAL(clicked()), this, SLOT(CtrPlanes()));
    connect(m_Controls.button_man3_clipseg, SIGNAL(clicked()), this, SLOT(ClipperImage()));
    connect(m_Controls.button_auto1_savelabels, SIGNAL(clicked()), this, SLOT(SaveLabels()));
    connect(m_Controls.button_auto2_clippers, SIGNAL(clicked()), this, SLOT(ClipperBalls()));
    connect(m_Controls.button_auto3_spacing, SIGNAL(clicked()), this, SLOT(InterPvSpacing()));
    connect(m_Controls.slider, SIGNAL(valueChanged(int)), this, SLOT(CtrPlanesPlacer()));
    connect(m_Controls.spinBox, SIGNAL(valueChanged(double)), this, SLOT(CtrPlanesPlacer()));
    connect(m_Controls.comboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(CtrLinesSelector(int)));

    // Display correct buttons
    automaticPipeline = AtrialFibresClipperView::isAutomatic;
    SetAutomaticModeButtons(automaticPipeline);
    SetManualModeButtons(!automaticPipeline);
    corridorMax = 5;
    corridorCount = 0;

    //Create GUI widgets
    inputs = new QDialog(0,0);
    m_Labels.setupUi(inputs);
    connect(m_Labels.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_Labels.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));

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
        clipper = std::unique_ptr<CemrgAtriaClipper>(new CemrgAtriaClipper(directory, surface));
        Visualiser();
    }
    m_Controls.slider->setEnabled(false);
    m_Controls.spinBox->setEnabled(false);
    m_Controls.comboBox->setEnabled(false);
}

void AtrialFibresClipperView::SetFocus() {
    m_Controls.button_man1_ctrlines->setFocus();
}

void AtrialFibresClipperView::OnSelectionChanged(
        berry::IWorkbenchPart::Pointer /*src*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

AtrialFibresClipperView::~AtrialFibresClipperView() {
    inputs->deleteLater();
}

void AtrialFibresClipperView::SetDirectoryFile(const QString directory, const QString fileName, const bool isAutomatic) {
    AtrialFibresClipperView::fileName = fileName;
    AtrialFibresClipperView::directory = directory;
    AtrialFibresClipperView::isAutomatic = isAutomatic;
}

void AtrialFibresClipperView::iniPreSurf() {
    //Find the selected node
    QString path = AtrialFibresClipperView::directory + mitk::IOUtil::GetDirectorySeparator() + AtrialFibresClipperView::fileName;
    mitk::Surface::Pointer shell = mitk::IOUtil::Load<mitk::Surface>(path.toStdString());
    // mitk::Surface::Pointer shell = CemrgCommonUtils::LoadVTKMesh(path.toStdString());
    // CemrgCommonUtils::FlipXYPlane(shell, "", "");

    vtkSmartPointer<vtkCellDataToPointData> cell_to_point = vtkSmartPointer<vtkCellDataToPointData>::New();
    cell_to_point->SetInputData(shell->GetVtkPolyData());
    cell_to_point->PassCellDataOn();
    cell_to_point->Update();
    shell->SetVtkPolyData(cell_to_point->GetPolyDataOutput());

    vtkSmartPointer<vtkPolyData> pd = cell_to_point->GetPolyDataOutput();
    vtkFloatArray *pointScalars = vtkFloatArray::New();
    vtkFloatArray *cellScalars = vtkFloatArray::New();

    pointScalars = vtkFloatArray::SafeDownCast(pd->GetPointData()->GetScalars());
    cellScalars = vtkFloatArray::SafeDownCast(pd->GetCellData()->GetScalars());
    double s;
    for (vtkIdType vId = 0; vId < pd->GetNumberOfPoints() ; vId++) {
        s = pointScalars->GetTuple1(vId);

        if (std::floor(s)!=s){ // scalar value is not a category
            int s2;
            vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
            cellIds->Initialize();

            pd->GetPointCells(vId, cellIds);
            s2 = cellScalars->GetTuple1(cellIds->GetId(0));
            pointScalars->SetTuple1(vId, s2);
        }
    }
    shell->GetVtkPolyData()->GetPointData()->SetScalars(pointScalars);
    mitk::IOUtil::Save(shell, path.toStdString());
    surface = shell;

}

// Manual pipeline
void AtrialFibresClipperView::CtrLines() {

    if (clipper->GetCentreLines().size() == 0 && pickedSeedLabels.size() == 0) {

        QMessageBox::warning(NULL, "Attention", "Please maske sure you have selected seeds!");
        return;

    } else {

        //Retrieve centrelines and labels
        if (clipper->GetCentreLines().size() != 0)
            QMessageBox::information(NULL, "Attention", "You are using precomputed centrelines!");

        this->BusyCursorOn();
        for (unsigned int i=0; i<clipper->GetCentreLines().size(); i++) {
            if (i == 0) pickedSeedLabels.clear();
            vtkSmartPointer<vtkPolyData> pd = clipper->GetCentreLines().at(i)->GetOutput();
            vtkSmartPointer<vtkIntArray> lb = vtkIntArray::SafeDownCast(pd->GetFieldData()->GetAbstractArray("PickedSeedLabels"));
            pickedSeedLabels.push_back(lb->GetValue(0));
        }//_for
        bool successful = clipper->ComputeCtrLines(pickedSeedLabels, pickedSeedIds, m_Controls.checkBox->isChecked());
        this->BusyCursorOff();

        //Check for failure
        if (!successful) {
            QMessageBox::critical(NULL, "Attention", "Computation of Centrelines Failed!");
            return;
        }//_if
    }//_if

    //Set surface opacity
    renderer->RemoveAllViewProps();
    Visualiser(0.5);

    //Setup shortcut keys
    std::string sk;
    vtkSmartPointer<vtkTextActor> txtActor = vtkSmartPointer<vtkTextActor>::New();
    sk = "A: revert to automatic\nO: change opacity\nArrows: move planes\nSpace: add clip point\nDelete: remove clip point";
    txtActor->SetInput(sk.c_str());
    txtActor->GetTextProperty()->SetFontSize(14);
    txtActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
    renderer->AddActor2D(txtActor);

    //Create a mapper and actor for centre lines
    std::vector<vtkSmartPointer<vtkvmtkPolyDataCenterlines>> ctrLines = clipper->GetCentreLines();
    for (unsigned int i=0; i<ctrLines.size(); i++) {
        vtkSmartPointer<vtkPolyDataMapper> linesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        linesMapper->SetInputData(ctrLines.at(i)->GetOutput());
        linesMapper->ScalarVisibilityOff();
        vtkSmartPointer<vtkActor> linesActor = vtkSmartPointer<vtkActor>::New();
        linesActor->SetMapper(linesMapper);
        linesActor->GetProperty()->SetOpacity(1.0);
        linesActor->GetProperty()->SetColor(1,0,0);
        renderer->AddActor(linesActor);
    }//_for
    m_Controls.widget_1->GetRenderWindow()->Render();

    //Adjust controllers
    m_Controls.button_man1_ctrlines->setEnabled(false);
    m_Controls.checkBox->setEnabled(false);
}

void AtrialFibresClipperView::CtrPlanes() {

    if (clipper->GetCentreLines().size() == 0 || pickedSeedLabels.size() == 0) {

        QMessageBox::warning(NULL, "Attention", "Please maske sure you have computed centre lines!");
        return;

    } else {

        this->BusyCursorOn();
        bool successful = clipper->ComputeCtrLinesClippers(pickedSeedLabels);
        this->BusyCursorOff();

        //Check for failure
        if (!successful) {
            QMessageBox::critical(NULL, "Attention", "Computation of Clipper Planes Failed!");
            return;
        }//_if
    }//_if

    std::vector<vtkSmartPointer<vtkRegularPolygonSource>> ctrPlanes = clipper->GetCentreLinePolyPlanes();
    std::vector<vtkSmartPointer<vtkvmtkPolyDataCenterlines>> ctrLines = clipper->GetCentreLines();
    vtkSmartPointer<vtkRegularPolygonSource> ctrPlane = ctrPlanes.at(0);
    vtkSmartPointer<vtkvmtkPolyDataCenterlines> ctrLine = ctrLines.at(0);
    int index = ctrLine->GetOutput()->FindPoint(ctrPlane->GetCenter());

    //Create a mapper and actor for centre lines clippers
    for (unsigned int i=0; i<ctrPlanes.size(); i++) {
        vtkSmartPointer<vtkPolyDataMapper> clipperMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        clipperMapper->SetInputConnection(ctrPlanes.at(i)->GetOutputPort());
        clipperMapper->ScalarVisibilityOff();
        vtkSmartPointer<vtkActor> clipperActor = vtkSmartPointer<vtkActor>::New();
        clipperActor->SetMapper(clipperMapper);
        clipperActor->GetProperty()->SetOpacity(1.0);
        clipperActor->GetProperty()->SetColor(1,1,0);
        renderer->AddActor(clipperActor);
        clipperActors.push_back(clipperActor);
        QString comboText = "DEFAULT";
        if (pickedSeedLabels.at(i) == 11)
            comboText = "LEFT SUPERIOR PV";
        else if (pickedSeedLabels.at(i) == 12)
            comboText = "LEFT MIDDLE PV";
        else if (pickedSeedLabels.at(i) == 13)
            comboText = "LEFT INFERIOR PV";
        else if (pickedSeedLabels.at(i) == 14)
            comboText = "LEFT COMMON PV";
        else if (pickedSeedLabels.at(i) == 15)
            comboText = "RIGHT SUPERIOR PV";
        else if (pickedSeedLabels.at(i) == 16)
            comboText = "RIGHT MIDDLE PV";
        else if (pickedSeedLabels.at(i) == 17)
            comboText = "RIGHT INFERIOR PV";
        else if (pickedSeedLabels.at(i) == 18)
            comboText = "RIGHT COMMON PV";
        else if (pickedSeedLabels.at(i) == 19)
            comboText = "APPENDAGE CUT";
        else if (pickedSeedLabels.at(i) == 20)
            comboText = "APPENDAGE UNCUT";
        m_Controls.comboBox->insertItem(i, comboText);
    }//_for
    m_Controls.widget_1->GetRenderWindow()->Render();

    //Adjust controllers
    m_Controls.comboBox->setCurrentIndex(0);
    m_Controls.slider->setValue(index);
    m_Controls.slider->setEnabled(true);
    m_Controls.spinBox->setEnabled(true);
    m_Controls.comboBox->setEnabled(true);
    m_Controls.button_man2_clippers->setEnabled(false);
}

void AtrialFibresClipperView::ClipperImage() {

    if (this->GetDataManagerSelection().empty()) {
        QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation to clip!");
        return;
    } else if (clipper->GetCentreLines().size() == 0) {
        QMessageBox::warning(NULL, "Attention", "Please maske sure you have computed centre lines!");
        return;
    } else if (clipper->GetCentreLinePolyPlanes().size() == 0) {
        QMessageBox::warning(NULL, "Attention", "Please maske sure you have computed clipper planes!");
        return;
    } else {

        //Check for selection of segmentation image
        QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
        mitk::DataNode::Pointer segNode = nodes.at(0);
        mitk::BaseData::Pointer data = segNode->GetData();
        if (data) {

            //Test if this data item is an image
            mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
            if (image) {

                //Check if the right segmentation
                bool morphAnalysis = false;
                if (segNode->GetName().compare(fileName.left(fileName.lastIndexOf(QChar('.'))).toStdString()) == 0) {
                    int reply1 = QMessageBox::question(
                                NULL, "Question", "Are you sure you want to clip the same segmentation used for creating the visualised mesh?",
                                QMessageBox::Yes, QMessageBox::No);
                    if (reply1 == QMessageBox::No) {
                        QMessageBox::warning(NULL, "Attention", "Please select the correct segmentation to clip!");
                        return;
                    }//_if_no
                    int reply2 = QMessageBox::question(
                                NULL, "Question", "Do you want to save PV clipping information for morphological analysis?",
                                QMessageBox::Yes, QMessageBox::No);
                    if (reply2 == QMessageBox::Yes) {
                        QMessageBox::information(NULL, "Attention", "Labelled bloodpool image and ostia information will be saved!");
                        morphAnalysis = true;
                    }//_if_yes
                }//_if

                this->BusyCursorOn();
                this->GetDataStorage()->Remove(segNode);
                clipper->ClipVeinsImage(pickedSeedLabels, image, morphAnalysis ? true : false);
                this->BusyCursorOff();
                QMessageBox::information(NULL, "Attention", "Segmentation is now clipped!");

            } else {
                QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation!");
                return;
            }//_image
        } else {
            QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation!");
            return;
        }//_data
    }//_if

    //Add to storage
    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(clipper->GetClippedSegImage());
    node->SetName("PVeinsCroppedImage");
    this->GetDataStorage()->Add(node);

    //Adjust controllers
    m_Controls.slider->setEnabled(false);
    m_Controls.spinBox->setEnabled(false);
    m_Controls.comboBox->setEnabled(false);
    m_Controls.button_man3_clipseg->setEnabled(false);
}

void AtrialFibresClipperView::CtrPlanesPlacer() {

    if (m_Controls.comboBox->count() == 0)
        return;

    //Retrieve cutting planes
    std::vector<vtkSmartPointer<vtkRegularPolygonSource>> ctrPlanes = clipper->GetCentreLinePolyPlanes();
    vtkSmartPointer<vtkRegularPolygonSource> ctrPlane = ctrPlanes.at(m_Controls.comboBox->currentIndex());

    //Parameters recalculations to adjust
    int position = m_Controls.slider->value();
    double value = m_Controls.spinBox->value();
    int indexBox = m_Controls.comboBox->currentIndex();

    clipper->SetRadiusAdjustment(value);
    clipper->CalcParamsOfPlane(ctrPlane, indexBox, position);

    ctrPlane->Update();
    m_Controls.widget_1->GetRenderWindow()->Render();
}

void AtrialFibresClipperView::CtrLinesSelector(int index) {

    if (m_Controls.comboBox->count() == 0)
        return;

    //Find the right position
    vtkSmartPointer<vtkPolyData> line = clipper->GetCentreLines().at(index)->GetOutput();
    vtkSmartPointer<vtkRegularPolygonSource> ctrPlane = clipper->GetCentreLinePolyPlanes().at(index);
    int position = line->FindPoint(ctrPlane->GetCenter());

    //Find the right radius
    vtkSmartPointer<vtkDoubleArray> radii = vtkSmartPointer<vtkDoubleArray>::New();
    radii = vtkDoubleArray::SafeDownCast(line->GetPointData()->GetArray("MaximumInscribedSphereRadius"));
    double adjust = ctrPlane->GetRadius() / radii->GetValue(position);

    //Adjust highlighted clipper
    for (unsigned int i=0; i<clipperActors.size(); i++)
        clipperActors.at(i)->GetProperty()->SetOpacity(i==(unsigned)index ? 1.0 : 0.1);

    //Set corresponding label
    int type = clipper->GetManualType().at(index);
    if (type == 0) {
        m_Controls.label->setText(" Automatic ");
        m_Controls.label->setStyleSheet("QLabel {border-width:1px; border-color:black; border-radius:10px; background-color:yellow;}");
    } else if (type == 1) {
        m_Controls.label->setText(" Manual--1 ");
        m_Controls.label->setStyleSheet("QLabel {border-width:1px; border-color:black; border-radius:10px; background-color:lime;}");
    } else {
        m_Controls.label->setText(" Manual--2 ");
        m_Controls.label->setStyleSheet("QLabel {border-width:1px; border-color:black; border-radius:10px; background-color:red;}");
    }//_if_type

    //Set controllers
    m_Controls.slider->setRange(0, line->GetNumberOfPoints()-1);
    m_Controls.slider->setValue(position);
    m_Controls.spinBox->setValue(adjust);
    pickedCutterSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
    m_Controls.widget_1->GetRenderWindow()->Render();
}

// Automatic pipeline
void AtrialFibresClipperView::SaveLabels(){
    if (pickedSeedIds->GetNumberOfIds()==0) {
        std::string msg = "Select the corresponding veins and assign labels to them.";
        QMessageBox::warning(NULL, "Attention - No points selected.", msg.c_str());
        MITK_WARN << msg;

        return;
    }

    MITK_INFO << "[SaveLabels] Saving labels to file.";
    QString prodPath = directory + mitk::IOUtil::GetDirectorySeparator();
    ofstream fileLabels, fileIds, fileIgnoreIds, fileDiscardIds;

    fileLabels.open((prodPath + "prodSeedLabels.txt").toStdString());
    fileIds.open((prodPath + "prodSeedIds.txt").toStdString());
    fileIgnoreIds.open((prodPath + "prodIgnoreSeedIds.txt").toStdString());
    fileDiscardIds.open((prodPath + "prodDiscardSeedIds.txt").toStdString());

    // 14=ignore, 18=discard
    for (unsigned int i=0; i<pickedSeedLabels.size(); i++){
        if(pickedSeedLabels.at(i)==14){ //ignore
            fileIgnoreIds << pickedSeedIds->GetId(i) << "\n";
        } else if(pickedSeedLabels.at(i)==18){ // discard
            fileDiscardIds << pickedSeedIds->GetId(i) << "\n";
        } else{
            fileLabels << pickedSeedLabels.at(i) << "\n";
            fileIds << pickedSeedIds->GetId(i) << "\n";
        }
    }
    fileLabels.close();
    fileIds.close();

}

void AtrialFibresClipperView::ClipperBalls(){
    if(pickedSeedLabels.size()==0){
        MITK_WARN << "[ClipperBalls] showing clippers error";
        return;
    }
    // output file:
    QString outp = AtrialFibresClipperView::directory + mitk::IOUtil::GetDirectorySeparator() +
    outp += "clipper_radius.txt";

    std::ofstream fo(outp.toStdString());

    for (unsigned int i=0; i<pickedSeedLabels.size(); i++){
        if(pickedSeedLabels.at(i)!=14 && pickedSeedLabels.at(i)!=18 && pickedSeedLabels.at(i)!=19){
            double* point = surface->GetVtkPolyData()->GetPoint(pickedSeedIds->GetId(i));
            double radius = 9;
            bool radiusAccepted=false;

            while(!radiusAccepted){
                vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
                sphereSource->SetCenter(point[0], point[1], point[2]);
                sphereSource->SetRadius(radius);
                sphereSource->SetPhiResolution(40);
                sphereSource->SetThetaResolution(40);
                sphereSource->Update();

                VisualiseSphere(sphereSource->GetOutput());
                m_Controls.widget_1->GetRenderWindow()->Render();

                QDialog* inputs = new QDialog(0,0);
                m_UIRadius.setupUi(inputs);
                connect(m_UIRadius.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
                connect(m_UIRadius.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
                int dialogCode = inputs->exec();

                //Act on dialog return code
                if (dialogCode == QDialog::Accepted) {
                    bool ok1;
                    radius = m_UIRadius.lineEdit_radius->text().toDouble(&ok1);
                    if(!ok1) radius = 9;

                    int reply = QMessageBox::question(NULL, "Question", "Is this radius OK?", QMessageBox::Yes, QMessageBox::No);
                    radiusAccepted = (reply == QMessageBox::Yes);
                }
            }
            fo << pickedSeedIds->GetId(i) << ", " << radius << std::endl;
        }
    }


}

void AtrialFibresClipperView::InterPvSpacing(){
    if (corridorSeedIds->GetNumberOfIds()==0) {
        std::string msg = "Select the corresponding veins and assign labels to them.";
        QMessageBox::warning(NULL, "Attention - No points selected.", msg.c_str());
        MITK_WARN << msg;

        return;
    }
    MITK_INFO << "[InterPvSpacing] Calculating inter PV spacing";
    int numCorridorIds = this->corridorSeedIds->GetNumberOfIds();
    std::vector<int> idVectors;
    for(int i=0; i<numCorridorIds; i++){
        idVectors.push_back(this->corridorSeedIds->GetId(i));
    }

    MITK_INFO << "[InterPvSpacing] Creating shortest path and corridor";
    QString prodPathOut = AtrialFibresClipperView::directory + mitk::IOUtil::GetDirectorySeparator();
    // csadv parameters
    bool circleCorridor = false;
    int thickness = 3;
    std::string lrpre = "";

    csadv = std::unique_ptr<CemrgScarAdvanced>(new CemrgScarAdvanced());
    csadv->SetOutputFileName((prodPathOut+"corridor.csv").toStdString());
    csadv->SetOutputPath(prodPathOut.toStdString());
    csadv->SetInputData(surface->GetVtkPolyData());
    csadv->SetWeightedCorridorBool(false);
    csadv->SetLeftRightPrefix(lrpre);
    csadv->CorridorFromPointList(idVectors, circleCorridor);

    MITK_INFO << "[InterPvSpacing] Setting values in corridor to atrium body label";

    QString path2corridor = prodPathOut + "corridor.csv";
    std::ifstream fi(path2corridor.toStdString());
    MITK_INFO << ("[InterPvSpacing] Opened file :" + path2corridor).toStdString();

    vtkFloatArray *pointScalars = vtkFloatArray::New();
    vtkFloatArray *cellScalars = vtkFloatArray::New();
    pointScalars = vtkFloatArray::SafeDownCast(surface->GetVtkPolyData()->GetPointData()->GetScalars());
    cellScalars = vtkFloatArray::SafeDownCast(surface->GetVtkPolyData()->GetCellData()->GetScalars());

    std::string line, header;
    std::getline(fi, header);

    while(std::getline(fi, line)){
        QString qline =  QString::fromStdString(line);
        vtkIdType vId;
        double vScalar;
        vId = qline.section(',', 1, 1).toInt();
        vScalar = qline.section(',', 6, 6).toDouble();

        if (vScalar != 1){
            pointScalars->SetTuple1(vId, 1);

            vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
            cellIds->Initialize();

            surface->GetVtkPolyData()->GetPointCells(vId, cellIds);
            for (vtkIdType ix = 0; ix < cellIds->GetNumberOfIds() ; ix++) {
                cellScalars->SetTuple1(cellIds->GetId(ix), 1);
            }
        }
    }
    fi.close();

    surface->GetVtkPolyData()->GetPointData()->SetScalars(pointScalars);
    surface->GetVtkPolyData()->GetCellData()->SetScalars(cellScalars);

    MITK_INFO << "[InterPvSpacing] Set new scalar vector into surface.";
    mitk::IOUtil::Save(surface, (prodPathOut+AtrialFibresClipperView::fileName).toStdString());
    MITK_INFO << "[InterPvSpacing] Saved surface";

    Visualiser();

}

void AtrialFibresClipperView::Visualiser(double opacity) {

    MITK_INFO << "Visualiser";
    double max_scalar=-1, min_scalar=1e9,s;
    vtkFloatArray *scalars = vtkFloatArray::New();
    vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
    scalars = vtkFloatArray::SafeDownCast(surface->GetVtkPolyData()->GetPointData()->GetScalars());
    for (vtkIdType i=0;i<surface->GetVtkPolyData()->GetNumberOfPoints();i++) {
        s = scalars->GetTuple1(i);
        if (s > max_scalar)
            max_scalar = s;
        if (s < min_scalar)
            min_scalar = s;
    }
    this->maxScalar = max_scalar;
    this->minScalar = min_scalar;

    SphereSourceVisualiser(pickedLineSeeds);
    SphereSourceVisualiser(corridorLineSeeds, "0.0,0.0,1.0");

    //Create a mapper and actor for surface
    vtkSmartPointer<vtkPolyDataMapper> surfMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    surfMapper->SetInputData(surface->GetVtkPolyData());
    surfMapper->SetScalarRange(min_scalar, max_scalar);
    surfMapper->SetScalarModeToUseCellData();
    surfMapper->ScalarVisibilityOn();
    lut->SetTableRange(min_scalar, max_scalar);
    lut->SetHueRange(0.3, 0.0);  // this is the way_neighbourhood_size you tell which colors you want to be displayed.
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

void AtrialFibresClipperView::VisualiseSphere(vtkSmartPointer<vtkPolyData> pd){
    vtkSmartPointer<vtkPolyDataMapper> surfMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor> surfActor = vtkSmartPointer<vtkActor>::New();

    surfMapper->SetInputData(pd);
    surfActor->SetMapper(surfMapper);
    surfActor->GetProperty()->SetOpacity(0.5);
    renderer->AddActor(surfActor);
}

void AtrialFibresClipperView::SphereSourceVisualiser(vtkSmartPointer<vtkPolyData> pointSources, QString colour, double scaleFactor){
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

void AtrialFibresClipperView::PickCallBack(bool pvCorridor) {

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
    if(!pvCorridor){
        pickedSeedIds->InsertNextId(pickedSeedId);
        pickedLineSeeds->GetPoints()->InsertNextPoint(point);
        pickedLineSeeds->Modified();
    } else{
        corridorSeedIds->InsertNextId(pickedSeedId);
        corridorLineSeeds->GetPoints()->InsertNextPoint(point);
        corridorLineSeeds->Modified();
    }

    m_Controls.widget_1->GetRenderWindow()->Render();
}

void AtrialFibresClipperView::ManualCutterCallBack() {

    if (m_Controls.label->text() != " Manual--2 ") {

        SphereSourceVisualiser(pickedCutterSeeds);

        //Adjust labels
        m_Controls.label->setText(" Manual--2 ");
        m_Controls.label->setStyleSheet("QLabel {border-width:1px; border-color:black; border-radius:10px; background-color:red;}");
    }//_if

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
    if (pickedSeedId == -1)
        pickedSeedId = pickedCellPointIds->GetId(0);

    double* point = surface->GetVtkPolyData()->GetPoint(pickedSeedId);
    pickedCutterSeeds->GetPoints()->InsertNextPoint(point);
    pickedCutterSeeds->Modified();
    m_Controls.widget_1->GetRenderWindow()->Render();
}

void AtrialFibresClipperView::KeyCallBackFunc(vtkObject*, long unsigned int, void* ClientData, void*) {

    AtrialFibresClipperView* self;
    self = reinterpret_cast<AtrialFibresClipperView*>(ClientData);
    std::string key = self->interactor->GetKeySym();

    if (self->IsPointSelectionControlsAvailable()) {

        if (key == "space") {
            //Ask the labels
            self->PickCallBack();
            self->UserSelectPvLabel();

        } else if (key == "Delete") {

            //Clean up last dropped seed point
            vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
            vtkSmartPointer<vtkPoints> points = self->pickedLineSeeds->GetPoints();
            for (int i=0; i<points->GetNumberOfPoints()-1; i++){
                newPoints->InsertNextPoint(points->GetPoint(i));
            }
            self->pickedLineSeeds->SetPoints(newPoints);
            vtkSmartPointer<vtkIdList> newPickedSeedIds = vtkSmartPointer<vtkIdList>::New();
            newPickedSeedIds->Initialize();
            vtkSmartPointer<vtkIdList> pickedSeedIds = self->pickedSeedIds;
            for (int i=0; i<pickedSeedIds->GetNumberOfIds()-1; i++){
                newPickedSeedIds->InsertNextId(pickedSeedIds->GetId(i));
            }
            self->pickedSeedIds = newPickedSeedIds;

            if (self->pickedSeedLabels.empty() == false) {
                int radioButtonNumber = self->pickedSeedLabels.back() - 10;
                if (radioButtonNumber == 1)
                self->m_Labels.radioBtn_LSPV->setEnabled(true);
                else if (radioButtonNumber == 3)
                self->m_Labels.radioBtn_LIPV->setEnabled(true);
                else if (radioButtonNumber == 4)
                self->m_Labels.radioBtn_ignore->setEnabled(true);
                else if (radioButtonNumber == 5)
                self->m_Labels.radioBtn_RSPV->setEnabled(true);
                else if (radioButtonNumber == 7)
                self->m_Labels.radioBtn_RIPV->setEnabled(true);
                else if (radioButtonNumber == 8)
                self->m_Labels.radioBtn_discard->setEnabled(true);
                else if (radioButtonNumber == 9)
                self->m_Labels.radioBtn_LAA->setEnabled(true);

                self->pickedSeedLabels.pop_back();
            }//_if

            self->m_Controls.widget_1->GetRenderWindow()->Render();
        } else if (key == "X" || key == "x"){
            if(self->automaticPipeline){
                if(self->corridorCount<self->corridorMax){
                    bool pvCorridor = true;
                    self->PickCallBack(pvCorridor);
                    self->corridorCount++;
                } else{
                    std::string msg = "Select less than 5 points for inter PV corridor!\nPress V to fix another vein.";
                    QMessageBox::warning(NULL, "Attention - too many corridor points", msg.c_str());
                    MITK_INFO << msg;
                }
            }
        } else if (key == "D" || key == "d"){
            if(self->automaticPipeline){
                //Clean up last dropped seed point
                vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
                vtkSmartPointer<vtkPoints> points = self->corridorLineSeeds->GetPoints();
                for (int i=0; i<points->GetNumberOfPoints()-1; i++){
                    newPoints->InsertNextPoint(points->GetPoint(i));
                }
                self->corridorLineSeeds->SetPoints(newPoints);
                vtkSmartPointer<vtkIdList> newPickedSeedIds = vtkSmartPointer<vtkIdList>::New();
                newPickedSeedIds->Initialize();
                vtkSmartPointer<vtkIdList> pickedSeedIds = self->corridorSeedIds;
                for (int i=0; i<pickedSeedIds->GetNumberOfIds()-1; i++){
                    newPickedSeedIds->InsertNextId(pickedSeedIds->GetId(i));
                }
                self->corridorSeedIds = newPickedSeedIds;
                self->corridorCount--;
            }
        } else if(key == "V"){ // new veins
                if(self->automaticPipeline){
                    int reply1 = QMessageBox::question(NULL, "Question - different PV",
                        "Do you want to fix another pair of veins?", QMessageBox::Yes, QMessageBox::No);
                    if(reply1==QMessageBox::Yes){
                    self->corridorCount=0;
                    self->corridorSeedIds = vtkSmartPointer<vtkIdList>::New();
                    self->corridorSeedIds->Initialize();
                    self->corridorLineSeeds = vtkSmartPointer<vtkPolyData>::New();
                    self->corridorLineSeeds->Initialize();
                    self->corridorLineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
                }
            }
        } //_if_space

    } else if (self->IsClipperManualControlsAvailable()) {

        double adjustments[2];
        int idxClipper = self->m_Controls.comboBox->currentIndex();
        adjustments[0] = self->clipper->GetMClipperAngles()[idxClipper][0];
        adjustments[1] = self->clipper->GetMClipperAngles()[idxClipper][1];

        if (key == "Up") {

            self->m_Controls.label->setText(" Manual--1 ");
            QString style = "QLabel {border-width:1px; border-color:black; border-radius:10px; background-color:lime;}";
            self->m_Controls.label->setStyleSheet(style);
            if (self->clipper->GetMClipperAngles()[idxClipper][0] >= 180.0)
            adjustments[0] = 0.0;
            else
            adjustments[0] = self->clipper->GetMClipperAngles()[idxClipper][0] + 0.1;
            self->clipperActors.at(idxClipper)->GetProperty()->SetColor(0,1,0);
            self->clipper->SetMClipperAngles(adjustments, idxClipper);
            self->CtrPlanesPlacer();

        } else if (key == "Down") {

            self->m_Controls.label->setText(" Manual--1 ");
            QString style = "QLabel {border-width:1px; border-color:black; border-radius:10px; background-color:lime;}";
            self->m_Controls.label->setStyleSheet(style);
            if (self->clipper->GetMClipperAngles()[idxClipper][0] <= 0.0)
            adjustments[0] = 180.0;
            else
            adjustments[0] = self->clipper->GetMClipperAngles()[idxClipper][0] - 0.1;
            self->clipperActors.at(idxClipper)->GetProperty()->SetColor(0,1,0);
            self->clipper->SetMClipperAngles(adjustments, idxClipper);
            self->CtrPlanesPlacer();

        } else if (key == "Right") {

            self->m_Controls.label->setText(" Manual--1 ");
            QString style = "QLabel {border-width:1px; border-color:black; border-radius:10px; background-color:lime;}";
            self->m_Controls.label->setStyleSheet(style);
            if (self->clipper->GetMClipperAngles()[idxClipper][1] >= 360.0)
            adjustments[1] = 0.0;
            else
            adjustments[1] = self->clipper->GetMClipperAngles()[idxClipper][1] + 0.1;
            self->clipperActors.at(idxClipper)->GetProperty()->SetColor(0,1,0);
            self->clipper->SetMClipperAngles(adjustments, idxClipper);
            self->CtrPlanesPlacer();

        } else if (key == "Left") {

            self->m_Controls.label->setText(" Manual--1 ");
            QString style = "QLabel {border-width:1px; border-color:black; border-radius:10px; background-color:lime;}";
            self->m_Controls.label->setStyleSheet(style);
            if (self->clipper->GetMClipperAngles()[idxClipper][1] <= 0.0)
            adjustments[1] = 360.0;
            else
            adjustments[1] = self->clipper->GetMClipperAngles()[idxClipper][1] - 0.1;
            self->clipperActors.at(idxClipper)->GetProperty()->SetColor(0,1,0);
            self->clipper->SetMClipperAngles(adjustments, idxClipper);
            self->CtrPlanesPlacer();

        } else if (key == "space") {

            self->clipperActors.at(idxClipper)->GetProperty()->SetColor(1,0,0);
            self->ManualCutterCallBack();
            self->clipper->SetMClipperSeeds(self->pickedCutterSeeds, self->m_Controls.comboBox->currentIndex());
            self->CtrPlanesPlacer();

        } else if (key == "Delete") {

            vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
            vtkSmartPointer<vtkPoints> points = self->pickedCutterSeeds->GetPoints();
            for (int i=0; i<points->GetNumberOfPoints()-1; i++)
            newPoints->InsertNextPoint(points->GetPoint(i));
            self->pickedCutterSeeds->SetPoints(newPoints);
            self->clipper->SetMClipperSeeds(self->pickedCutterSeeds, self->m_Controls.comboBox->currentIndex());
            self->CtrPlanesPlacer();

        } else if (key == "a" || key == "A") {

            self->m_Controls.label->setText(" Automatic ");
            QString style = "QLabel {border-width:1px; border-color:black; border-radius:10px; background-color:yellow;}";
            self->m_Controls.label->setStyleSheet(style);
            self->clipper->SetToAutomaticClipperMode(idxClipper);
            self->clipperActors.at(idxClipper)->GetProperty()->SetColor(1,1,0);
            self->CtrPlanesPlacer();

        } else if (key == "o" || key == "O") {

            if (self->surfActor->GetProperty()->GetOpacity() == 1.0)
            self->surfActor->GetProperty()->SetOpacity(0.5);
            else {
                self->surfActor->GetProperty()->SetOpacity(1.0);
                for (unsigned int i=0; i<self->clipperActors.size(); i++)
                self->clipperActors.at(i)->GetProperty()->SetOpacity(1.0);
            }//_if
            self->m_Controls.widget_1->GetRenderWindow()->Render();

        }//_if_key

    } else if (!self->IsPointSelectionControlsAvailable()) {

        if (key == "r" || key =="R") {

            //Clear renderer
            self->renderer->RemoveAllViewProps();
            self->clipperActors.clear();
            self->pickedSeedLabels.clear();
            self->InitialisePickerObjects();
            self->Visualiser();

            //Reset clipper
            if (self->m_Controls.button_man2_clippers->isEnabled() && self->m_Controls.button_man3_clipseg->isEnabled())
            self->clipper = std::unique_ptr<CemrgAtriaClipper>(new CemrgAtriaClipper(self->directory, self->surface));
            else
            self->clipper->ResetCtrLinesClippingPlanes();

            //Reset controls
            self->SetAutomaticModeButtons(self->automaticPipeline);
            self->SetManualModeButtons(self->automaticPipeline);
            self->m_Controls.comboBox->clear();
            self->m_Controls.slider->setRange(0,2);
            self->m_Controls.slider->setValue(0);
            self->m_Controls.spinBox->setValue(2.0);
            self->m_Controls.slider->setEnabled(false);
            self->m_Controls.spinBox->setEnabled(false);
            self->m_Controls.comboBox->setEnabled(false);
            self->m_Controls.checkBox->setEnabled(true);
            self->m_Labels.radioBtn_LSPV->setEnabled(true);
            self->m_Labels.radioBtn_LIPV->setEnabled(true);
            self->m_Labels.radioBtn_ignore->setEnabled(true);
            self->m_Labels.radioBtn_RSPV->setEnabled(true);
            self->m_Labels.radioBtn_RIPV->setEnabled(true);
            self->m_Labels.radioBtn_discard->setEnabled(true);
            self->m_Labels.radioBtn_LAA->setEnabled(true);

            //Setup shortcuts
            vtkSmartPointer<vtkTextActor> txtActor = vtkSmartPointer<vtkTextActor>::New();
            std::string shortcuts = self->GetShortcuts();
            txtActor->SetInput(shortcuts.c_str());
            txtActor->GetTextProperty()->SetFontSize(14);
            txtActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
            self->renderer->AddActor2D(txtActor);

        }//_if
    }//_if_main

}

// helper functions
void AtrialFibresClipperView::InitialisePickerObjects(){
    pickedSeedIds = vtkSmartPointer<vtkIdList>::New();
    pickedSeedIds->Initialize();
    pickedLineSeeds = vtkSmartPointer<vtkPolyData>::New();
    pickedLineSeeds->Initialize();
    pickedLineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
    if(automaticPipeline){
        corridorCount=0;
        corridorSeedIds = vtkSmartPointer<vtkIdList>::New();
        corridorSeedIds->Initialize();
        corridorLineSeeds = vtkSmartPointer<vtkPolyData>::New();
        corridorLineSeeds->Initialize();
        corridorLineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
    }
    pickedCutterSeeds = vtkSmartPointer<vtkPolyData>::New();
    pickedCutterSeeds->Initialize();
    pickedCutterSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());
}

bool AtrialFibresClipperView::IsPointSelectionControlsAvailable(){
    bool res = false;
    if(automaticPipeline){
        res = m_Controls.button_auto1_savelabels->isEnabled();
    } else{
        res = m_Controls.button_man1_ctrlines->isEnabled();
    }
    return res;
}

bool AtrialFibresClipperView::IsClipperManualControlsAvailable(){
    bool res = false;
    if(!automaticPipeline){
        res = clipper->GetCentreLinePolyPlanes().size() != 0 && !m_Controls.button_man1_ctrlines->isEnabled() && m_Controls.button_man3_clipseg->isEnabled();
    }
    return res;
}

std::string AtrialFibresClipperView::GetShortcuts(){
    std::string res = "";
    if(automaticPipeline){
        res += "POINT SELECTION:\n\tSpace: select landmark\n\tDelete: remove point";
        res += "\nINTER PV CORRECTION:\n\tX: add seed point\n\tD: remove point";
    } else{
        res += "R: reset centrelines\nSpace: add seed point\nDelete: remove seed point";
    }
    return res;
}

void AtrialFibresClipperView::SetManualModeButtons(bool b){
    m_Controls.button_man1_ctrlines->setVisible(b);
    m_Controls.button_man2_clippers->setVisible(b);
    m_Controls.button_man3_clipseg->setVisible(b);
    m_Controls.label_man->setVisible(b);
}

void AtrialFibresClipperView::SetAutomaticModeButtons(bool b){
    m_Controls.button_auto1_savelabels->setVisible(b);
    m_Controls.button_auto2_clippers->setVisible(b);
    m_Controls.button_auto3_spacing->setVisible(b);
    m_Controls.label_auto->setVisible(b);
}

void AtrialFibresClipperView::UserSelectPvLabel(){
    int dialogCode = inputs->exec();
    QRect screenGeometry = QApplication::desktop()->screenGeometry();
    int x = (screenGeometry.width() - inputs->width()) / 2;
    int y = (screenGeometry.height() - inputs->height()) / 2;
    inputs->move(x,y);

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {

        if (m_Labels.radioBtn_LSPV->isChecked()) {
            pickedSeedLabels.push_back(11); // LSPV
            m_Labels.radioBtn_LSPV->setEnabled(false);
        } else if (m_Labels.radioBtn_LIPV->isChecked()) {
            pickedSeedLabels.push_back(13); // LIPV
            m_Labels.radioBtn_LIPV->setEnabled(false);
        } else if (m_Labels.radioBtn_ignore->isChecked()) {
            pickedSeedLabels.push_back(14); // Ignore
        } else if (m_Labels.radioBtn_RSPV->isChecked()) {
            pickedSeedLabels.push_back(15); // RSPV
            m_Labels.radioBtn_RSPV->setEnabled(false);
        } else if (m_Labels.radioBtn_RIPV->isChecked()) {
            pickedSeedLabels.push_back(17); // RIPV
            m_Labels.radioBtn_RIPV->setEnabled(false);
        } else if (m_Labels.radioBtn_discard->isChecked()) {
            pickedSeedLabels.push_back(18); // Discard
        } else if (m_Labels.radioBtn_LAA->isChecked()) {
            pickedSeedLabels.push_back(19); // LAAP_1
            m_Labels.radioBtn_LAA->setEnabled(false);
        } else
        pickedSeedLabels.push_back(21);
        m_Labels.radioBtn_default->setChecked(true);

    } else if (dialogCode == QDialog::Rejected) {
        inputs->close();
    }//_if
}

void AtrialFibresClipperView::IgnoreLabel(std::vector<int> ignoredIds){
    MITK_INFO << "[IgnoreLabel] Ignoring labels in shell";
    QString prodPathOut = AtrialFibresClipperView::directory + mitk::IOUtil::GetDirectorySeparator();

    vtkFloatArray *cellScalars = vtkFloatArray::New();
    vtkFloatArray *pointScalars = vtkFloatArray::New();
    cellScalars = vtkFloatArray::SafeDownCast(surface->GetVtkPolyData()->GetCellData()->GetScalars());
    pointScalars = vtkFloatArray::SafeDownCast(surface->GetVtkPolyData()->GetPointData()->GetScalars());

    double ptLabel;
    std::cout << "Ignoring: ";
    for (vtkIdType labelIdx = 0; labelIdx < ignoredIds.size(); labelIdx++) {
        ptLabel = pointScalars->GetTuple1(labelIdx);
        std::cout << ptLabel << "," ;
        double s;
        for (vtkIdType vId = 0; vId < surface->GetVtkPolyData()->GetNumberOfPoints(); vId++) {
            s = pointScalars->GetTuple1(vId);

            if(s==ptLabel){
                pointScalars->SetTuple1(vId, 1);

                vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
                cellIds->Initialize();
                surface->GetVtkPolyData()->GetPointCells(vId, cellIds);

                for (vtkIdType ix = 0; ix < cellIds->GetNumberOfIds() ; ix++) {
                    cellScalars->SetTuple1(cellIds->GetId(ix), 1);
                }
            }
        }
    }

    surface->GetVtkPolyData()->GetPointData()->SetScalars(pointScalars);
    surface->GetVtkPolyData()->GetCellData()->SetScalars(cellScalars);

    MITK_INFO << "[IgnoreLabel] Set new scalar vector into surface.";
    mitk::IOUtil::Save(surface, (prodPathOut+AtrialFibresClipperView::fileName).toStdString());
    MITK_INFO << "[IgnoreLabel] Saved surface";

}

void AtrialFibresClipperView::DiscardUnwantedLabel(std::vector<int> discardedIds){
    vtkFloatArray *pointScalars = vtkFloatArray::New();
    pointScalars = vtkFloatArray::SafeDownCast(surface->GetVtkPolyData()->GetPointData()->GetScalars());

    double ptLabel;
    for (vtkIdType labelIdx = 0; labelIdx < discardedIds.size(); labelIdx++) {
        ptLabel = pointScalars->GetTuple1(labelIdx);
    }
}

void AtrialFibresClipperView::PrintCorridorIds(){
    if (corridorSeedIds->GetNumberOfIds() > 0) {
        int numCorridor = corridorSeedIds->GetNumberOfIds();
        std::cout << "[PriontCorridorIds] Points in corridor:" << '\n';
        for (int ix = 0; ix < numCorridor; ix++) {
            std::cout << "[" << ix << "] " << corridorSeedIds->GetId(ix) << '\n';
        }
    } else {
        std::cout << "[PriontCorridorIds] No points in corridor" << '\n';
    }
}
