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

//Blueberry
#include <berryIWorkbenchPage.h>
#include <berryFileEditorInput.h>

//Qmitk
#include <QmitkIOUtil.h>
#include <mitkProgressBar.h>
#include <mitkIDataStorageService.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkDataStorageEditorInput.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkBoundingObject.h>
#include <mitkCuboid.h>
#include <mitkAffineImageCropperInteractor.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkUnstructuredGrid.h>
#include "mitkPluginActivator.h"
#include "WallThicknessCalculationsView.h"
#include "WallThicknessCalculationsClipperView.h"

//Micro services
#include <usModuleRegistry.h>
#include <arpa/inet.h>

//VTK
#include <vtkDecimatePro.h>
#include <vtkFieldData.h>

//ITK
#include <itkAddImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkThresholdImageFilter.h>

//Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>

//MyCemrgLib
#include <CemrgImageUtils.h>
#include <CemrgCommandLine.h>
#include <CemrgMeasure.h>


const std::string WallThicknessCalculationsView::VIEW_ID = "my.cemrgproject.views.wallthicknesscalculationsview";

void WallThicknessCalculationsView::CreateQtPartControl(QWidget *parent) {
    
    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, SIGNAL(clicked()), this, SLOT(LoadDICOM()));
    connect(m_Controls.button_2, SIGNAL(clicked()), this, SLOT(ProcessIMGS()));
    connect(m_Controls.button_3, SIGNAL(clicked()), this, SLOT(SegmentIMGS()));
    connect(m_Controls.button_4, SIGNAL(clicked()), this, SLOT(ClipperPV()));
    connect(m_Controls.button_5, SIGNAL(clicked()), this, SLOT(MorphologyAnalysis()));
    connect(m_Controls.button_6, SIGNAL(clicked()), this, SLOT(ThicknessAnalysis()));
    connect(m_Controls.button_r, SIGNAL(clicked()), this, SLOT(Reset()));
    
    //Set visibility of buttons
    m_Controls.button_2_1->setVisible(false);
    m_Controls.button_2_2->setVisible(false);
    m_Controls.button_2_3->setVisible(false);
    m_Controls.button_2_4->setVisible(false);
    connect(m_Controls.button_2_1, SIGNAL(clicked()), this, SLOT(ConvertNII()));
    connect(m_Controls.button_2_2, SIGNAL(clicked()), this, SLOT(CropIMGS()));
    connect(m_Controls.button_2_3, SIGNAL(clicked()), this, SLOT(ResampIMGS()));
    connect(m_Controls.button_2_4, SIGNAL(clicked()), this, SLOT(ApplyFilter()));
    m_Controls.button_3_1->setVisible(false);
    m_Controls.button_3_2->setVisible(false);
    m_Controls.button_3_3->setVisible(false);
    connect(m_Controls.button_3_1, SIGNAL(clicked()), this, SLOT(SelectROI()));
    connect(m_Controls.button_3_2, SIGNAL(clicked()), this, SLOT(SelectLandmarks()));
    connect(m_Controls.button_3_3, SIGNAL(clicked()), this, SLOT(CombineSegs()));
    m_Controls.button_6_1->setVisible(false);
    m_Controls.button_6_2->setVisible(false);
    connect(m_Controls.button_6_1, SIGNAL(clicked()), this, SLOT(ConvertNRRD()));
    connect(m_Controls.button_6_2, SIGNAL(clicked()), this, SLOT(ThicknessCalculator()));
}

void WallThicknessCalculationsView::SetFocus() {

    m_Controls.button_1->setFocus();
}

void WallThicknessCalculationsView::LoadDICOM() {
    
    //Use MITK DICOM editor
    QString editor_id = "org.mitk.editors.dicomeditor";
    berry::IEditorInput::Pointer input(new berry::FileEditorInput(QString()));
    this->GetSite()->GetPage()->OpenEditor(input, editor_id);
}

void WallThicknessCalculationsView::ProcessIMGS() {
    
    //Toggle visibility of buttons
    if (m_Controls.button_2_1->isVisible()) {
        m_Controls.button_2_1->setVisible(false);
        m_Controls.button_2_2->setVisible(false);
        m_Controls.button_2_3->setVisible(false);
        m_Controls.button_2_4->setVisible(false);
    } else {
        m_Controls.button_2_1->setVisible(true);
        m_Controls.button_2_2->setVisible(true);
        m_Controls.button_2_3->setVisible(true);
        m_Controls.button_2_4->setVisible(true);
    }
}

void WallThicknessCalculationsView::ConvertNII() {
    
    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() < 1) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please load and select images from the Data Manager before starting this step!");
        return;
    }//_if
    
    //Ask the user for a dir to store data
    if (directory.isEmpty()) {
        directory = QFileDialog::getExistingDirectory(
                    NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
                    QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            return;
        }//_if
    }
    
    //Generic Conversion to nii
    int ctr = 0;
    QString path;
    this->BusyCursorOn();
    mitk::ProgressBar::GetInstance()->AddStepsToDo(nodes.size());
    foreach (mitk::DataNode::Pointer node, nodes) {
        mitk::BaseData::Pointer data = node->GetData();
        if (data) {
            //Test if this data item is an image
            mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
            if (image) {
                path = directory + mitk::IOUtil::GetDirectorySeparator() + "dcm-" + QString::number(ctr++) + ".nii";
                mitk::IOUtil::Save(image, path.toStdString());
                this->GetDataStorage()->Remove(node);
            } else
                return;
        } else
            return;
        mitk::ProgressBar::GetInstance()->Progress();
    }//_for
    nodes.clear();
    this->BusyCursorOff();
    
    //Load first item
    ctr = 0;
    path = directory + mitk::IOUtil::GetDirectorySeparator() + "dcm-" + QString::number(ctr) + ".nii";
    mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
}

void WallThicknessCalculationsView::CropIMGS() {
    
    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please select an image from the Data Manager to perform cropping!");
        return;
    }//_if
    
    //Check to cut now or not
    if (m_Controls.button_2_2->text() == QString::fromStdString("Are you done?")) {
        
        QString path;
        //Ask the user for a dir to locate data
        if (directory.isEmpty()) {
            directory = QFileDialog::getExistingDirectory(
                        NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
                        QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
            if (directory.isEmpty() || directory.simplified().contains(" ")) {
                QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
                directory = QString();
                return;
            }//_if
        }
        
        //Cut selected image
        this->BusyCursorOn();
        mitk::ProgressBar::GetInstance()->AddStepsToDo(2);
        mitk::Image::Pointer outputImage = CemrgImageUtils::CropImage();
        path = directory + mitk::IOUtil::GetDirectorySeparator() + CemrgImageUtils::GetImageNode()->GetName().c_str() + ".nii";
        mitk::IOUtil::Save(outputImage, path.toStdString());
        this->BusyCursorOff();
        
        //Update datastorage
        AddToStorage(CemrgImageUtils::GetImageNode()->GetName(), outputImage);
        this->GetDataStorage()->Remove(CemrgImageUtils::GetImageNode());
        this->GetDataStorage()->Remove(CemrgImageUtils::GetCuttingNode());
        mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
        
        m_Controls.button_2_2->setText("Crop Images");
        return;
    }//_if
    
    //Prepare cutting cuboid
    mitk::Cuboid::Pointer cuttingObject = mitk::Cuboid::New();
    mitk::DataNode::Pointer cuttingNode = mitk::DataNode::New();
    cuttingNode->SetData(cuttingObject);
    cuttingNode->SetProperty("opacity", mitk::FloatProperty::New(0.4));
    cuttingNode->SetProperty("color", mitk::ColorProperty::New(1.0, 1.0, 0.0));
    cuttingNode->SetProperty("name", mitk::StringProperty::New("Cropper"));
    this->GetDataStorage()->Add(cuttingNode);
    
    //Mouse interactions
    mitk::AffineImageCropperInteractor::Pointer affineDataInteractor = mitk::AffineImageCropperInteractor::New();
    affineDataInteractor->LoadStateMachine("ClippingPlaneInteraction3D.xml", us::ModuleRegistry::GetModule("MitkDataTypesExt"));
    affineDataInteractor->SetEventConfig("CropperDeformationConfig.xml", us::ModuleRegistry::GetModule("MitkDataTypesExt"));
    affineDataInteractor->SetDataNode(cuttingNode);
    cuttingNode->SetBoolProperty("pickable", true);
    
    //Fit the cuboid to the image
    mitk::Image::Pointer imageToCut;
    mitk::BoundingObject::Pointer cuttingCube;
    mitk::DataNode::Pointer imageNode = nodes.at(0);
    mitk::BaseData::Pointer data = imageNode->GetData();
    cuttingCube = dynamic_cast<mitk::BoundingObject*>(cuttingNode->GetData());
    if (data) {
        //Test if this data item is an image
        imageToCut = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (imageToCut)
            cuttingCube->FitGeometry(imageToCut->GetGeometry());
        else return;
    } else return;
    
    //To be used for actual cutting
    CemrgImageUtils::SetImageToCut(imageToCut);
    CemrgImageUtils::SetCuttingCube(cuttingCube);
    CemrgImageUtils::SetImageNode(imageNode);
    CemrgImageUtils::SetCuttingNode(cuttingNode);
    m_Controls.button_2_2->setText("Are you done?");
}

void WallThicknessCalculationsView::ResampIMGS() {
    
    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please select an image from the Data Manager to perform downsampling!");
        return;
    }
    
    //Ask the user for a dir to store data
    if (directory.isEmpty()) {
        directory = QFileDialog::getExistingDirectory(
                    NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
                    QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            return;
        }//_if
    }
    
    //Find the selected node
    QString path;
    mitk::DataNode::Pointer imgNode = nodes.at(0);
    mitk::BaseData::Pointer data = imgNode->GetData();
    if (data) {
        
        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image) {
            
            bool ok;
            int factor = QInputDialog::getInt(NULL, tr("Downsampling"), tr("By factor of:"), 3, 1, 5, 1, &ok);
            if (ok) {
                
                //Downsample selected image
                this->BusyCursorOn();
                mitk::ProgressBar::GetInstance()->AddStepsToDo(1);
                mitk::Image::Pointer outputImage = CemrgImageUtils::Downsample(image, factor);
                path = directory + mitk::IOUtil::GetDirectorySeparator() + imgNode->GetName().c_str() + ".nii";
                mitk::IOUtil::Save(outputImage, path.toStdString());
                this->BusyCursorOff();
                //Update datastorage
                AddToStorage(imgNode->GetName(), outputImage);
                this->GetDataStorage()->Remove(imgNode);
                mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());

            }//_if
        } else
            return;
    } else
        return;
}

void WallThicknessCalculationsView::ApplyFilter() {
    
}

void WallThicknessCalculationsView::SegmentIMGS() {
    
    //Toggle visibility of buttons
    if (m_Controls.button_3_1->isVisible()) {
        m_Controls.button_3_1->setVisible(false);
        m_Controls.button_3_2->setVisible(false);
        m_Controls.button_3_3->setVisible(false);
        return;
    } else {
        m_Controls.button_3_1->setVisible(true);
        m_Controls.button_3_2->setVisible(true);
        m_Controls.button_3_3->setVisible(true);
    }//_if
    
    int reply = QMessageBox::question(
                NULL, "Question", "Do you have a segmentation to load?",
                QMessageBox::Yes, QMessageBox::No);
    
    if (reply == QMessageBox::Yes) {
        
        QString path = QFileDialog::getOpenFileName(
                    NULL, "Open Segmentation file", mitk::IOUtil::GetProgramPath().c_str(),
                    QmitkIOUtil::GetFileOpenFilterString());
        if (path.isEmpty())
            return;
        mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
        mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());

        //Restore image name
        char sep = mitk::IOUtil::GetDirectorySeparator();
        fileName = path.mid(path.lastIndexOf(sep) + 1);
        
    } else {
        //Show the plugin
        QMessageBox::information(
                    NULL, "Attention",
                    "After finalising the segmentation, save the image and load back using this same button!");
        this->GetSite()->GetPage()->ShowView("org.mitk.views.segmentation");
    }//_if
}

void WallThicknessCalculationsView::SelectROI() {
    
}

void WallThicknessCalculationsView::SelectLandmarks() {
    
}

void WallThicknessCalculationsView::CombineSegs() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 2) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please select both segmentations from the Data Manager to create a surface for clipping!");
        return;
    }//_if

    //Ask the user for a dir to store data
    if (directory.isEmpty()) {
        directory = QFileDialog::getExistingDirectory(
                    NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
                    QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            return;
        }//_if
    }//_if

    //Find the selected nodes
    mitk::DataNode::Pointer segNode_0 = nodes.at(0);
    mitk::BaseData::Pointer data_0 = segNode_0->GetData();
    mitk::DataNode::Pointer segNode_1 = nodes.at(1);
    mitk::BaseData::Pointer data_1 = segNode_1->GetData();
    if (data_0 && data_1) {
        //Test if this data item is an image
        mitk::Image::Pointer image_0 = dynamic_cast<mitk::Image*>(data_0.GetPointer());
        mitk::Image::Pointer image_1 = dynamic_cast<mitk::Image*>(data_1.GetPointer());
        if (image_0 && image_1) {

            //Add images
            typedef itk::Image<short, 3> ImageType;
            typedef itk::AddImageFilter<ImageType, ImageType, ImageType> AddFilterType;
            QString path = directory + mitk::IOUtil::GetDirectorySeparator() + "segmentation.nii";
            //Cast seg to ITK format
            ImageType::Pointer itkImage_0 = ImageType::New();
            ImageType::Pointer itkImage_1 = ImageType::New();
            CastToItkImage(image_0, itkImage_0);
            CastToItkImage(image_1, itkImage_1);
            AddFilterType::Pointer addFilter = AddFilterType::New();
            addFilter->SetInput1(itkImage_0);
            addFilter->SetInput2(itkImage_1);
            addFilter->UpdateLargestPossibleRegion();
            typedef itk::ImageRegionIteratorWithIndex<ImageType> ItType;
            ItType itLbl(addFilter->GetOutput(), addFilter->GetOutput()->GetRequestedRegion());
            for (itLbl.GoToBegin(); !itLbl.IsAtEnd(); ++itLbl)
                if ((int)itLbl.Get() != 0)
                    itLbl.Set(1);

            //Save result
            mitk::Image::Pointer resImage = mitk::ImportItkImage(addFilter->GetOutput())->Clone();
            mitk::IOUtil::Save(resImage, path.toStdString());
            AddToStorage("segmentation", resImage);
            fileName = "segmentation.nii";
        } else
            return;
    } else
        return;
}

void WallThicknessCalculationsView::ClipperPV() {
    
    //Ask the user for a dir to store data
    if (directory.isEmpty()) {
        directory = QFileDialog::getExistingDirectory(
                    NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
                    QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            return;
        }//_if
    }

    //Show the plugin
    this->GetSite()->GetPage()->ResetPerspective();
    WallThicknessCalculationsClipperView::SetDirectoryFile(directory, fileName);
    this->GetSite()->GetPage()->ShowView("my.cemrgproject.views.wallthicknesscalculationsclipperview");
}

void WallThicknessCalculationsView::MorphologyAnalysis() {

    //Ask the user for a dir to store data
    if (directory.isEmpty()) {
        directory = QFileDialog::getExistingDirectory(
                    NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
                    QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            return;
        }//_if
    }//_if

    try {

        QString path = directory + mitk::IOUtil::GetDirectorySeparator() + "AnalyticBloodpool.nii";
        mitk::Image::Pointer analyticImage = mitk::IOUtil::LoadImage(path.toStdString());

        if (analyticImage) {

            //Loop through labelled image
            typedef itk::Image<short, 3> ImageType;
            typedef itk::ImageRegionIteratorWithIndex<ImageType> ItType;
            ImageType::Pointer analyticItkImage = ImageType::New();
            CastToItkImage(analyticImage, analyticItkImage);
            ItType itLbl(analyticItkImage, analyticItkImage->GetRequestedRegion());
            for (itLbl.GoToBegin(); !itLbl.IsAtEnd(); ++itLbl)
                if ((int)itLbl.Get() == APPENDAGECUT || (int)itLbl.Get() == APPENDAGEUNCUT)
                    itLbl.Set(0);

            //Relabel the components to separate bloodpool and appendage based on the number of voxels
            typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> ConnectedComponentImageFilterType;
            ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
            connected->SetInput(analyticItkImage);
            connected->Update();
            typedef itk::RelabelComponentImageFilter<ImageType, ImageType> RelabelFilterType;
            RelabelFilterType::Pointer relabeler = RelabelFilterType::New();
            relabeler->SetInput(connected->GetOutput());
            relabeler->Update();

            //Keep the selected labels
            const int bpLabel = 1;
            const int apLabel = 2;
            typedef itk::ThresholdImageFilter<ImageType> ThresholdImageFilterType;
            ThresholdImageFilterType::Pointer thresholdFilter1 = ThresholdImageFilterType::New();
            thresholdFilter1->SetInput(analyticItkImage);
            thresholdFilter1->ThresholdOutside(bpLabel, bpLabel);
            thresholdFilter1->SetOutsideValue(0);
            ThresholdImageFilterType::Pointer thresholdFilter2 = ThresholdImageFilterType::New();
            thresholdFilter2->SetInput(relabeler->GetOutput());
            thresholdFilter2->ThresholdOutside(apLabel, apLabel);
            thresholdFilter2->SetOutsideValue(0);

            //Import to MITK image
            mitk::Image::Pointer bp = mitk::ImportItkImage(thresholdFilter1->GetOutput());
            mitk::Image::Pointer ap = mitk::ImportItkImage(thresholdFilter2->GetOutput());

            //Ask for user input to set the parameters
            QDialog* inputs = new QDialog(0,0);
            m_UIMeshing.setupUi(inputs);
            connect(m_UIMeshing.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
            connect(m_UIMeshing.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
            int dialogCode = inputs->exec();

            //Act on dialog return code
            if (dialogCode == QDialog::Accepted) {

                bool ok1, ok2, ok3, ok4, ok5;
                int iter = m_UIMeshing.lineEdit_1->text().toInt(&ok1);
                float th = m_UIMeshing.lineEdit_2->text().toFloat(&ok2);
                int blur = m_UIMeshing.lineEdit_3->text().toInt(&ok3);
                int smth = m_UIMeshing.lineEdit_4->text().toInt(&ok4);
                float ds = m_UIMeshing.lineEdit_5->text().toFloat(&ok5);

                //Set default values
                if (!ok1 || !ok2 || !ok3 || !ok4 || !ok5)
                    QMessageBox::warning(NULL, "Attention", "Reverting to default parameters!");
                if(!ok1) iter = 1;
                if(!ok2) th   = 0.5;
                if(!ok3) blur = 0;
                if(!ok4) smth = 0;
                if(!ok5) ds   = 0.99;
                //_if

                this->BusyCursorOn();
                QString path1 = directory + mitk::IOUtil::GetDirectorySeparator() + "bp.nii";
                mitk::IOUtil::Save(bp, path1.toStdString());
                mitk::ProgressBar::GetInstance()->AddStepsToDo(4);
                std::unique_ptr<CemrgCommandLine> cmd1(new CemrgCommandLine());
                QString output = cmd1->ExecuteSurf(directory, path1, iter, th, blur, smth);
                QMessageBox::information(NULL, "Attention", "Command Line Operations (Bloodpool) Finished!");
                //Decimate the mesh to visualise
                mitk::Surface::Pointer shell = mitk::IOUtil::LoadSurface(output.toStdString());
                vtkSmartPointer<vtkDecimatePro> deci = vtkSmartPointer<vtkDecimatePro>::New();
                deci->SetInputData(shell->GetVtkPolyData());
                deci->SetTargetReduction(ds);
                deci->PreserveTopologyOn();
                deci->Update();
                shell->SetVtkPolyData(deci->GetOutput());
                mitk::Surface::Pointer surfLA = shell->Clone();
                QString path2 = directory + mitk::IOUtil::GetDirectorySeparator() + "ap.nii";
                mitk::IOUtil::Save(ap, path2.toStdString());
                mitk::ProgressBar::GetInstance()->AddStepsToDo(4);
                std::unique_ptr<CemrgCommandLine> cmd2(new CemrgCommandLine());
                output = cmd2->ExecuteSurf(directory, path2, iter, th, blur, smth);
                QMessageBox::information(NULL, "Attention", "Command Line Operations (Appendage) Finished!");
                //Decimate the mesh to visualise
                shell = mitk::IOUtil::LoadSurface(output.toStdString());
                deci->SetInputData(shell->GetVtkPolyData());
                deci->SetTargetReduction(ds);
                deci->PreserveTopologyOn();
                deci->Update();
                shell->SetVtkPolyData(deci->GetOutput());
                mitk::Surface::Pointer surfAP = shell->Clone();
                this->BusyCursorOff();

                //Volume and surface calculations
                std::unique_ptr<CemrgMeasure> morphAnal = std::unique_ptr<CemrgMeasure>(new CemrgMeasure());
                double surfceLA = morphAnal->calcSurfaceMesh(surfLA);
                double volumeLA = morphAnal->calcVolumeMesh(surfLA);
                double surfceAP = morphAnal->calcSurfaceMesh(surfAP);
                double volumeAP = morphAnal->calcVolumeMesh(surfAP);

                //Store in text file
                ofstream morphResult;
                QString morphPath = directory + mitk::IOUtil::GetDirectorySeparator() + "morphResults.txt";
                morphResult.open(morphPath.toStdString(), std::ios_base::app);
                morphResult << "SA" << " " << surfceLA << "\n";
                morphResult << "VA" << " " << volumeLA << "\n";
                morphResult << "SP" << " " << surfceAP << "\n";
                morphResult << "VP" << " " << volumeAP << "\n";
                morphResult.close();

                //Tidy up data files
                bool ok;
                QString morphFile = QInputDialog::getText(NULL, tr("Save As"), tr("File Name:"), QLineEdit::Normal, "morphResults.txt", &ok);
                if (ok && !morphFile.isEmpty() && morphFile.endsWith(".txt")) {
                    QString morphNewPath = directory + mitk::IOUtil::GetDirectorySeparator() + morphFile;
                    int result = rename(morphPath.toStdString().c_str(), morphNewPath.toStdString().c_str());
                    if (result == 0)
                        QMessageBox::information(NULL, "Attention", "File name was changed successfully!");
                    else
                        QMessageBox::warning(NULL, "Attention", "Wrong file name input. Saved with default name!");
                } else
                    QMessageBox::warning(NULL, "Attention", "Wrong file name input. Saved with default name!");
                inputs->deleteLater();

            } else if(dialogCode == QDialog::Rejected) {
                inputs->close();
                inputs->deleteLater();
                this->GetSite()->GetPage()->ResetPerspective();
                return;
            }//_if
        }//_if

    } catch (...) {
        QMessageBox::critical(NULL, "Attention", "Cropped segmentation image was not found!");
        return;
    }//_try
}

void WallThicknessCalculationsView::ThicknessAnalysis() {

    //Toggle visibility of buttons
    if (m_Controls.button_6_1->isVisible()) {
        m_Controls.button_6_1->setVisible(false);
        m_Controls.button_6_2->setVisible(false);
        return;
    } else {
        m_Controls.button_6_1->setVisible(true);
        m_Controls.button_6_2->setVisible(true);
    }//_if
}

void WallThicknessCalculationsView::ConvertNRRD() {

    //Ask the user for a dir to store data
    if (directory.isEmpty()) {
        directory = QFileDialog::getExistingDirectory(
                    NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
                    QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            return;
        }//_if
    }

    try {

        bool ok;
        QString path = directory + mitk::IOUtil::GetDirectorySeparator() + "PVeinsCroppedImage.nii";
        mitk::Image::Pointer clippedImage = mitk::IOUtil::LoadImage(path.toStdString());
        QString tmpFileName = QInputDialog::getText(NULL, tr("Save Segmentation As"), tr("File Name:"), QLineEdit::Normal, ".nrrd", &ok);

        if (ok && !tmpFileName.isEmpty() && tmpFileName.endsWith(".nrrd") && tmpFileName != ".nrrd") {
            path = directory + mitk::IOUtil::GetDirectorySeparator() + tmpFileName;
            mitk::IOUtil::Save(clippedImage, path.toStdString());
            QMessageBox::information(NULL, "Attention", "Clipped Segmentation was successfully converted!");
            fileName = tmpFileName;
        } else {
            QMessageBox::warning(NULL, "Attention", "Please type a file name with the right extension (i.e. .nrrd)!");
            return;
        }//_fileName

    } catch (...) {
        QMessageBox::critical(NULL, "Attention", "Cropped segmentation image was not found!");
        return;
    }//_try
}

void WallThicknessCalculationsView::Browse() {
    
    QString para = "";
    para = QFileDialog::getOpenFileName(
                NULL, "Open text file containing parameters",
                directory, QmitkIOUtil::GetFileOpenFilterString());
    m_Thickness.lineEdit_1->setText(para);
}

void WallThicknessCalculationsView::ThicknessCalculator() {
    
    //Ask the user for a dir to store data
    if (directory.isEmpty()) {
        directory = QFileDialog::getExistingDirectory(
                    NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
                    QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            return;
        }//_if
    }
    
    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please select a segmentation from the Data Manager to calculate wall thickness!");
        return;
    }
    
    //Find the selected node
    QString templatePath;
    mitk::DataNode::Pointer segNode = nodes.at(0);
    mitk::BaseData::Pointer data = segNode->GetData();
    if (data) {
        
        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image) {

            try {

                //Convert image to right type
                mitk::Image::Pointer cvrImage = mitk::Image::New();
                itk::Image<unsigned char,3>::Pointer itkImage = itk::Image<unsigned char,3>::New();
                mitk::CastToItkImage(image, itkImage);
                mitk::CastToMitkImage(itkImage, cvrImage);
                
                //Access image volume
                mitk::ImagePixelReadAccessor<unsigned char,3> readAccess(cvrImage);
                unsigned char* pv = (unsigned char*)readAccess.GetData();

                //Image information
                const char* imgType = "unsigned fixed";
                int bitLength = sizeof(unsigned char) * 8;
                char scale[20]; sprintf(scale, "SCALE=2**0");
                unsigned int* dimensions = cvrImage->GetDimensions();
                double spacing[3]; cvrImage->GetGeometry()->GetSpacing().ToArray(spacing);
                const char* cpu = (htonl(47) == 47) ? "sun" : "decm"; //Big or little endian

                //Write to binary file
                FILE * pFile;
                std::string path = (directory + mitk::IOUtil::GetDirectorySeparator() + "converted.inr").toStdString();
                pFile = fopen(path.c_str(), "wb");

                //Prepare header of inr file
                std::ostringstream oss;
                oss << "#INRIMAGE-4#{" << "\n";
                oss << "XDIM=" << dimensions[0] << "\n";
                oss << "YDIM=" << dimensions[1] << "\n";
                oss << "ZDIM=" << dimensions[2] << "\n";
                oss << "VDIM=" << 1 << "\n";
                oss << "TYPE=" << imgType << "\n";
                oss << "PIXSIZE=" << bitLength <<" bits\n";
                oss << scale << "\n";
                oss << "CPU=" << cpu << "\n";
                oss << "VX=" << spacing[0] << "\n";
                oss << "VY=" << spacing[1] << "\n";
                oss << "VZ=" << spacing[2] << "\n";

                //Write content of header
                size_t hdrLength = oss.str().length();
                const void* b = oss.str().data();
                char* bufA = (char*)b;
                fwrite(bufA, 1, hdrLength, pFile);

                //Write filler of header
                std::size_t pos = oss.str().length();
                pos = pos % 256;
                if (pos > 252) {
                    for (std::size_t i=pos; i<256; i++)
                        fwrite("\n", 1, 1, pFile);
                    pos = 0;
                }//_if
                char bufB[257];
                bufB[0] = '\0';
                for (std::size_t i=pos; i < 252; i++)
                    strcat(bufB, "\n");
                strcat(bufB, "##}\n");
                fwrite(bufB, 1, strlen(bufB), pFile);

                //Write content of image
                int dims = dimensions[0] * dimensions[1] * dimensions[2];
                std::size_t sizeImg = dims * 1 * sizeof(unsigned char);
                char* data = (char*)pv;
                fwrite(data, 1, sizeImg, pFile);
                fclose(pFile);
                
                //Ask for user input to set the parameters
                QDialog* inputs = new QDialog(0,0);
                m_Thickness.setupUi(inputs);
                connect(m_Thickness.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
                connect(m_Thickness.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
                connect(m_Thickness.pushButton_1, SIGNAL(clicked()), this, SLOT(Browse()));
                int dialogCode = inputs->exec();
                
                //Act on dialog return code
                if (dialogCode == QDialog::Accepted) {
                    
                    //Checking input files
                    templatePath = m_Thickness.lineEdit_1->text();
                    
                    //Absolute path
                    if (templatePath.isEmpty()) {
                        QString aPath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "M3DLib";
#if defined(__APPLE__)
                        aPath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                                mitk::IOUtil::GetDirectorySeparator() + QString("M3DLib");
#endif
                        QMessageBox::warning(NULL, "Attention", "Reverting to default parameter file!");
                        templatePath = aPath + mitk::IOUtil::GetDirectorySeparator() + "param-template.par";
                    }//_if

                    //FileName checks
                    QString meshName = "CGALMesh";
                    if (fileName.endsWith(".nrrd"))
                        meshName = fileName.left(fileName.length() - 5);
                    
                    //Run Mesh3DTool
                    this->BusyCursorOn();
                    mitk::ProgressBar::GetInstance()->AddStepsToDo(1);
                    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
                    cmd->ExecuteCreateCGALMesh(directory, meshName, templatePath);
                    QMessageBox::information(NULL, "Attention", "Command Line Operations Finished!");
                    this->BusyCursorOff();
                    
                } else if(dialogCode == QDialog::Rejected) {
                    inputs->close();
                    inputs->deleteLater();
                }//_if
                
            } catch(mitk::Exception& e) {
                return;
            }//_try
        } else
            return;
    } else
        return;
}

void WallThicknessCalculationsView::Reset() {
    
    try {
        
        ctkPluginContext* context = mitk::PluginActivator::getContext();
        mitk::IDataStorageService* dss = 0;
        ctkServiceReference dsRef = context->getServiceReference<mitk::IDataStorageService>();
        
        if (dsRef)
            dss = context->getService<mitk::IDataStorageService>(dsRef);
        
        if (!dss) {
            MITK_WARN << "IDataStorageService service not available. Unable to close project.";
            context->ungetService(dsRef);
            return;
        }
        
        mitk::IDataStorageReference::Pointer dataStorageRef = dss->GetActiveDataStorage();
        if (dataStorageRef.IsNull()) {
            //No active data storage set (i.e. not editor with a DataStorageEditorInput is active).
            dataStorageRef = dss->GetDefaultDataStorage();
        }
        
        mitk::DataStorage::Pointer dataStorage = dataStorageRef->GetDataStorage();
        if (dataStorage.IsNull()) {
            MITK_WARN << "No data storage available. Cannot close project.";
            return;
        }
        
        //Check if we got the default datastorage and if there is anything else then helper object in the storage
        if (dataStorageRef->IsDefault() && dataStorage->GetSubset(
                    mitk::NodePredicateNot::New(
                        mitk::NodePredicateProperty::New("helper object", mitk::BoolProperty::New(true))))->empty())
            return;
        
        //Remove everything
        mitk::DataStorage::SetOfObjects::ConstPointer nodesToRemove = dataStorage->GetAll();
        dataStorage->Remove(nodesToRemove);
        
        //Remove the datastorage from the data storage service
        dss->RemoveDataStorageReference(dataStorageRef);
        
        //Close all editors with this data storage as input
        mitk::DataStorageEditorInput::Pointer dsInput(new mitk::DataStorageEditorInput(dataStorageRef));
        QList<berry::IEditorReference::Pointer> dsEditors =
                this->GetSite()->GetPage()->FindEditors(dsInput, QString(), berry::IWorkbenchPage::MATCH_INPUT);
        
        if (!dsEditors.empty()) {
            QList<berry::IEditorReference::Pointer> editorsToClose = dsEditors;
            this->GetSite()->GetPage()->CloseEditors(editorsToClose, false);
        }
        
    } catch (std::exception& e) {
        
        MITK_ERROR << "Exception caught during closing project: " << e.what();
        QMessageBox::warning(NULL, "Error", QString("An error occurred during Close Project: %1").arg(e.what()));
    }//_try
    
    //Clear project directory
    fileName = "";
    directory.clear();
    m_Controls.button_2_2->setText("Crop Images");
    this->GetSite()->GetPage()->ResetPerspective();
}

void WallThicknessCalculationsView::AddToStorage(std::string nodeName, mitk::BaseData* data) {
    
    if (!data) return;
    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(data);
    node->SetName(nodeName);
    this->GetDataStorage()->Add(node);
}

mitk::Surface::Pointer WallThicknessCalculationsView::ReadVTKMesh(std::string meshPath) {
    
    try {

        //Load the mesh
        mitk::Surface::Pointer surface = mitk::IOUtil::LoadSurface(meshPath);
        vtkSmartPointer<vtkPolyData> pd = surface->GetVtkPolyData();

        //Prepare points for MITK visualisation
        double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;
        for (int i=0; i<pd->GetNumberOfPoints(); i++) {
            double* point = pd->GetPoint(i);
            point[0] = -point[0];
            point[1] = -point[1];
            pd->GetPoints()->SetPoint(i, point);
            //Find mins and maxs
            if (i==0) {
                Xmin = point[0];
                Xmax = point[0];
                Ymin = point[1];
                Ymax = point[1];
                Zmin = point[2];
                Zmax = point[2];
            } else {
                if (point[0]<Xmin) Xmin = point[0];
                if (point[0]>Xmax) Xmax = point[0];
                if (point[1]<Ymin) Ymin = point[1];
                if (point[1]>Ymax) Ymax = point[1];
                if (point[2]<Zmin) Zmin = point[2];
                if (point[2]>Zmax) Zmax = point[2];
            }//_if
        }//_for

        double bounds[6] = {Xmin, Xmax, Ymin, Ymax, Zmin, Zmax};
        surface->GetGeometry()->SetBounds(bounds);
        return surface;

    } catch (...) {
        return mitk::Surface::New();
    }//_catch
}
