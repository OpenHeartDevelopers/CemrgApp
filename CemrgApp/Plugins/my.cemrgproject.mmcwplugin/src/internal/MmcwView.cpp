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
 * Motion Measurement of Cardiac Wall (MMCW) Plugin for MITK
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbench.h>
#include <berryIWorkbenchWindow.h>
#include <berryIWorkbenchPage.h>
#include <berryFileEditorInput.h>

// Qmitk
#include <QmitkIOUtil.h>
#include <mitkCoreObjectFactory.h>
#include <mitkIOUtil.h>
#include <mitkDataNode.h>
#include <mitkCuboid.h>
#include <mitkAffineImageCropperInteractor.h>
#include <mitkIDataStorageService.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkDataStorageEditorInput.h>
#include <mitkBaseData.h>
#include <mitkProgressBar.h>
#include <mitkImageReadAccessor.h>
#include "mitkPluginActivator.h"
#include "MmcwView.h"
#include "MmcwViewPlot.h"

//micro services
#include <usModuleRegistry.h>

// VTK
#include <vtkPolyData.h>

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QSignalMapper>
#include <QInputDialog>

// MyCemrgLib
#include <CemrgCommandLine.h>
#include <CemrgImageUtils.h>
#include <numeric>


const std::string MmcwView::VIEW_ID = "my.cemrgproject.views.mmcwview";

void MmcwView::CreateQtPartControl( QWidget *parent ) {
    
    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, SIGNAL(clicked()), this, SLOT(LoadDICOM()));
    connect(m_Controls.button_2, SIGNAL(clicked()), this, SLOT(ProcessIMGS()));
    connect(m_Controls.button_3, SIGNAL(clicked()), this, SLOT(SegmentIMGS()));
    connect(m_Controls.button_4, SIGNAL(clicked()), this, SLOT(CreateSurf()));
    connect(m_Controls.button_5, SIGNAL(clicked()), this, SLOT(TrackingIMGS()));
    connect(m_Controls.button_6, SIGNAL(clicked()), this, SLOT(LandmarkSelection()));
    connect(m_Controls.button_7, SIGNAL(clicked()), this, SLOT(Plotting()));
    connect(m_Controls.button_r, SIGNAL(clicked()), this, SLOT(Reset()));
    
    //Set visibility of buttons
    m_Controls.button_2_1->setVisible(false);
    m_Controls.button_2_2->setVisible(false);
    m_Controls.button_2_3->setVisible(false);
    connect(m_Controls.button_2_1, SIGNAL(clicked()), this, SLOT(ConvertNII()));
    connect(m_Controls.button_2_2, SIGNAL(clicked()), this, SLOT(CropinIMGS()));
    connect(m_Controls.button_2_3, SIGNAL(clicked()), this, SLOT(ResampIMGS()));
    
    //Set visibility of buttons
    m_Controls.button_5_1->setVisible(false);
    m_Controls.button_5_2->setVisible(false);
    m_Controls.button_5_3->setVisible(false);
    connect(m_Controls.button_5_1, SIGNAL(clicked()), this, SLOT(Tracking()));
    connect(m_Controls.button_5_2, SIGNAL(clicked()), this, SLOT(Applying()));
    connect(m_Controls.button_5_3, SIGNAL(clicked()), this, SLOT(Demoings()));

    //Temporal resolution
    timePoints = 0;
}

/**
 * @brief MmcwView::SetFocus
 */
void MmcwView::SetFocus() {
    
    m_Controls.button_1->setFocus();
}

void MmcwView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*source*/,
                                  const QList<mitk::DataNode::Pointer>& /*nodes*/) {
    
}

/**
 * @brief MmcwView::LoadDICOM
 */
void MmcwView::LoadDICOM() {
    
    //Use MITK DICOM editor
    QString editor_id = "org.mitk.editors.dicomeditor";
    berry::IEditorInput::Pointer input(new berry::FileEditorInput(QString()));
    this->GetSite()->GetPage()->OpenEditor(input, editor_id);
}

/**
 * @brief MmcwView::ProcessIMGS
 */
void MmcwView::ProcessIMGS() {
    
    //Toggle visibility of buttons
    if (m_Controls.button_2_1->isVisible()) {
        m_Controls.button_2_1->setVisible(false);
        m_Controls.button_2_2->setVisible(false);
        m_Controls.button_2_3->setVisible(false);
    } else {
        m_Controls.button_2_1->setVisible(true);
        m_Controls.button_2_2->setVisible(true);
        m_Controls.button_2_3->setVisible(true);
    }
}

/**
 * @brief MmcwView::ConvertNII
 */
void MmcwView::ConvertNII() {

    //Check for temporal resolution
    bool ok = true;
    if (timePoints == 0)
        timePoints = QInputDialog::getInt(NULL, tr("Time Points"), tr("Resolution:"), 10, 1, 100, 1, &ok);
    if (!ok) {
        QMessageBox::warning(NULL, "Attention", "Enter a correct value for the temporal resolution!");
        timePoints = 0;
        return;
    }//_if

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != timePoints) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please load and select all images from the Data Manager before starting this step!");
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

    //Order dicoms based on their cycle stages
    std::vector<int> indexNodes;
    std::string seriesDescription;
    foreach (mitk::DataNode::Pointer node, nodes) {
        node->GetData()->GetPropertyList()->GetStringProperty("dicom.series.SeriesDescription", seriesDescription);
        if (seriesDescription.find("90.0%")      != seriesDescription.npos) indexNodes.push_back(9);
        else if (seriesDescription.find("80.0%") != seriesDescription.npos) indexNodes.push_back(8);
        else if (seriesDescription.find("70.0%") != seriesDescription.npos) indexNodes.push_back(7);
        else if (seriesDescription.find("60.0%") != seriesDescription.npos) indexNodes.push_back(6);
        else if (seriesDescription.find("50.0%") != seriesDescription.npos) indexNodes.push_back(5);
        else if (seriesDescription.find("40.0%") != seriesDescription.npos) indexNodes.push_back(4);
        else if (seriesDescription.find("30.0%") != seriesDescription.npos) indexNodes.push_back(3);
        else if (seriesDescription.find("20.0%") != seriesDescription.npos) indexNodes.push_back(2);
        else if (seriesDescription.find("10.0%") != seriesDescription.npos) indexNodes.push_back(1);
        else if (seriesDescription.find("0.0%")  != seriesDescription.npos) indexNodes.push_back(0);
    }//_for
    //Sort indexes based on comparing values
    std::vector<int> index(indexNodes.size());
    std::iota(index.begin(), index.end(), 0);
    std::sort(index.begin(), index.end(), [&](int i1, int i2) {return indexNodes[i1]<indexNodes[i2];});

    //Warning for cases when order is not found
    size_t length1 = nodes.size();
    size_t length2 = indexNodes.size();
    if (length1 != length2) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Cannot find the order of images automatically. Revert to user order and selections in the data manager!");
        index.resize(nodes.size());
        std::iota(index.begin(), index.end(), 0);
    }//_if
    
    //Convert to Nifti
    int ctr = 0;
    QString path;
    this->BusyCursorOn();
    mitk::ProgressBar::GetInstance()->AddStepsToDo(index.size());
    foreach (int idx, index) {
        mitk::BaseData::Pointer data = nodes.at(idx)->GetData();
        if (data) {
            //Test if this data item is an image
            mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
            if (image) {
                path = directory + mitk::IOUtil::GetDirectorySeparator() + "dcm-" + QString::number(ctr++) + ".nii";
                mitk::IOUtil::Save(image, path.toStdString());
                this->GetDataStorage()->Remove(nodes.at(idx));
            } else
                return;
        } else
            return;
        mitk::ProgressBar::GetInstance()->Progress();
    }//for
    nodes.clear();
    this->BusyCursorOff();
    
    //Load first item
    ctr = 0;
    path = directory + mitk::IOUtil::GetDirectorySeparator() + "dcm-" + QString::number(ctr) + ".nii";
    mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
}

/**
 * @brief QmitkAbstractView::CropinIMGS
 */
void MmcwView::CropinIMGS() {
    
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

        //Check for temporal resolution
        bool ok = true;
        if (timePoints == 0)
            timePoints = QInputDialog::getInt(NULL, tr("Time Points"), tr("Resolution:"), 10, 1, 100, 1, &ok);
        if (!ok) {
            QMessageBox::warning(NULL, "Attention", "Enter a correct value for the temporal resolution!");
            timePoints = 0;
            return;
        }//_if
        
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
        
        //Cut rest of images
        int reply = QMessageBox::question(
                    NULL, "Question", "Would you like to automate cropping of other images in the cycle?",
                    QMessageBox::Yes, QMessageBox::No);
        if (reply == QMessageBox::Yes) {
            
            this->BusyCursorOn();
            mitk::ProgressBar::GetInstance()->AddStepsToDo((timePoints-1)*2);
            for (int i=1; i<timePoints; i++) {
                
                mitk::Image::Pointer inputImage;
                path = directory + mitk::IOUtil::GetDirectorySeparator() + "dcm-" + QString::number(i) + ".nii";
                try {
                    inputImage = dynamic_cast<mitk::Image*>(
                                mitk::IOUtil::Load(path.toStdString()).front().GetPointer());
                } catch(const std::exception& e) {
                    continue;
                }//_try

                //Setup cropper
                CemrgImageUtils::SetImageToCut(inputImage);
                outputImage = CemrgImageUtils::CropImage();
                mitk::IOUtil::Save(outputImage, path.toStdString());

            }//_for
            this->BusyCursorOff();
        }//_if
        
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

/**
 * @brief MmcwView::ResampIMGS
 */
void MmcwView::ResampIMGS() {
    
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

    //Check for temporal resolution
    bool ok = true;
    if (timePoints == 0)
        timePoints = QInputDialog::getInt(NULL, tr("Time Points"), tr("Resolution:"), 10, 1, 100, 1, &ok);
    if (!ok) {
        QMessageBox::warning(NULL, "Attention", "Enter a correct value for the temporal resolution!");
        timePoints = 0;
        return;
    }//_if
    
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
                
                //Downsample rest of images
                int reply = QMessageBox::question(
                            NULL, "Question", "Would you like to automate downsampling of other images in the cycle?",
                            QMessageBox::Yes, QMessageBox::No);
                
                if (reply == QMessageBox::Yes) {
                    
                    this->BusyCursorOn();
                    mitk::ProgressBar::GetInstance()->AddStepsToDo(9);
                    for (int i=1; i<timePoints; i++) {
                        
                        mitk::Image::Pointer inputImage;
                        path = directory + mitk::IOUtil::GetDirectorySeparator() + "dcm-" + QString::number(i) + ".nii";
                        try {
                            inputImage = dynamic_cast<mitk::Image*>(
                                        mitk::IOUtil::Load(path.toStdString()).front().GetPointer());
                        } catch(const std::exception& e) {
                            continue;
                        }//_try
                        
                        //Setup sampler
                        outputImage = CemrgImageUtils::Downsample(inputImage, factor);
                        mitk::IOUtil::Save(outputImage, path.toStdString());

                    }//_for
                    this->BusyCursorOff();
                }//_if
            }//_if
        } else
            return;
    } else
        return;
}

/**
 * @brief MmcwView::SegmentIMGS
 */
void MmcwView::SegmentIMGS() {
    
    int reply = QMessageBox::question(
                NULL, "Question", "Do you have a segmentation to load?", QMessageBox::Yes, QMessageBox::No);
    if (reply == QMessageBox::Yes) {
        
        QString path = QFileDialog::getOpenFileName(
                    NULL, "Open Segmentation file", mitk::IOUtil::GetProgramPath().c_str(), QmitkIOUtil::GetFileOpenFilterString());
        if (path.isEmpty())
            return;

        int cropReply = QMessageBox::question(
                    NULL, "Question", "Do you want to crop your segmentation?", QMessageBox::Yes, QMessageBox::No);
        if (cropReply == QMessageBox::Yes) {

            this->BusyCursorOn();
            mitk::Image::Pointer segImage = mitk::IOUtil::LoadImage(path.toStdString());
            CemrgImageUtils::SetImageToCut(segImage);
            mitk::Image::Pointer outImage = CemrgImageUtils::CropImage();
            if (outImage.IsNull()) {
                QMessageBox::critical(NULL, "Attention", "There is no previously used cutter to use now!");
                this->BusyCursorOff();
                return;
            }//_cropper
            mitk::DataNode::Pointer node = mitk::DataNode::New();
            node->SetData(outImage);
            node->SetName("CroppedSegmentation");
            this->GetDataStorage()->Add(node);
            this->BusyCursorOff();
            return;

        }//_crop

        mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
        mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
        
    } else {

        //Show the plugin
        this->GetSite()->GetPage()->ShowView("org.mitk.views.segmentation");

    }//_loadSegmentation
}

/**
 * @brief MmcwView::CreateSurf
 */
void MmcwView::CreateSurf() {
    
    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please select a segmentation from the Data Manager to create a surface!");
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
    mitk::DataNode::Pointer segNode = nodes.at(0);
    mitk::BaseData::Pointer data = segNode->GetData();
    if (data) {
        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image) {
            path = directory + mitk::IOUtil::GetDirectorySeparator() + "segmentation.nii";
            mitk::IOUtil::Save(image, path.toStdString());
            this->GetDataStorage()->Remove(segNode);
        } else
            return;
    } else
        return;
    
    //Ask for user input to set the parameters
    QDialog* inputs = new QDialog(0,0);
    m_UIMeshing.setupUi(inputs);
    connect(m_UIMeshing.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_UIMeshing.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    int dialogCode = inputs->exec();
    
    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {
        
        bool ok1, ok2, ok3, ok4;
        int iter = m_UIMeshing.lineEdit_1->text().toInt(&ok1);
        float th = m_UIMeshing.lineEdit_2->text().toFloat(&ok2);
        int blur = m_UIMeshing.lineEdit_3->text().toInt(&ok3);
        int smth = m_UIMeshing.lineEdit_4->text().toInt(&ok4);
        
        //Set default values
        if (!ok1 || !ok2 || !ok3 || !ok4)
            QMessageBox::warning(NULL, "Attention", "Reverting to default parameters!");
        if(!ok1) iter = 1;
        if(!ok2) th   = 0.5;
        if(!ok3) blur = 3;
        if(!ok4) smth = 40;
        //_if
        
        this->BusyCursorOn();
        mitk::ProgressBar::GetInstance()->AddStepsToDo(4);
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        QString output = cmd->ExecuteSurf(directory, path, iter, th, blur, smth);
        QMessageBox::information(NULL, "Attention", "Command Line Operations Finished!");
        this->BusyCursorOff();
        
        //Add the mesh to storage
        AddToStorage("Segmented Mesh", ReadVTKMesh(output.toStdString()));
        mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
        inputs->deleteLater();
        
    } else if(dialogCode == QDialog::Rejected) {
        inputs->close();
        inputs->deleteLater();
    }//_if
}

void MmcwView::TrackingIMGS() {
    
    //Toggle visibility of buttons
    if (m_Controls.button_5_1->isVisible()) {
        m_Controls.button_5_1->setVisible(false);
        m_Controls.button_5_2->setVisible(false);
        m_Controls.button_5_3->setVisible(false);
    } else {
        m_Controls.button_5_1->setVisible(true);
        m_Controls.button_5_2->setVisible(true);
        m_Controls.button_5_3->setVisible(true);
    }
}

void MmcwView::BrowseT(const QString& buttDir) {
    
    QString time, para = "";
    QString buttID = buttDir.left(1);
    QString direct = buttDir.right(buttDir.size()-1);
    
    //Load target, time and parameter files
    switch (buttID.toInt()) {
    case 1:
        time = QFileDialog::getOpenFileName(
                    NULL, "Open text file containing time points of source images",
                    direct, QmitkIOUtil::GetFileOpenFilterString());
        m_UITracking.lineEdit_1->setText(time);
        break;
    case 2:
        para = QFileDialog::getOpenFileName(
                    NULL, "Open text file containing parameters",
                    direct, QmitkIOUtil::GetFileOpenFilterString());
        m_UITracking.lineEdit_2->setText(para);
        break;
    }//_switch
}

void MmcwView::BrowseA(const QString& buttDir) {
    
    QString input, dofin = "";
    QString buttID = buttDir.left(1);
    QString direct = buttDir.right(buttDir.size()-1);
    
    //Load input mesh, dofin file
    switch (buttID.toInt()) {
    case 1:
        input = QFileDialog::getOpenFileName(
                    NULL, "Open the input mesh",
                    direct, QmitkIOUtil::GetFileOpenFilterString());
        m_UIApplying.lineEdit_1->setText(input);
        break;
    case 2:
        dofin = QFileDialog::getOpenFileName(
                    NULL, "Open the transformation file",
                    direct, QmitkIOUtil::GetFileOpenFilterString());
        m_UIApplying.lineEdit_3->setText(dofin);
        break;
    }//_switch
}

void MmcwView::Tracking() {
    
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

    //Check for temporal resolution
    bool ok = true;
    if (timePoints == 0)
        timePoints = QInputDialog::getInt(NULL, tr("Time Points"), tr("Resolution:"), 10, 1, 100, 1, &ok);
    if (!ok) {
        QMessageBox::warning(NULL, "Attention", "Enter a correct value for the temporal resolution!");
        timePoints = 0;
        return;
    }//_if
    
    //Ask for user input to set the parameters
    QDialog* inputs = new QDialog(0,0);
    QSignalMapper* signalMapper = new QSignalMapper(this);
    
    m_UITracking.setupUi(inputs);
    connect(m_UITracking.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_UITracking.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    connect(m_UITracking.pushButton_1, SIGNAL(clicked()), signalMapper, SLOT(map()));
    connect(m_UITracking.pushButton_2, SIGNAL(clicked()), signalMapper, SLOT(map()));
    signalMapper->setMapping(m_UITracking.pushButton_1, "1"+directory);
    signalMapper->setMapping(m_UITracking.pushButton_2, "2"+directory);
    connect(signalMapper, SIGNAL(mapped(QString)), this, SLOT(BrowseT(const QString&)));
    
    int dialogCode = inputs->exec();
    
    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {
        
        QString time = m_UITracking.lineEdit_1->text();
        QString para = m_UITracking.lineEdit_2->text();
        
        //Checking input files
        if (time.isEmpty() || para.isEmpty())
            QMessageBox::warning(NULL, "Attention", "Reverting to default time or parameter file!");
        
        if (time.isEmpty()) {
            ofstream file;
            //Absolute path
            QString aPath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
            aPath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                    mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                    mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
            file.open(aPath.toStdString() + mitk::IOUtil::GetDirectorySeparator() + "imgTimes.lst");
            file << directory << mitk::IOUtil::GetDirectorySeparator() << "dcm- .nii" << endl;
            for (int i=0; i<timePoints; i++)
                file << i << " " << i*10 << endl;
            file.close();
            time = aPath + mitk::IOUtil::GetDirectorySeparator() + "imgTimes.lst";
        }//_if
        
        //Commandline execution
        this->BusyCursorOn();
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->ExecuteTracking(directory, time, para);
        QMessageBox::information(NULL, "Attention", "Command Line Operations Finished!");
        this->BusyCursorOff();
        signalMapper->deleteLater();
        inputs->deleteLater();
        
    } else if(dialogCode == QDialog::Rejected) {
        inputs->close();
        signalMapper->deleteLater();
        inputs->deleteLater();
    }//_if
}

void MmcwView::Applying() {
    
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

    //Check for temporal resolution
    bool ok = true;
    if (timePoints == 0)
        timePoints = QInputDialog::getInt(NULL, tr("Time Points"), tr("Resolution:"), 10, 1, 100, 1, &ok);
    if (!ok) {
        QMessageBox::warning(NULL, "Attention", "Enter a correct value for the temporal resolution!");
        timePoints = 0;
        return;
    }//_if
    
    //Ask for user input to set the parameters
    QDialog* inputs = new QDialog(0,0);
    QSignalMapper* signalMapper = new QSignalMapper(this);
    
    m_UIApplying.setupUi(inputs);
    connect(m_UIApplying.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_UIApplying.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    connect(m_UIApplying.pushButton_1, SIGNAL(clicked()), signalMapper, SLOT(map()));
    connect(m_UIApplying.pushButton_2, SIGNAL(clicked()), signalMapper, SLOT(map()));
    signalMapper->setMapping(m_UIApplying.pushButton_1, "1"+directory);
    signalMapper->setMapping(m_UIApplying.pushButton_2, "2"+directory);
    connect(signalMapper, SIGNAL(mapped(QString)), this, SLOT(BrowseA(const QString&)));
    m_UIApplying.lineEdit_4->setPlaceholderText("Enter Number of Frames (default = " + QString::number(timePoints) + ")");
    
    int dialogCode = inputs->exec();
    
    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {
        
        //Load input mesh, dofin file
        QString input = m_UIApplying.lineEdit_1->text();
        QString dofin = m_UIApplying.lineEdit_3->text();
        
        //Load initial time, number of frames, smoothness
        bool ok1, ok2;
        int smoothness;
        double iniTime = m_UIApplying.lineEdit_2->text().toDouble(&ok1);
        int frames = m_UIApplying.lineEdit_4->text().toInt(&ok2);
        if (m_UIApplying.radioButton_1->isChecked()) smoothness = 1;
        else if (m_UIApplying.radioButton_2->isChecked()) smoothness = 2;
        else if (m_UIApplying.radioButton_3->isChecked()) smoothness = 5;
        
        //Checking input values
        if (input.isEmpty() || dofin.isEmpty()) {
            QMessageBox::warning(NULL, "Attention", "Please select the input files to apply tracking!");
            signalMapper->deleteLater();
            inputs->deleteLater();
            return;
        }
        if (!ok1 || !ok2)
            QMessageBox::warning(NULL, "Attention", "Reverting to default parameters!");
        if(!ok1) iniTime = 0;
        if(!ok2) frames = timePoints;
        //_if
        
        //Commandline execution
        this->BusyCursorOn();
        mitk::ProgressBar::GetInstance()->AddStepsToDo(frames*smoothness);
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->ExecuteApplying(directory, input, iniTime, dofin, frames, smoothness);
        QMessageBox::information(NULL, "Attention", "Command Line Operations Finished!");
        MmcwViewPlot::SetNoFrames(frames, smoothness);
        this->BusyCursorOff();
        signalMapper->deleteLater();
        inputs->deleteLater();
        
    } else if (dialogCode == QDialog::Rejected) {
        inputs->close();
        signalMapper->deleteLater();
        inputs->deleteLater();
    }//_if
}

void MmcwView::Demoings() {
    
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

    //Check for temporal resolution
    bool ok = true;
    if (timePoints == 0)
        timePoints = QInputDialog::getInt(NULL, tr("Time Points"), tr("Resolution:"), 10, 1, 100, 1, &ok);
    if (!ok) {
        QMessageBox::warning(NULL, "Attention", "Enter a correct value for the temporal resolution!");
        timePoints = 0;
        return;
    }//_if
    
    //Read all images and meshes from the project directory
    QString path;
    mitk::Image::Pointer img3D;
    mitk::Image::Pointer img4D = mitk::Image::New();
    mitk::Surface::Pointer sur3D;
    mitk::Surface::Pointer sur4D = mitk::Surface::New();
    
    this->BusyCursorOn();
    mitk::ProgressBar::GetInstance()->AddStepsToDo(timePoints);
    
    for (int tS=0; tS<timePoints; tS++) {
        
        //Image
        path = directory + mitk::IOUtil::GetDirectorySeparator() + "dcm-" + QString::number(tS) + ".nii";
        img3D = mitk::IOUtil::LoadImage(path.toStdString());
        //Initialise
        if (tS==0) {
            mitk::ImageDescriptor::Pointer dsc = img3D->GetImageDescriptor();
            img4D->Initialize(dsc->GetChannelDescriptor(0).GetPixelType(), *img3D->GetGeometry(), dsc->GetNumberOfChannels(), timePoints);
        }//_if
        img4D->SetVolume(mitk::ImageReadAccessor(img3D).GetData(), tS);
        
        //Mesh
        path = directory + mitk::IOUtil::GetDirectorySeparator() + "transformed-" + QString::number(tS) + ".vtk";
        sur3D = this->ReadVTKMesh(path.toStdString());
        sur4D->SetVtkPolyData(sur3D->GetVtkPolyData(), tS);
        
        mitk::ProgressBar::GetInstance()->Progress();
    }
    
    for(int i=0; i<timePoints; i++)
        sur4D->GetGeometry(i)->SetBounds(sur3D->GetGeometry()->GetBounds());
    this->BusyCursorOff();
    
    AddToStorage("4DImage", img4D);
    AddToStorage("4DMesh", sur4D);
    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
    
    //Show the plugin
    this->GetSite()->GetPage()->ShowView("org.mitk.views.imagenavigator");
}

void MmcwView::LandmarkSelection() {
    
    //Show the plugin
    QMessageBox::information(
                NULL, "Attention",
                "Please select 6 points in order:\n\n1 on the Apex\n3 on the Mitral Valve surface\n2 on the Right Ventricle cusps\n");
    this->GetSite()->GetPage()->ShowView("org.mitk.views.pointsetinteraction");
}

void MmcwView::Plotting() {
    
    //Show the plugin
    MmcwViewPlot::SetDirectory(directory);
    this->GetSite()->GetPage()->ResetPerspective();
    this->GetSite()->GetPage()->ShowView("my.cemrgproject.views.mmcwviewplot");
}

void MmcwView::Reset() {
    
    try {
        
        this->GetSite()->GetPage()->ResetPerspective();
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
    timePoints = 0;
    directory.clear();
    m_Controls.button_2_2->setText("Crop Images");
}

void MmcwView::AddToStorage(std::string nodeName, mitk::BaseData* data) {
    
    if (!data) return;
    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(data);
    node->SetName(nodeName);
    this->GetDataStorage()->Add(node);
}

mitk::Surface::Pointer MmcwView::ReadVTKMesh(std::string meshPath) {
    
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
}
