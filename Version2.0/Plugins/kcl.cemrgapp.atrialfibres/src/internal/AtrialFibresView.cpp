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
#include <berryFileEditorInput.h>

//Qmitk
#include <mitkImage.h>
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
#include <mitkManualSegmentationToSurfaceFilter.h>
#include "kcl_cemrgapp_atrialfibres_Activator.h"
#include "AtrialFibresView.h"
#include "AtrialFibresClipperView.h"

//Micro services
#include <usModuleRegistry.h>
#ifdef _WIN32
// _WIN32 = we're in windows
#include <winsock2.h>
#else
// or linux/mac
#include <arpa/inet.h>
#endif

//VTK
#include <vtkFieldData.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkTimerLog.h>
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

//Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>

//CemrgAppModule
#include <CemrgCommonUtils.h>
#include <CemrgCommandLine.h>
#include <CemrgMeasure.h>

const std::string AtrialFibresView::VIEW_ID = "org.mitk.views.atrialfibresview"; // org.mitk.views.atrialfibresview

void AtrialFibresView::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, SIGNAL(clicked()), this, SLOT(LoadDICOM()));
    connect(m_Controls.button_2, SIGNAL(clicked()), this, SLOT(ProcessIMGS()));
    // Segmentation to Labelled Mesh pipeline
    connect(m_Controls.button_3_segmentation, SIGNAL(clicked()), this, SLOT(SegmentIMGS()));
    connect(m_Controls.button_4_labelmesh, SIGNAL(clicked()), this, SLOT(CreateLabelledMesh()));
    connect(m_Controls.button_5_pvclipper, SIGNAL(clicked()), this, SLOT(SelectLandmarks())); // pv clipper
    connect(m_Controls.button_6_mvclipper, SIGNAL(clicked()), this, SLOT(ClipperMV()));
    // Labelled Mesh to UAC
    connect(m_Controls.button_x_meshtools, SIGNAL(clicked()), this, SLOT(MeshingOptions()));
    connect(m_Controls.button_y_calculateUac, SIGNAL(clicked()), this, SLOT(UacCalculation()));
    connect(m_Controls.button_z_refineUac, SIGNAL(clicked()), this, SLOT(UacMeshRefinement()));

    connect(m_Controls.button_r, SIGNAL(clicked()), this, SLOT(Reset()));

    //Set visibility of buttons
    m_Controls.button_2_1->setVisible(false);
    connect(m_Controls.button_2_1, SIGNAL(clicked()), this, SLOT(ConvertNII()));

    // Set default variables
    tagName = "tag-segmentation";
    askedAboutAutoPipeline = false;
}

void AtrialFibresView::SetFocus() {
    m_Controls.button_1->setFocus();
}

void AtrialFibresView::OnSelectionChanged(
        berry::IWorkbenchPart::Pointer /*source*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

void AtrialFibresView::LoadDICOM() {
    int reply1 = QMessageBox::No;
#if defined(__APPLE__)
    MITK_INFO << "Ask user about alternative DICOM reader";
    reply1 = QMessageBox::question(NULL, "Question", "Use alternative DICOM reader?",
                        QMessageBox::Yes, QMessageBox::No);
#endif

    if (reply1 == QMessageBox::Yes) {

        QString dicomFolder = QFileDialog::getExistingDirectory(NULL, "Open folder with DICOMs.", mitk::IOUtil::GetProgramPath().c_str(), QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        QString tmpNiftiFolder = cmd->DockerDicom2Nifti(dicomFolder);

        if (tmpNiftiFolder.compare("ERROR_IN_PROCESSING") != 0) {

            // add results in NIIs folder to Data Manager
            MITK_INFO << ("Conversion succesful. Intermediate NII folder: " + tmpNiftiFolder).toStdString();
            QMessageBox::information(NULL, "Information", "Conversion successful, press the Process Images button to continue.");
            QDir niftiFolder(tmpNiftiFolder);
            QStringList niftiFiles = niftiFolder.entryList();

            if (niftiFiles.size()>0) {

                QString thisFile, path;
                for(int ix=0; ix<niftiFiles.size(); ix++) {

                    // load here files
                    thisFile = niftiFiles.at(ix);
                    if (thisFile.contains(".nii", Qt::CaseSensitive)) {
                        if (thisFile.contains("lge", Qt::CaseInsensitive) ||  thisFile.contains("mra", Qt::CaseInsensitive)) {

                            path = niftiFolder.absolutePath() + mitk::IOUtil::GetDirectorySeparator() + thisFile;
                            mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(path.toStdString());
                            std::string key = "dicom.series.SeriesDescription";
                            mitk::DataStorage::SetOfObjects::Pointer set = mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
                            set->Begin().Value()->GetData()->GetPropertyList()->SetStringProperty(key.c_str(), thisFile.left(thisFile.length()-4).toStdString().c_str());

                        }//_if
                    }//_if
                }//_for

            } else {

                MITK_WARN << "Problem with conversion.";
                QMessageBox::warning(NULL, "Attention", "Problem with alternative conversion. Try MITK Dicom editor?");
                return;

            }//_if
        }//_if

    } else {

        MITK_INFO << "Using MITK DICOM editor";
        QString editor_id = "org.mitk.editors.dicomeditor";
        berry::IEditorInput::Pointer input(new berry::FileEditorInput(QString()));
        this->GetSite()->GetPage()->OpenEditor(input, editor_id);

    }//_if
}

void AtrialFibresView::ProcessIMGS() {
    //Toggle visibility of buttons
    if (m_Controls.button_2_1->isVisible()) {
        m_Controls.button_2_1->setVisible(false);
    } else {
        m_Controls.button_2_1->setVisible(true);
    }
}

void AtrialFibresView::ConvertNII() {
    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() < 1) {
        MITK_WARN << "load and select images from the Data Manager before starting this step!";
        QMessageBox::warning(NULL, "Attention",
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
    bool successfulNitfi, resampleImage, reorientToRAI;
    resampleImage = false;
    reorientToRAI = true;

    this->BusyCursorOn();
    mitk::ProgressBar::GetInstance()->AddStepsToDo(nodes.size());
    foreach (mitk::DataNode::Pointer node, nodes) {
        path = directory + mitk::IOUtil::GetDirectorySeparator() + "dcm-" + QString::number(ctr++) + ".nii";
        successfulNitfi = CemrgCommonUtils::ConvertToNifti(node->GetData(), path, resampleImage, reorientToRAI);
        if (successfulNitfi) {
            this->GetDataStorage()->Remove(node);
        } else {
            mitk::ProgressBar::GetInstance()->Progress(nodes.size());
            return;
        }
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


// Segmentation to Labelled Mesh pipeline
void AtrialFibresView::SegmentIMGS() {
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    int replyAuto = QMessageBox::question(NULL, "Question", "Do you want an automatic segmentation?",
            QMessageBox::Yes, QMessageBox::No);

    if(replyAuto == QMessageBox::Yes){
        //Check for selection of image
        QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
        if (nodes.size() != 1) {
            QMessageBox::warning(NULL, "Attention",
                    "Please select the CEMRA images from the Data Manager to segment!");
            return;
        }//_if

        //Test if this data item is an image
        QString mraPath;
        mitk::BaseData::Pointer data = nodes.at(0)->GetData();
        if (data) {
            mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
            if (image) {

                this->BusyCursorOn();
                mitk::ProgressBar::GetInstance()->AddStepsToDo(2);

                //CNN prediction
                mraPath = directory + mitk::IOUtil::GetDirectorySeparator() + "test.nii";
                mitk::IOUtil::Save(image, mraPath.toStdString());
                std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
                cmd->SetUseDockerContainers(true);
                cnnPath = cmd->DockerCemrgNetPrediction(mraPath);
                mitk::ProgressBar::GetInstance()->Progress();

                //Clean prediction
                ImageType::Pointer orgSegImage = ImageType::New();
                mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>(cnnPath.toStdString()), orgSegImage);
                ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
                connected->SetInput(orgSegImage);
                connected->Update();
                LabelShapeKeepNObjImgFilterType::Pointer lblShpKpNObjImgFltr = LabelShapeKeepNObjImgFilterType::New();
                lblShpKpNObjImgFltr->SetInput(connected->GetOutput());
                lblShpKpNObjImgFltr->SetBackgroundValue(0);
                lblShpKpNObjImgFltr->SetNumberOfObjects(1);
                lblShpKpNObjImgFltr->SetAttribute(LabelShapeKeepNObjImgFilterType::LabelObjectType::NUMBER_OF_PIXELS);
                lblShpKpNObjImgFltr->Update();
                mitk::Image::Pointer segImage = mitk::Image::New();
                mitk::CastToMitkImage(lblShpKpNObjImgFltr->GetOutput(), segImage);
                cnnPath = directory + mitk::IOUtil::GetDirectorySeparator() + "LA.nii";
                mitk::IOUtil::Save(segImage, cnnPath.toStdString());
                mitk::IOUtil::Load(cnnPath.toStdString(), *this->GetDataStorage());
                remove(mraPath.toStdString().c_str());
                fileName = "LA.nii";
                QMessageBox::information(NULL, "Attention", "Command Line Operations Finished!");
                mitk::ProgressBar::GetInstance()->Progress();
                this->BusyCursorOff();

            } else
                QMessageBox::warning(NULL, "Attention", "Please select a CEMRA to segment!");
        }//_if_data
    } else {
        int replyLoad =  QMessageBox::question(NULL, "Question", "Do you have a segmentation to load?",
                QMessageBox::Yes, QMessageBox::No);

        if (replyLoad == QMessageBox::Yes) {
            QString path = QFileDialog::getOpenFileName(NULL, "Open Segmentation file",
                                                        directory.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());
            if (path.isEmpty()) return;
            mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
            mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());

            //Restore image name
            QFileInfo fullPathInfo(path);
            fileName = fullPathInfo.fileName();
        } else {
            //Show segmentation plugin
            this->GetSite()->GetPage()->ShowView("org.mitk.views.segmentation");
        } //_if_q_loadSegmentation
    }//_if_q_automatic
}

void AtrialFibresView::CreateLabelledMesh(){
    if(cnnPath.isEmpty()){
        if (!RequestProjectDirectoryFromUser()) return; // checks if the directory has been set

        QString cnnPath = QFileDialog::getOpenFileName(NULL, "Open Automatic Segmentation file",
                        directory.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());

        QFileInfo fi(cnnPath);
        if(directory.compare(fi.absolutePath())==0){
            MITK_INFO << ("[ATTENTION] Changing working directory to: " + fi.absolutePath());
            directory = fi.absolutePath();
        }
    }

    atrium = std::unique_ptr<CemrgAtrialTools>(new CemrgAtrialTools());;
    atrium->SetDebugModeOn();
    atrium->SetWorkingDirectory(directory);
    ImageType::Pointer tagAtrium = atrium->CleanAutomaticSegmentation(directory);
    ImageType::Pointer tagAtriumAuto = atrium->AssignAutomaticLabels(tagAtrium, directory);

    double th=0.5;
    double bl=1;
    double smth=3;
    double ds=0.5;
    atrium->GetSurfaceWithTags(tagAtriumAuto, directory, tagName+".vtk", th, bl, smth, ds);
}

// pv clipper
void AtrialFibresView::SelectLandmarks(){
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    //Show the plugin
    this->GetSite()->GetPage()->ResetPerspective();
    AtrialFibresClipperView::SetDirectoryFile(directory, tagName+".vtk", automaticPipeline);
    this->GetSite()->GetPage()->ShowView("org.mitk.views.atrialfibresclipperview");
}

// mv clipper
void AtrialFibresView::ClipperMV(){
    MITK_INFO << "Clipping Mitral Valve";
    if(automaticPipeline){
        bool oldDebug = atrium->Debugging();
        atrium->SetDebugModeOn();
        atrium->ClipMitralValveAuto(directory, "prodMVI.nii", tagName+".vtk");
        atrium->SetDebugMode(oldDebug);
    } else{

    }
}

// Labelled Mesh to UAC
void AtrialFibresView::MeshingOptions(){

}

void AtrialFibresView::UacCalculation(){

}

void AtrialFibresView::UacMeshRefinement(){

}

void AtrialFibresView::Reset() {

    try {

        ctkPluginContext* context = mitk::kcl_cemrgapp_atrialfibres_Activator::getContext();
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
    directory.clear();
    this->GetSite()->GetPage()->ResetPerspective();
}

// helper functions
bool AtrialFibresView::RequestProjectDirectoryFromUser() {

    bool succesfulAssignment = true;

    //Ask the user for a dir to store data
    if (directory.isEmpty()) {

        MITK_INFO << "Directory is empty. Requesting user for directory.";
        directory = QFileDialog::getExistingDirectory(
                    NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
                    QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);
        MITK_INFO << ("Directory selected:" + directory).toStdString();
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            MITK_WARN << "Please select a project directory with no spaces in the path!";
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            succesfulAssignment = false;
        }//_if

    } else {
        MITK_INFO << ("Project directory already set: " + directory).toStdString();
    }//_if

    if(!askedAboutAutoPipeline){
        int reply1 = QMessageBox::question(NULL, "Question",
        "Automatic Pipeline?",
        QMessageBox::Yes, QMessageBox::No);
        SetAutomaticPipeline(reply1==QMessageBox::Yes);
        askedAboutAutoPipeline = true;
    }

    return succesfulAssignment;
}
