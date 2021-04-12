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
#include "AtrialFibresLandmarksView.h"

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
#include <vtkDecimatePro.h>

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
#include <QDirIterator>

//CemrgAppModule
#include <CemrgCommonUtils.h>
#include <CemrgCommandLine.h>
#include <CemrgMeasure.h>
#include <CemrgScar3D.h>

const std::string AtrialFibresView::VIEW_ID = "org.mitk.views.atrialfibresview"; // org.mitk.views.atrialfibresview

void AtrialFibresView::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, SIGNAL(clicked()), this, SLOT(LoadDICOM()));
    connect(m_Controls.button_2, SIGNAL(clicked()), this, SLOT(ProcessIMGS()));
    // Segmentation to Labelled Mesh pipeline
    connect(m_Controls.button_3_imanalysis, SIGNAL(clicked()), this, SLOT(AnalysisChoice()));
    connect(m_Controls.button_auto4_meshpreproc, SIGNAL(clicked()), this, SLOT(MeshPreprocessing()));
    connect(m_Controls.button_man4_segmentation, SIGNAL(clicked()), this, SLOT(SegmentIMGS()));
    connect(m_Controls.button_auto5_clipPV, SIGNAL(clicked()), this, SLOT(ClipperPV()));
    connect(m_Controls.button_man5_idPV, SIGNAL(clicked()), this, SLOT(IdentifyPV())); // pv clipper
    connect(m_Controls.button_man6_labelmesh, SIGNAL(clicked()), this, SLOT(CreateLabelledMesh()));
    connect(m_Controls.button_man7_clipMV, SIGNAL(clicked()), this, SLOT(ClipperMV()));
    connect(m_Controls.button_man8_clipPV, SIGNAL(clicked()), this, SLOT(ClipperPV()));
    connect(m_Controls.button_0_landmarks, SIGNAL(clicked()), this, SLOT(SelectLandmarks()));

    // Labelled Mesh to UAC
    connect(m_Controls.button_x_meshtools, SIGNAL(clicked()), this, SLOT(MeshingOptions()));
    connect(m_Controls.button_y_vtk2carp, SIGNAL(clicked()), this, SLOT(ConvertFormat()));
    connect(m_Controls.button_0_calculateUac, SIGNAL(clicked()), this, SLOT(UacCalculation()));
    connect(m_Controls.button_0_refineUac, SIGNAL(clicked()), this, SLOT(UacMeshRefinement()));

    // Scar analysis
    connect(m_Controls.button_z_scar, SIGNAL(clicked()), this, SLOT(ScarProjection()));

    connect(m_Controls.button_r, SIGNAL(clicked()), this, SLOT(Reset()));

    //Set visibility of buttons
    m_Controls.button_2_1->setVisible(false);
    connect(m_Controls.button_2_1, SIGNAL(clicked()), this, SLOT(ConvertNII()));

    // Set default variables
    tagName = "labelled";
    refinedSuffix = "-refined";
    analysisOnLge = false;
    uiRemesh_vtk2carp = true;
    uiRemesh_extractParts = false;
    atrium = std::unique_ptr<CemrgAtrialTools>(new CemrgAtrialTools());
    atrium->SetDebugModeOff();

    SetManualModeButtonsOff();
    SetAutomaticModeButtonsOff();
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
    reply1 = Ask("Question", "Use alternative DICOM reader?");
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

    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    //Order dicoms based on their type
    std::vector<int> indexNodes;
    std::vector<std::string> seriesDscrps;
    foreach (mitk::DataNode::Pointer node, nodes) {

        std::string seriesDescription;
        node->GetData()->GetPropertyList()->GetStringProperty("dicom.series.SeriesDescription", seriesDescription);

        if (seriesDescription.find("LGE")      != seriesDescription.npos) indexNodes.push_back(0);
        else if (seriesDescription.find("MRA") != seriesDescription.npos) indexNodes.push_back(1);

        //Trim whitespaces
        seriesDescription = QString::fromStdString(seriesDescription).replace(")","").toStdString();
        seriesDescription = QString::fromStdString(seriesDescription).replace("(","").toStdString();
        seriesDescription = QString::fromStdString(seriesDescription).simplified().replace(" ","").toStdString();
        seriesDscrps.push_back(seriesDescription);
    }//_for

    //Sort indexes based on comparing values
    std::vector<int> index(indexNodes.size());
    std::iota(index.begin(), index.end(), 0);
    std::sort(index.begin(), index.end(), [&](int i1, int i2) {return indexNodes[i1]<indexNodes[i2];});

    //Warning for cases when type is not found
    size_t length1 = nodes.size();
    size_t length2 = indexNodes.size();
    bool test = std::adjacent_find(indexNodes.begin(), indexNodes.end(), std::not_equal_to<int>()) == indexNodes.end();
    if (length1 != length2 || test) {
        QMessageBox::warning(NULL, "Attention",
                             "Cannot find the type of images automatically. Revert to user order and selections in the data manager: LGE at the top, then CEMRA at the bottom!");
        index.resize(nodes.size());
        std::iota(index.begin(), index.end(), 0);
    }//_if

    //Convert to Nifti
    int ctr = 0;
    QString prodPath, type;
    bool successfulNitfi, resampleImage, reorientToRAI;
    resampleImage = true;
    reorientToRAI = true;

    this->BusyCursorOn();
    mitk::ProgressBar::GetInstance()->AddStepsToDo(index.size());
    foreach (int idx, index) {
        type = (ctr==0) ? "LGE":"MRA";
        prodPath = directory + mitk::IOUtil::GetDirectorySeparator() + "dcm-" + type + "-" + seriesDscrps.at(idx).c_str() + ".nii";
        successfulNitfi = CemrgCommonUtils::ConvertToNifti(nodes.at(idx)->GetData(), prodPath, resampleImage, reorientToRAI);
        if (successfulNitfi) {
            this->GetDataStorage()->Remove(nodes.at(idx));
            std::string key = "dicom.series.SeriesDescription";
            mitk::DataStorage::SetOfObjects::Pointer set = mitk::IOUtil::Load(prodPath.toStdString(), *this->GetDataStorage());
            set->Begin().Value()->GetData()->GetPropertyList()->SetStringProperty(key.c_str(), seriesDscrps.at(idx).c_str());
            ctr++;
        } else {
            mitk::ProgressBar::GetInstance()->Progress(index.size());
            return;
        }//_if
        mitk::ProgressBar::GetInstance()->Progress();
    }//for
    nodes.clear();
    this->BusyCursorOff();

    MITK_INFO << "Loading all items";
    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
}

void AtrialFibresView::AnalysisChoice(){
    int reply1 = Ask("Question", "Use the semi-automatic pipeline?");

    SetAutomaticPipeline(reply1==QMessageBox::Yes);

    if(reply1==QMessageBox::Yes){
        MITK_INFO<<"Automatic pipeline";
        if (!RequestProjectDirectoryFromUser()) return; // checks if the directory has been set
        SetAutomaticModeButtonsOn();

        int reply = Ask("Check for previous output!","Do you want to skip steps 1 to 3?");
        if(reply==QMessageBox::No){
            MITK_INFO << "Performing Automatic analysis on CemrgNet prediction.";
            AutomaticAnalysis();
        } else if(reply==QMessageBox::Yes){
            MITK_INFO << "Skipping Neural network prediction. Checking for vtk file.";
            if(!QFile::exists(directory+mitk::IOUtil::GetDirectorySeparator()+tagName+".vtk")){
                std::string msg = "Labelled mesh not found. Perform Step 3 again.";
                QMessageBox::warning(NULL, "File not found", msg.c_str());
                MITK_INFO << msg;
                return;
            }
        }
    } else{
        MITK_INFO << "Setting up manual analysis.";
        SetManualModeButtonsOn();
    }
}

// Automatic pipeline
void AtrialFibresView::AutomaticAnalysis(){
    QString prodPath = directory + mitk::IOUtil::GetDirectorySeparator();
    if(cnnPath.isEmpty()){
        int reply = Ask("Question", "Do you have a previous automatic segmentation?");

        if(reply==QMessageBox::Yes){
            cnnPath = QFileDialog::getOpenFileName(NULL, "Open Automatic Segmentation file",
            directory.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());
            MITK_INFO << ("Loaded file: " + cnnPath).toStdString();
        } else{
            QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
            if (nodes.size() > 1) {
                QMessageBox::warning(NULL, "Attention",
                "Please select the CEMRA image from the Data Manager to segment!");
                return;
            }

            QString mraPath = prodPath + "test.nii";
            if(nodes.size() == 0){
                MITK_INFO << "Searching MRA files";

                QString imagePath = QFileDialog::getOpenFileName(NULL, "Open MRA image file",
                directory.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());

                if(!QFile::exists(mraPath)){
                    if(!QFile::copy(imagePath, mraPath)){
                        MITK_ERROR << "MRA could not be prepared for automatic segmentation.";
                        return;
                    }
                }
            } else if(nodes.size() == 1){
                mitk::BaseData::Pointer data = nodes.at(0)->GetData();
                bool successfulImageLoading = false;
                if (data) {
                    mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
                    if (image) {
                        mitk::IOUtil::Save(image, mraPath.toStdString());
                        successfulImageLoading = true;
                    }
                }

                if(!successfulImageLoading){
                    MITK_ERROR << "MRA image could not be loaded";
                    return;
                }
            }

            MITK_INFO << "[AUTOMATIC_PIPELINE] CemrgNet (auto segmentation)";
            std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
            cmd->SetUseDockerContainers(true);
            cnnPath = cmd->DockerCemrgNetPrediction(mraPath);
            remove(mraPath.toStdString().c_str());
        }
    }
    CemrgCommonUtils::RoundPixelValues(cnnPath);
    QFileInfo fi(cnnPath);
    if(directory.compare(fi.absolutePath())==0){
        MITK_INFO << ("[ATTENTION] Changing working directory to: " + fi.absolutePath());
        directory = fi.absolutePath();
    }
    MITK_INFO << ("Directory: "+ directory + " segmentation name: " + fi.baseName()).toStdString();

    MITK_INFO << "[AUTOMATIC_PIPELINE] Clean Segmentation";
    ImageType::Pointer segImage = atrium->CleanAutomaticSegmentation(directory, (fi.baseName()+".nii"), "prodClean.nii");
    cnnPath = prodPath + "LA.nii";
    QString laregPath = prodPath + "LA-reg.nii";
    mitk::IOUtil::Save(mitk::ImportItkImage(segImage), cnnPath.toStdString());

    int replyLAreg = Ask("Register to LGE", "Move to LGE space and prepare for scar projection?");
    if(replyLAreg==QMessageBox::Yes){
        MITK_INFO << "search for lge and mra";
        QString lgePath = GetFilePath("dcm-LGE", ".nii");
        QString mraPath = GetFilePath("dcm-MRA", ".nii");

        if(lgePath.isEmpty()){
            lgePath = QFileDialog::getOpenFileName(NULL, "Open LGE image file",
                directory.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());
        }
        if(mraPath.isEmpty()){
            mraPath = QFileDialog::getOpenFileName(NULL, "Open MRA image file",
                directory.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());
        }
        // register and transform
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->ExecuteRegistration(directory, lgePath, mraPath); // rigid.dof is the default name
        cmd->ExecuteTransformation(directory, cnnPath, laregPath);
        //update cnnPath

        segImage = atrium->LoadImage(laregPath);
        cnnPath = laregPath;
        analysisOnLge = true;
        tagName += "-reg";
    }

    mitk::IOUtil::Load(cnnPath.toStdString(), *this->GetDataStorage());

    MITK_INFO << "[AUTOMATIC_PIPELINE] Assign PV labels automatically";
    ImageType::Pointer tagAtriumAuto = atrium->AssignAutomaticLabels(segImage, directory);

    MITK_INFO << "[AUTOMATIC_PIPELINE] Create Mesh";
    //Ask for user input to set the parameters
    bool userInputsAccepted = GetUserMeshingInputs();

    if(userInputsAccepted){

        atrium->ProjectTagsOnSurface(tagAtriumAuto, directory, tagName+".vtk",
            uiMesh_th, uiMesh_bl, uiMesh_smth, uiMesh_ds, true);

        MITK_INFO << "[AUTOMATIC_PIPELINE] Add the mesh to storage";
        QString path = prodPath + tagName + ".vtk";

        std::cout << "Path to load: " << path.toStdString() <<'\n';
        std::cout << "tagName: " << tagName.toStdString() << '\n';
        mitk::Surface::Pointer surface = mitk::IOUtil::Load<mitk::Surface>(path.toStdString());
        // mitk::Surface::Pointer surface = CemrgCommonUtils::LoadVTKMesh(path.toStdString());
        // CemrgCommonUtils::FlipXYPlane(surface, "");

        std::string meshName = tagName.toStdString() + "-Mesh";
        CemrgCommonUtils::AddToStorage(surface, meshName, this->GetDataStorage());
    }

}

void AtrialFibresView::MeshPreprocessing(){
    MITK_INFO << "[MeshPreprocessing] ";
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    //Show the plugin
    this->GetSite()->GetPage()->ResetPerspective();
    AtrialFibresClipperView::SetDirectoryFile(directory, tagName+".vtk", automaticPipeline);
    this->GetSite()->GetPage()->ShowView("org.mitk.views.atrialfibresclipperview");
}

// Manual pipeline
void AtrialFibresView::SegmentIMGS() {
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.
    QString prodPath = directory + mitk::IOUtil::GetDirectorySeparator();

    int replyAuto = Ask("Question", "Do you want an automatic segmentation?");
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
                mraPath = prodPath + "test.nii";
                mitk::IOUtil::Save(image, mraPath.toStdString());
                std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
                cmd->SetUseDockerContainers(true);
                if(cnnPath.isEmpty()){
                    cnnPath = cmd->DockerCemrgNetPrediction(mraPath);
                }
                CemrgCommonUtils::RoundPixelValues(cnnPath);
                mitk::ProgressBar::GetInstance()->Progress();

                //Clean prediction
                ImageType::Pointer segImage = atrium->RemoveNoiseFromAutomaticSegmentation(directory);
                cnnPath = prodPath + "LA.nii";

                mitk::IOUtil::Save(mitk::ImportItkImage(segImage), cnnPath.toStdString());
                mitk::IOUtil::Load(cnnPath.toStdString(), *this->GetDataStorage());
                remove(mraPath.toStdString().c_str());

                QMessageBox::information(NULL, "Attention", "Command Line Operations Finished!");
                mitk::ProgressBar::GetInstance()->Progress();
                this->BusyCursorOff();

            } else
                QMessageBox::warning(NULL, "Attention", "Please select a CEMRA to segment!");
        }//_if_data
    } else {
        int replyLoad =  Ask("Question", "Do you have a segmentation to load?");

        if (replyLoad == QMessageBox::Yes) {
            QString path = QFileDialog::getOpenFileName(NULL, "Open Segmentation file",
                            directory.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());

            if (path.isEmpty()) return;

            mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
            mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());

            //Restore image name
            QFileInfo fullPathInfo(path);
        } else {
            //Show segmentation plugin
            this->GetSite()->GetPage()->ShowView("org.mitk.views.segmentation");
        } //_if_q_loadSegmentation
    }//_if_q_automatic
}

void AtrialFibresView::IdentifyPV(){
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.
    QString path =  directory + mitk::IOUtil::GetDirectorySeparator() + "segmentation.vtk";

    if(!QFileInfo::exists(path)){
        // Select segmentation from data manager
        QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
        if (nodes.size() != 1) {
            QMessageBox::warning(NULL, "Attention",
                "Please select the loaded or created segmentation to create a surface!");
            return;
        }

        //Find the selected node
        mitk::DataNode::Pointer segNode = nodes.at(0);
        mitk::BaseData::Pointer data = segNode->GetData();
        if (data) {
            mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
            if (image) {

                MITK_INFO << "Mesh and save as 'segmentation.vtk'";
                ImageType::Pointer segImage = ImageType::New();
                mitk::CastToItkImage(image, segImage);
                bool userInputsAccepted = GetUserMeshingInputs();
                if(userInputsAccepted){
                    atrium->SurfSegmentation(segImage, directory, "segmentation.vtk", uiMesh_th, uiMesh_bl, uiMesh_smth, uiMesh_ds);

                    //Add the mesh to storage
                    std::string meshName = segNode->GetName() + "-Mesh";
                    CemrgCommonUtils::AddToStorage(
                        CemrgCommonUtils::LoadVTKMesh(path.toStdString()), meshName, this->GetDataStorage());
                }
            }
        }
    }

    MITK_INFO << "Loading org.mitk.views.atrialfibresclipperview";
    this->GetSite()->GetPage()->ResetPerspective();
    AtrialFibresClipperView::SetDirectoryFile(directory, "segmentation.vtk", automaticPipeline);
    this->GetSite()->GetPage()->ShowView("org.mitk.views.atrialfibresclipperview");
}

void AtrialFibresView::CreateLabelledMesh(){

}

void AtrialFibresView::ClipperMV(){
    if(automaticPipeline){
        MITK_INFO << "[ClipperMV] Clipping mitral valve";
        atrium->ClipMitralValveAuto(directory, "prodMVI.nii", tagName+".vtk");
    } else {

    }

}

void AtrialFibresView::ClipperPV(){
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.#
    if (!LoadSurfaceChecks()) return; // Surface was not loaded and user could not find file.

    QString prodPath = directory + mitk::IOUtil::GetDirectorySeparator();

    if(automaticPipeline){
        MITK_INFO << "[ClipperPV] clipping PVs from automatic pipeline.";

        QString path = prodPath + tagName + ".vtk";

        mitk::Surface::Pointer surface = mitk::IOUtil::Load<mitk::Surface>(path.toStdString());
        path = prodPath + "prodClipperIDsAndRadii.txt";
        if(QFile::exists(path)){

            MITK_INFO << "[ClipperPV] Reading in centre IDs and radii";
            std::ifstream fi;
            fi.open((path.toStdString()));

            int numPts;
            fi >> numPts;

            for (int ix = 0; ix < numPts; ix++) {
                int ptId;
                double radius, x_c, y_c, z_c;

                fi >> ptId >> x_c >> y_c >> z_c >> radius;

                QString spherePath = prodPath + "pvClipper_" + QString::number(ptId) + ".vtk";
                std::cout << "Read points: ID:" << ptId << " C=[" << x_c << " " << y_c << " " << z_c << "] R=" << radius <<'\n';

                surface = CemrgCommonUtils::ClipWithSphere(surface, x_c, y_c, z_c, radius, spherePath);
            }
            fi.close();

            // save surface
            path = prodPath + tagName + ".vtk";
            mitk::IOUtil::Save(surface, path.toStdString());
            atrium->SetSurface(path);

            SetTagNameFromPath(path);
            ClipperMV();

            QString correctLabels = prodPath + "prodSeedLabels.txt";
            QString naiveLabels = prodPath + "prodNaiveSeedLabels.txt";
            atrium->SetSurfaceLabels(correctLabels, naiveLabels);

            QMessageBox::information(NULL, "Attention", "Clipping of PV and MV finished");

        } else{
            QMessageBox::warning(NULL, "Warning", "Radii file not found");
            return;
        }

    } else {
        MITK_INFO << "[ClipperPV] clipping PVs from manual pipeline.";
    }
}

// Labelled Mesh to UAC
void AtrialFibresView::SelectLandmarks(){
    MITK_INFO << "[MeshPreprocessing] ";
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.
    if (!LoadSurfaceChecks()) return;

    //Show the plugin
    this->GetSite()->GetPage()->ResetPerspective();
    AtrialFibresLandmarksView::SetDirectoryFile(directory, tagName+".vtk");
    this->GetSite()->GetPage()->ShowView("org.mitk.views.atrialfibreslandmarksview");
}

void AtrialFibresView::MeshingOptions(){
    if (!RequestProjectDirectoryFromUser()) return;
    if (!LoadSurfaceChecks()) return;

    QString prodPath = directory + mitk::IOUtil::GetDirectorySeparator();

    MITK_INFO << "[MeshingOptions] Remeshing";
    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);
    bool userInputsAccepted = GetUserRemeshingInputs();
    if(userInputsAccepted){
        QString refinedPath = cmd->DockerRemeshSurface(directory, tagName, tagName+refinedSuffix, uiRemesh_max, uiRemesh_min, uiRemesh_avrg, uiRemesh_surfcorr);

        atrium->SetSurface(refinedPath);

        if(uiRemesh_extractParts){
            MITK_INFO << "[MeshingOptions] save individual vtks for: RS, RI, LS, LI and LAA";
            atrium->ExtractLabelFromShell(directory, atrium->LSPV(), "LSPV");
            atrium->ExtractLabelFromShell(directory, atrium->LIPV(), "LIPV");
            atrium->ExtractLabelFromShell(directory, atrium->RSPV(), "RSPV");
            atrium->ExtractLabelFromShell(directory, atrium->RIPV(), "RIPV");
            atrium->ExtractLabelFromShell(directory, atrium->LAAP(), "LAAP");
        }

        MITK_INFO << "[MeshingOptions] finished";
    }
}


void AtrialFibresView::UacCalculation(){
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    QString pathRoughLandmark = LandmarkFilesCreated("prodRoughLandmarks", "ROUGH");
    QString pathRefinedLandmark = LandmarkFilesCreated("prodRefinedLandmarks", "REFINED");

    if(pathRoughLandmark.compare("FILE_NOT_FOUND")==0 || pathRefinedLandmark.compare("FILE_NOT_FOUND")==0){
        return;
    }

}

void AtrialFibresView::UacMeshRefinement(){

}

void AtrialFibresView::ConvertFormat(){
    if (!RequestProjectDirectoryFromUser()) return;

    QString meshPath = QFileDialog::getOpenFileName(NULL, "Open Mesh file",
        directory.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());

    QString prodPath = directory + mitk::IOUtil::GetDirectorySeparator();

    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);

    QFileInfo fi(meshPath);
    bool userInputsAccepted = GetUserConvertFormatInputs(fi.baseName(), fi.suffix());
    if(userInputsAccepted){
        MITK_INFO << "[ConvertFormat] Creating CARP output in microns";
        cmd->DockerConvertMeshFormat(directory, fi.baseName(), "vtk", uiFormat_outName, uiFormat_outExt, uiFormat_scale);
    }
}

void AtrialFibresView::ScarProjection(){

    if (!RequestProjectDirectoryFromUser()) return;
    QString prodPath = directory + mitk::IOUtil::GetDirectorySeparator();

    if(!analysisOnLge){
        tagName += "-reg";
        if(!LoadSurfaceChecks()){
            QMessageBox::information(NULL, "Attention", "Repeat the pipeline on the LGE image!");
            return;
        }
    }
    atrium->SetSurface(prodPath + tagName + ".vtk");

    MITK_INFO << "[SCAR_PROJECTION][1] get lge path";
    QString lgePath = GetFilePath("dcm-LGE", ".nii");
    if(lgePath.isEmpty()){
        lgePath = QFileDialog::getOpenFileName(NULL, "Open LGE image file",
            directory.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());
    }
    typedef itk::Image<float, 3> FloatImageType;
    FloatImageType::Pointer lgeFloat = FloatImageType::New();
    mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>(lgePath.toStdString()), lgeFloat);
    mitk::Image::Pointer lge = mitk::ImportItkImage(lgeFloat);

    MITK_INFO << "[SCAR_PROJECTION][2] copy prodClean file to PVeinsCroppedImage";
    QString cleanPath = prodPath + "prodClean.nii";
    QString segPath = prodPath + "PVeinsCroppedImage.nii";
    if(!QFile::copy(cleanPath, segPath)){
        ImageType::Pointer im = atrium->RemoveNoiseFromAutomaticSegmentation(directory);
        mitk::IOUtil::Save(mitk::ImportItkImage(im), cleanPath.toStdString());
        mitk::IOUtil::Save(mitk::ImportItkImage(im), segPath.toStdString());
    }
    atrium->ResampleSegmentationLabelToImage(segPath, lgePath);

    MITK_INFO << "[SCAR_PROJECTION][3] convert tagname to segmentation.vtk";
    QString segvtkName = "segmentation-reg.vtk";
    CemrgCommonUtils::FlipXYPlane(atrium->GetSurface(), directory, segvtkName);
    // mitk::Surface::s mitk::IOUtil::Load<mitk::Surface>

    bool userInputAccepted = GetUserScarProjectionInputs();
    if(userInputAccepted){
        MITK_INFO << "[SCAR_PROJECTION][4] scar projection";
        std::unique_ptr<CemrgScar3D> scar(new CemrgScar3D());
        scar->SetMinStep(uiScar_minStep);
        scar->SetMaxStep(uiScar_maxStep);
        scar->SetMethodType(uiScar_projectionMethod);

        scar->SetScarSegImage(mitk::IOUtil::Load<mitk::Image>(segPath.toStdString()));
        mitk::Surface::Pointer scarShell = scar->Scar3D(directory.toStdString(), lge, (segvtkName).toStdString());

        QString scarPath = prodPath + "MaxScar.vtk";
        CemrgCommonUtils::SetCellDataToPointData(scarShell, scarPath, "scar");

        MITK_INFO << "[SCAR_PROJECTION][5] Thresholding";
        double mean = 0.0, stdv = 0.0;
        bool binarise=true;
        ImageType::Pointer segITK = atrium->LoadImage(segPath, binarise);
        mitk::Image::Pointer roiImage = atrium->ImErode(segITK);

        if(scar->CalculateMeanStd(lge, roiImage, mean, stdv)){
            MITK_INFO << "[SCAR_PROJECTION][6] Saving normalised shell and prodThresholds.txt";
            scar->SaveNormalisedScalars(mean, scarShell, (prodPath + "MaxScar_Normalised.vtk"));
            scar->PrintThresholdingResults(directory, uiScar_thresValues, uiScar_thresholdMethod, mean, stdv);
        } else{
            MITK_INFO << "[SCAR_PROJECTION][ERROR] Wrong dimensions between LGE and ROI Image";
        }

    }
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
        directory = QFileDialog::getExistingDirectory( NULL, "Open Project Directory",
            mitk::IOUtil::GetProgramPath().c_str(),QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);

        MITK_INFO << ("Directory selected:" + directory).toStdString();
        atrium->SetWorkingDirectory(directory);

        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            MITK_WARN << "Please select a project directory with no spaces in the path!";
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            succesfulAssignment = false;
        }//_if

    } else {
        MITK_INFO << ("Project directory already set: " + directory).toStdString();
    }//_if

    return succesfulAssignment;
}

bool AtrialFibresView::GetUserConvertFormatInputs(QString inname, QString inext){
    QDialog* inputs = new QDialog(0,0);
    bool userInputAccepted=false;
    m_UIFormat.setupUi(inputs);
    connect(m_UIFormat.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_UIFormat.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    m_UIFormat.label_1inname->setText(inname);
    m_UIFormat.label_2inext->setText(inext);
    int dialogCode = inputs->exec();

    if (dialogCode == QDialog::Accepted){
        uiFormat_outName = inname;
        uiFormat_outExt = m_UIFormat.comboBox->currentText();

        QString newoutname = m_UIFormat.lineEdit_1outname->text().simplified();
        if(!newoutname.isEmpty()){
            uiFormat_outName = newoutname;
        }

        bool ok1;
        uiFormat_scale = m_UIFormat.lineEdit_2scale->text().toInt(&ok1);
        if(!ok1){
            uiFormat_scale = 1000;
        }
        inputs->deleteLater();
        userInputAccepted=true;

    } else if (dialogCode == QDialog::Rejected) {
        inputs->close();
        inputs->deleteLater();
    }//_if

    return userInputAccepted;
}

bool AtrialFibresView::GetUserRemeshingInputs(){
    QDialog* inputs = new QDialog(0,0);
    bool userInputAccepted=false;
    m_UIRemesh.setupUi(inputs);
    connect(m_UIRemesh.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_UIRemesh.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    int dialogCode = inputs->exec();

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {
        uiRemesh_vtk2carp = m_UIRemesh.check_1vtk2carp->isChecked();
        uiRemesh_extractParts = m_UIRemesh.check_2extractlabels->isChecked();

        bool ok1, ok2, ok3, ok4;
        uiRemesh_max= m_UIRemesh.lineEdit_1max->text().toDouble(&ok1);
        uiRemesh_avrg= m_UIRemesh.lineEdit_2avrg->text().toDouble(&ok2);
        uiRemesh_min= m_UIRemesh.lineEdit_3min->text().toDouble(&ok3);
        uiRemesh_surfcorr= m_UIRemesh.lineEdit_4surfcorr->text().toDouble(&ok4);

        QString newsuffix = m_UIRemesh.lineEdit_5suffix->text().simplified();
        if(!newsuffix.isEmpty()){
            refinedSuffix = "-" + newsuffix;
        }

        //Set default values
        if (!ok1 || !ok2 || !ok3 || !ok4)
            QMessageBox::warning(NULL, "Attention", "Reverting to default parameters!");
        if (!ok1) uiRemesh_max=0.5;
        if (!ok2) uiRemesh_avrg=0.3;
        if (!ok3) uiRemesh_min=0.1;
        if (!ok4) uiRemesh_surfcorr=0.95;

        // checks on min,max and avrg
        if(uiRemesh_max < uiRemesh_min || uiRemesh_max < uiRemesh_avrg || uiRemesh_avrg < uiRemesh_min){
            QMessageBox::warning(NULL, "Attention", "Non consistent parameters, using defaults!");
            uiRemesh_max=0.5;
            uiRemesh_avrg=0.3;
            uiRemesh_min=0.1;
            uiRemesh_surfcorr=0.95;
        }

        inputs->deleteLater();
        userInputAccepted=true;

    } else if (dialogCode == QDialog::Rejected) {
        inputs->close();
        inputs->deleteLater();
    }//_if

    return userInputAccepted;
}

bool AtrialFibresView::GetUserMeshingInputs(){
    QDialog* inputs = new QDialog(0,0);
    bool userInputAccepted=false;
    m_UIMeshing.setupUi(inputs);
    connect(m_UIMeshing.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_UIMeshing.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    int dialogCode = inputs->exec();

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {

        bool ok1, ok2, ok3, ok4;
        uiMesh_th= m_UIMeshing.lineEdit_1->text().toDouble(&ok1);
        uiMesh_bl= m_UIMeshing.lineEdit_2->text().toDouble(&ok2);
        uiMesh_smth= m_UIMeshing.lineEdit_3->text().toDouble(&ok3);
        uiMesh_ds= m_UIMeshing.lineEdit_4->text().toDouble(&ok4);

        //Set default values
        if (!ok1 || !ok2 || !ok3 || !ok4)
            QMessageBox::warning(NULL, "Attention", "Reverting to default parameters!");
        if (!ok1) uiMesh_th=0.5;
        if (!ok2) uiMesh_bl=1.0;
        if (!ok3) uiMesh_smth=3;
        if (!ok4) uiMesh_ds=0.0;

        inputs->deleteLater();
        userInputAccepted=true;

    } else if (dialogCode == QDialog::Rejected) {
        inputs->close();
        inputs->deleteLater();
    }//_if

    return userInputAccepted;
}

bool AtrialFibresView::GetUserScarProjectionInputs(){
    bool userInputAccepted=false;

    QDialog* inputs = new QDialog(0,0);
    m_UIcemrgnet.setupUi(inputs);
    connect(m_UIcemrgnet.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_UIcemrgnet.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    int dialogCode = inputs->exec();

    QString meType_UI;
    QStringList separated_thresh_list;

    uiScar_minStep = -1;
    uiScar_maxStep = 3;
    uiScar_projectionMethod = 2;
    uiScar_thresholdMethod = 1;

    if (dialogCode == QDialog::Accepted) {

        MITK_INFO << "[UI] User inputs being selected.";
        MITK_INFO << "[UI] Intensity projection";
        bool ok1, ok2;
        uiScar_minStep = m_UIcemrgnet.minStep_lineEdit->text().toInt(&ok1);
        uiScar_maxStep = m_UIcemrgnet.maxStep_lineEdit->text().toInt(&ok2);
        if (!ok1) uiScar_minStep = -1;
        if (!ok2) uiScar_maxStep = 3;

        uiScar_projectionMethod = m_UIcemrgnet.maxProjection_radioButton->isChecked() ? 2 : 1;
        uiScar_thresholdMethod = m_UIcemrgnet.iir_radioButton->isChecked() ? 1 : 2;
        meType_UI = m_UIcemrgnet.maxProjection_radioButton->isChecked() ? "Max" : "Mean";

        MITK_INFO << ("[UI] Using: " + meType_UI + " Intensity projection.").toStdString();
        MITK_INFO << ("[UI] In/out values: (" + QString::number(uiScar_minStep) + ", " +
                      QString::number(uiScar_maxStep) + ")").toStdString();
        MITK_INFO << QString::number(uiScar_projectionMethod);
        MITK_INFO << "[UI] Thresholding information.";

        QString thresh_list, whichThresh;

        if (m_UIcemrgnet.iir_radioButton->isChecked()) { // IIR method
            whichThresh = "IIR";
            thresh_list = m_UIcemrgnet.iir_textEdit->toPlainText();
            separated_thresh_list << "0.97" << "1.16";
        } else if (m_UIcemrgnet.meanSD_radioButton->isChecked()) { // SDev method
            whichThresh = "MEAN+SD";
            thresh_list = m_UIcemrgnet.meanSD_textEdit->toPlainText();
            separated_thresh_list << "2.3" << "3.3";
        }

        MITK_INFO << ("[UI] Threshold: " + whichThresh).toStdString();
        MITK_INFO << ("[UI] Threshold list: " + thresh_list).toStdString();
        MITK_INFO << QString::number(uiScar_thresholdMethod);

        thresh_list.remove(" ", Qt::CaseSensitive);
        if (!thresh_list.isEmpty()) {

            MITK_INFO << "[UI] Creating list of thresholds";
            separated_thresh_list.removeLast();
            separated_thresh_list.removeLast();
            separated_thresh_list = thresh_list.split("," , QString::SkipEmptyParts);
            int listspaces = separated_thresh_list.removeAll(" ");
            int listduplicates = separated_thresh_list.removeDuplicates();
            separated_thresh_list.sort();
        }//_if

        double tryNumber;
        bool vOK;
        for(int ix=0; ix<separated_thresh_list.size(); ix++) {
            MITK_INFO << separated_thresh_list.at(ix);
            tryNumber = separated_thresh_list.at(ix).toDouble(&vOK);
            if (vOK) uiScar_thresValues.push_back(tryNumber);
        }

        inputs->close();
        inputs->deleteLater();
        userInputAccepted=true;

    } else if (dialogCode == QDialog::Rejected) {
        inputs->close();
        inputs->deleteLater();
    }//_if

    return userInputAccepted;
}

QString AtrialFibresView::LandmarkFilesCreated(QString defaultName, QString type){
    QString prodPath, res;
    prodPath = directory + mitk::IOUtil::GetDirectorySeparator();
    res = "FILE_NOT_FOUND";

    bool foundVtk, foundTxt;
    foundVtk = QFile::exists(prodPath+defaultName+".vtk");
    foundTxt = QFile::exists(prodPath+defaultName+".txt");

    MITK_INFO(foundVtk) << ("Found" + type + "file in VTK format").toStdString();
    MITK_INFO(foundTxt) << ("Found" + type + "file in TXT format").toStdString();

    std::string msg;
    if(!foundTxt && !foundVtk){
        msg = "File not found\n";
        msg += "Do you have a VTK or TXT for the: ";
        msg += (type.toStdString() + " landmarks?");

        int reply = Ask("File not found", msg);
        if(reply==QMessageBox::Yes){
            msg = ("Open " + type.toStdString() + " landmarks file");
            res = QFileDialog::getOpenFileName(NULL, msg.c_str(), directory.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());
        } else{
            msg = ("Use StepW button to Select Landmarks of type: " + type).toStdString();
            QMessageBox::information(NULL, "Attention", msg.c_str());
        }
    } else if(foundVtk) {
        MITK_INFO << "Loading file in VTK format";
        res = prodPath+defaultName+".vtk";
    } else{
        MITK_INFO << "Loading file in TXT format";
        res = prodPath+defaultName+".txt";
    }

    return res;
}

void AtrialFibresView::SetManualModeButtons(bool b){
    //Set visibility of buttons
    m_Controls.button_man4_segmentation->setVisible(b);
    m_Controls.button_man5_idPV->setVisible(b);
    m_Controls.button_man6_labelmesh->setVisible(b);
    m_Controls.button_man7_clipMV->setVisible(b);
    m_Controls.button_man8_clipPV->setVisible(b);


    m_Controls.button_0_landmarks->setText("    Step9: Select Landmarks");
    m_Controls.button_0_calculateUac->setText("    Step10: Calculate UAC");
    m_Controls.button_0_refineUac->setText("    Step11: Refine UAC");
}

void AtrialFibresView::SetAutomaticModeButtons(bool b){
    m_Controls.button_auto4_meshpreproc->setVisible(b);
    m_Controls.button_auto5_clipPV->setVisible(b);

    m_Controls.button_0_landmarks->setText("    Step6: Select Landmarks");
    m_Controls.button_0_calculateUac->setText("    Step7: Calculate UAC");
    m_Controls.button_0_refineUac->setText("    Step8: Refine UAC");
}

void AtrialFibresView::SetTagNameFromPath(QString path){
    QFileInfo fi(path);
    tagName = fi.baseName();
}

bool AtrialFibresView::LoadSurfaceChecks(){
    bool success = true;
    QString prodPath = directory + mitk::IOUtil::GetDirectorySeparator() + tagName + ".vtk";

    if(!QFile::exists(prodPath)){
        int reply1 = Ask("Surface file not found", "Do you have a surface file to load?");
        if(reply1==QMessageBox::Yes){
            UserLoadSurface();
        } else {
            success=false;
        }
    } else{
        std::string msg = ("Load automatically file called: [" + tagName + ".vtk]?").toStdString();
        int reply2 = Ask("Surface file to load", msg.c_str());
        if(reply2==QMessageBox::No){
            UserLoadSurface();
        }
    }

    MITK_INFO << ("[LoadSurfaceChecks] Loading surface: " + prodPath).toStdString();

    return success;
}

void AtrialFibresView::UserLoadSurface(){
    QString newpath = QFileDialog::getOpenFileName(NULL, "Select the surface file to load!", directory.toStdString().c_str(), QmitkIOUtil::GetFileOpenFilterString());
    SetTagNameFromPath(newpath);
}

int AtrialFibresView::Ask(std::string title, std::string msg){
    return QMessageBox::question(NULL, title.c_str(), msg.c_str(), QMessageBox::Yes, QMessageBox::No);
}

QString AtrialFibresView::GetFilePath(QString nameSubstring, QString extension){
    QString result = "";
    QDirIterator searchit(directory, QDirIterator::Subdirectories);
    while(searchit.hasNext()) {
        QFileInfo searchfinfo(searchit.next());
        if (searchfinfo.fileName().contains(extension, Qt::CaseSensitive)) {
            if (searchfinfo.fileName().contains((nameSubstring), Qt::CaseSensitive))
                result = searchfinfo.absoluteFilePath();
        }
    }//_while

    return result;
}
