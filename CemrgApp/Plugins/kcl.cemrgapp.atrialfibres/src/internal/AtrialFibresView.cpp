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
#include <mitkLog.h>

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
#include "AtrialFibresVisualiseView.h"

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
#include <vtkBooleanOperationPolyDataFilter.h>
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
#include <QDate>

//CemrgAppModule
#include <CemrgCommonUtils.h>
#include <CemrgCommandLine.h>
#include <CemrgMeasure.h>
#include <CemrgScar3D.h>

const std::string AtrialFibresView::VIEW_ID = "org.mitk.views.atrialfibresview"; // org.mitk.views.atrialfibresview

void AtrialFibresView::CreateQtPartControl(QWidget *parent) {

    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, SIGNAL(clicked()), this, SLOT(LoadDICOM()));
    connect(m_Controls.button_2, SIGNAL(clicked()), this, SLOT(ProcessIMGS()));

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
    connect(m_Controls.button_0_fibreMapUac, SIGNAL(clicked()), this, SLOT(UacFibreMapping()));

    // Scar analysis
    connect(m_Controls.button_z_scar, SIGNAL(clicked()), this, SLOT(ScarProjection()));

    connect(m_Controls.button_r, SIGNAL(clicked()), this, SLOT(Reset()));

    //Set visibility of buttons
    m_Controls.button_2_1->setVisible(false);
    connect(m_Controls.button_2_1, SIGNAL(clicked()), this, SLOT(ConvertNII()));

    m_Controls.button_man4_1_Lge->setVisible(false);
    connect(m_Controls.button_man4_1_Lge, SIGNAL(clicked()), this, SLOT(ManualLgeRegistration()));

    m_Controls.button_man4_2_postproc->setVisible(false);
    connect(m_Controls.button_man4_2_postproc, SIGNAL(clicked()), this, SLOT(SegmentationPostprocessing()));

    m_Controls.button_man7_1_landmarks->setVisible(false);
    m_Controls.button_man7_2_clipMV->setVisible(false);
    connect(m_Controls.button_man7_1_landmarks, SIGNAL(clicked()), this, SLOT(SelectMvLandmarks()));
    connect(m_Controls.button_man7_2_clipMV, SIGNAL(clicked()), this, SLOT(ClipMV()));

    m_Controls.button_0_1_uacRough->setVisible(false);
    m_Controls.button_0_2_uacRefined->setVisible(false);

    connect(m_Controls.button_0_1_uacRough, SIGNAL(clicked()), this, SLOT(UacCalculationRough()));
    connect(m_Controls.button_0_2_uacRefined, SIGNAL(clicked()), this, SLOT(UacCalculationRefined()));
    connect(m_Controls.button_0_3_checklabels, SIGNAL(clicked()), this, SLOT(UacCalculationVerifyLabels()));

    // Set default variables
    tagName = "Labelled";
    refinedSuffix = "-refined";
    resurfaceMesh = false;
    uiRemesh_isscalar = false;
    uiRemesh_extractParts = false;
    uiRemesh_cleanmesh = true;
    userHasSetMeshingParams = false;
    atrium = std::unique_ptr<CemrgAtrialTools>(new CemrgAtrialTools());
    // atrium->SetDebugModeOn();
    atrium->SetDebugModeOff();

    SetManualModeButtonsOff();
    SetAutomaticModeButtonsOff();
    m_Controls.button_man7_clipMV->setEnabled(false);
    m_Controls.button_man8_clipPV->setText("    Step8: Clip PVs and/or MV");

    uiSelector_pipeline = 0;
    uiSelector_imgauto_skipCemrgNet = false;
    uiSelector_imgauto_skipLabel = false;
    uiSelector_img_scar = false;
    uiSelector_man_useCemrgNet = false;
    SetLgeAnalysis(false);
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

                            path = niftiFolder.absolutePath() + "/" + thisFile;
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
        std::string msg = "Cannot find the type of images automatically.";
        msg += "Revert to user order and selections in the data manager:\n LGE at the top, then CEMRA at the bottom!";
        QMessageBox::warning(NULL, "Attention", msg.c_str());
        index.resize(nodes.size());
        std::iota(index.begin(), index.end(), 0);
    }//_if

    //Convert to Nifti
    int ctr = 0;
    QString prodPath, type;
    bool resampleImage, reorientToRAI;
    resampleImage = true;
    reorientToRAI = true;

    this->BusyCursorOn();
    mitk::ProgressBar::GetInstance()->AddStepsToDo(index.size());
    foreach (int idx, index) {
        type = (ctr==0) ? "LGE":"MRA";
        prodPath = directory + "/" + "dcm-" + type + "-" + seriesDscrps.at(idx).c_str() + ".nii";
        bool successfulNitfi = CemrgCommonUtils::ConvertToNifti(nodes.at(idx)->GetData(), prodPath, resampleImage, reorientToRAI);
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
    MITK_INFO << "[AnalysisChoice]";
    bool testRpdfu = RequestProjectDirectoryFromUser();
    MITK_INFO(testRpdfu) << "Should continue";
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.
    MITK_INFO << "After directory selection checkup";

    bool userInputsAccepted = GetUserAnalysisSelectorInputs();
    if(userInputsAccepted){
        // uiSelector_pipeline; // =0 (imgAuto), =1 (imgManual), =2 (surf)
        SetAutomaticPipeline(uiSelector_pipeline==0); // Image Pipeline (automatic)
        if(uiSelector_pipeline==0){
            MITK_INFO<<"[AnalysisChoice] Automatic pipeline";
            SetAutomaticModeButtonsOn();
            SetManualModeButtonsOff();
            m_Controls.button_3_imanalysis->setText("    Step3: Image Analysis");
            if(!uiSelector_imgauto_skipLabel){
                MITK_INFO << "[AnalysisChoice] Performing Automatic analysis on CemrgNet prediction.";
                AutomaticAnalysis();
            } else{
                MITK_INFO << "[AnalysisChoice] Skipping Neural network prediction and labelling. Checking for vtk file.";
                if(!LoadSurfaceChecks()) return;
            }
        } else if(uiSelector_pipeline==1){ // Image Pipeline (manual)
            MITK_INFO << "[AnalysisChoice] Setting up manual analysis.";
            SetManualModeButtonsOn();
            SetAutomaticModeButtonsOff();
            m_Controls.button_man4_segmentation->setEnabled(true);
            m_Controls.button_auto4_meshpreproc->setVisible(true);
            m_Controls.button_man4_2_postproc->setVisible(true);

        } else if(uiSelector_pipeline==2){ // Surface Pipeline (manual)
            MITK_INFO << "[AnalysisChoice] Analysis starting from surface";
            SetAutomaticPipeline(false);

            // Load surface mesh
            if(!LoadSurfaceChecks()) return;

            // Create fake segmentation image for labelling
            double origin[3] = {0, 0, 0};
            double spacing[3] = {1, 1, 1};
            CemrgCommonUtils::SaveImageFromSurfaceMesh(Path(tagName+".vtk"), origin, spacing);
            CemrgCommonUtils::SavePadImageWithConstant(Path(tagName+".nii"));

            mitk::Image::Pointer im = CemrgCommonUtils::ReturnBinarised(mitk::IOUtil::Load<mitk::Image>(StdStringPath(tagName+".nii")));
            // CemrgCommonUtils::Binarise(im);
            mitk::IOUtil::Save(im, StdStringPath(tagName+".nii"));
            CemrgCommonUtils::AddToStorage(im, tagName.toStdString(), this->GetDataStorage());

            SetManualModeButtonsOn();
            SetAutomaticModeButtonsOff();
            m_Controls.button_man4_segmentation->setEnabled(false);
            m_Controls.button_auto4_meshpreproc->setVisible(true);
            m_Controls.button_man4_2_postproc->setVisible(true);
        }
    }
}

// Automatic pipeline
void AtrialFibresView::AutomaticAnalysis(){
    MITK_INFO << "TIMELOG|AutomaticAnalysis| Start";
    QString prodPath = directory + "/";
    if(cnnPath.isEmpty()){
        if(uiSelector_imgauto_skipCemrgNet){
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
    mitk::IOUtil::Save(mitk::ImportItkImage(segImage), cnnPath.toStdString());
    cnnPath = UserIncludeLgeAnalysis(cnnPath, segImage);

    mitk::IOUtil::Load(cnnPath.toStdString(), *this->GetDataStorage());

    MITK_INFO << "[AUTOMATIC_PIPELINE] Assign PV labels automatically";
    ImageType::Pointer tagAtriumAuto = atrium->AssignAutomaticLabels(segImage, directory);

    MITK_INFO << "[AUTOMATIC_PIPELINE] Create Mesh";
    //Ask for user input to set the parameters
    bool userInputsAccepted = GetUserMeshingInputs();

    if(userInputsAccepted){
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->SetUseDockerContainers(true);

        cmd->ExecuteSurf(directory, Path("prodClean.nii"), uiMesh_th, uiMesh_bl, uiMesh_smth);
        atrium->ProjectTagsOnExistingSurface(tagAtriumAuto, directory, tagName+".vtk");

        MITK_INFO << "[AUTOMATIC_PIPELINE] Add the mesh to storage";
        QString path = prodPath + tagName + ".vtk";

        std::cout << "Path to load: " << path.toStdString() <<'\n';
        std::cout << "tagName: " << tagName.toStdString() << '\n';
        mitk::Surface::Pointer surface = mitk::IOUtil::Load<mitk::Surface>(path.toStdString());

        std::string meshName = tagName.toStdString() + "-Mesh";
        CemrgCommonUtils::AddToStorage(surface, meshName, this->GetDataStorage());

        std::string msg = "Created labelled mesh - Complete\n\n";
        msg += "Automatic labels were assigned to the mesh. \n";
        msg += "To set default labels: \n\n";
        msg += "+ Identify the LAA and PVs in Step4: Mesh Preprocessing, and set the clippers.\n";
        msg += "+ Click Step5: Clip PVs and MV (this step will set the correct labels)";
        QMessageBox::information(NULL, "Automatic Labelling complete", msg.c_str());
    }
    MITK_INFO << "TIMELOG|AutomaticAnalysis| End";
}

void AtrialFibresView::MeshPreprocessing(){
    MITK_INFO << "TIMELOG|MeshPreprocessing| Start";
    MITK_INFO << "[MeshPreprocessing] ";
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.
    if (!LoadSurfaceChecks()) return;

    //Show the plugin
    this->GetSite()->GetPage()->ResetPerspective();
    AtrialFibresClipperView::SetDirectoryFile(directory, tagName+".vtk", true);
    this->GetSite()->GetPage()->ShowView("org.mitk.views.atrialfibresclipperview");
}

// Manual pipeline
void AtrialFibresView::SegmentIMGS() {
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.
    QString prodPath = directory + "/";

    if(uiSelector_man_useCemrgNet){
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
                atrium->SaveImageToDisk(segImage, directory, "LA.nii");
                cnnPath = UserIncludeLgeAnalysis(Path("LA.nii"), segImage);

                mitk::Image::Pointer im = mitk::ImportItkImage(segImage);
                CemrgCommonUtils::Binarise(im);
                mitk::IOUtil::Save(im, cnnPath.toStdString());

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

            QFileInfo fi(path);
            if(!analysisOnLge && fi.baseName().contains("-reg")){
                std::string msg = "Registered segmentation " + fi.baseName().toStdString() + " found.";
                msg += "\nConsider LGE scar projection analysis?";
                int replyLgeAnalysis = Ask("Question", msg);
                if(replyLgeAnalysis==QMessageBox::Yes){
                    analysisOnLge = true;
                    tagName += "-reg";
                }
            }

            ImageType::Pointer segImage = atrium->LoadImage(path);
            path = UserIncludeLgeAnalysis(path, segImage);

            mitk::Image::Pointer im = mitk::ImportItkImage(segImage);
            CemrgCommonUtils::Binarise(im);
            mitk::IOUtil::Save(im, path.toStdString());

            mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
            mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
        } else {
            //Show segmentation plugin
            this->GetSite()->GetPage()->ShowView("org.mitk.views.segmentation");
            m_Controls.button_man4_1_Lge->setVisible(true);
        } //_if_q_loadSegmentation
    }//_if_q_automatic
}

void AtrialFibresView::ManualLgeRegistration(){
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1) {
        QMessageBox::warning(NULL, "Attention",
                "Please select the CEMRA images from the Data Manager to segment!");
        return;
    }//_if

    QString segPath, prodPath;
    prodPath = directory + "/";
    mitk::BaseData::Pointer data = nodes.at(0)->GetData();
    if (data) {
        mitk::Image::Pointer segImageMitk = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (segImageMitk) {

            this->BusyCursorOn();
            mitk::ProgressBar::GetInstance()->AddStepsToDo(2);

            segPath = prodPath + "LA.nii";
            mitk::IOUtil::Save(segImageMitk, segPath.toStdString());

            ImageType::Pointer segImage = ImageType::New();
            mitk::CastToItkImage(segImageMitk, segImage);

            QString segRegPath = UserIncludeLgeAnalysis(segPath, segImage);
            mitk::ProgressBar::GetInstance()->Progress();

            mitk::IOUtil::Save(mitk::ImportItkImage(segImage), segRegPath.toStdString());
            mitk::IOUtil::Load(segRegPath.toStdString(), *this->GetDataStorage());

            QMessageBox::information(NULL, "Attention", "Command Line Operations Finished!");
            mitk::ProgressBar::GetInstance()->Progress();
            this->BusyCursorOff();

        } else{
            QMessageBox::warning(NULL, "Attention", "Please select the segmentation!");
        }
    }//_if_data
}

void AtrialFibresView::SegmentationPostprocessing(){
    this->GetSite()->GetPage()->ShowView("org.mitk.views.segmentationutilities");
}

void AtrialFibresView::IdentifyPV(){
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.
    QString path =  directory + "/" + "segmentation.vtk";

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
                bool userInputsAccepted = GetUserMeshingInputs();
                if(userInputsAccepted){

                    MITK_INFO << "[IdentifyPV] Create clean segmentation";
                    mitk::IOUtil::Save(image, StdStringPath("prodClean.nii"));

                    MITK_INFO << "[IdentifyPV] Create surface file and projecting tags";
                    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
                    cmd->SetUseDockerContainers(true);

                    cmd->ExecuteSurf(directory, Path("prodClean.nii"), uiMesh_th, uiMesh_bl, uiMesh_smth);

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
    MITK_INFO << "TIMELOG|CreateLabelledMesh| Start";
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    if(!analysisOnLge){
        QString segRegPath = GetFilePath("LA-reg", ".nii");
        if(!segRegPath.isEmpty()){
            int reply = Ask("Registered segmentation found", "Consider a Scar Projection analysis?");
            SetLgeAnalysis(reply==QMessageBox::Yes);
        }
        tagName += analysisOnLge ? "-reg" : "";
    }

    QString prodPath =  directory + "/";

    ImageType::Pointer pveins = atrium->LoadImage(prodPath+"PVeinsLabelled.nii");
    if(!tagName.contains("Labelled")){
        std::string msg = "Changing working name from " + tagName.toStdString() + " to 'Labelled'";
        QMessageBox::information(NULL, "Attention", msg.c_str());
        tagName = "Labelled";
    }
    pveins = atrium->AssignOstiaLabelsToVeins(pveins, directory, tagName);

    MITK_INFO << "[CreateLabelledMesh] Create Mesh";
    //Ask for user input to set the parameters
    bool userInputsAccepted = GetUserMeshingInputs();

    if(userInputsAccepted){
        MITK_INFO << "[CreateLabelledMesh] Create clean segmentation";
        mitk::Image::Pointer segIm = mitk::Image::New();
        mitk::CastToMitkImage(pveins, segIm);
        segIm = CemrgCommonUtils::ReturnBinarised(segIm);
        mitk::IOUtil::Save(segIm, StdStringPath("prodClean.nii"));

        MITK_INFO << "[CreateLabelledMesh] Create surface file and projecting tags";
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->SetUseDockerContainers(true);

        cmd->ExecuteSurf(directory, Path("prodClean.nii"), uiMesh_th, uiMesh_bl, uiMesh_smth);
        atrium->ProjectTagsOnExistingSurface(pveins, directory, tagName+".vtk");

        MITK_INFO << "Add the mesh to storage";
        QString path = prodPath + tagName + ".vtk";

        std::cout << "Path to load: " << path.toStdString() <<'\n';
        std::cout << "tagName: " << tagName.toStdString() << '\n';
        mitk::Surface::Pointer surface = mitk::IOUtil::Load<mitk::Surface>(path.toStdString());

        std::string meshName = tagName.toStdString() + "-Mesh";
        CemrgCommonUtils::AddToStorage(surface, meshName, this->GetDataStorage());
    }
    MITK_INFO << "TIMELOG|CreateLabelledMesh| End";
}

void AtrialFibresView::ClipperMV(){
    if(automaticPipeline){
        MITK_INFO << "[ClipperMV] Clipping mitral valve";
        atrium->ClipMitralValveAuto(directory, "prodMVI.nii", tagName+".vtk");
    } else {
        if (m_Controls.button_man7_1_landmarks->isVisible()) {
            m_Controls.button_man7_1_landmarks->setVisible(false);
            m_Controls.button_man7_2_clipMV->setVisible(false);
            return;
        } else {
            MITK_INFO << "[ClipperMV] Reveal buttons for manual functionalities";
            m_Controls.button_man7_1_landmarks->setVisible(true);
            m_Controls.button_man7_2_clipMV->setVisible(true);
        }//_if
    }
}

void AtrialFibresView::SelectMvLandmarks(){
    if (m_Controls.button_man7_1_landmarks->text() == QString::fromStdString("Select Landmarks")) {

        //Show the plugin
        QMessageBox::information(NULL, "Attention", "Please select 3 points around the mitral valve!");
        this->GetSite()->GetPage()->ShowView("org.mitk.views.pointsetinteraction");
        m_Controls.button_man7_1_landmarks->setText("Display Clipper");

    } else {

        //Check for selection of points
        QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
        if (nodes.empty()) {
            QMessageBox::warning(NULL, "Attention", "Please select the pointsets from the Data Manager to clip the mitral valve!");
            return;
        }//_if

        //Check selection type
        mitk::DataNode::Pointer landMarks = nodes.front();
        mitk::PointSet::Pointer pointSet = dynamic_cast<mitk::PointSet*>(landMarks->GetData());
        if (!pointSet || pointSet->GetSize() != 3) {
            QMessageBox::warning(NULL, "Attention", "Please select landmarks with 3 points from the Data Manager to continue!");
            return;
        }//_if

        //If the path was chosen incorrectly returns
        if (!RequestProjectDirectoryFromUser()) return;

        //Reset the button
        m_Controls.button_man7_1_landmarks->setText("Select Landmarks");

        //Retrieve mean and distance of 3 points
        double centre[3];
        double radius = 0;
        radius = CemrgCommonUtils::GetSphereParametersFromLandmarks(pointSet, centre);
        double x_c = centre[0];
        double y_c = centre[1];
        double z_c = centre[2];

        //Create the clipper geometry
        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetCenter(x_c, y_c, z_c);
        sphereSource->SetRadius(radius);
        sphereSource->SetPhiResolution(40);
        sphereSource->SetThetaResolution(40);
        sphereSource->Update();
        mitk::Surface::Pointer mvClipper = mitk::Surface::New();
        mvClipper->SetVtkPolyData(sphereSource->GetOutput());

        //Adjust the data storage
        mitk::DataStorage::SetOfObjects::ConstPointer sob = this->GetDataStorage()->GetAll();
        for (mitk::DataStorage::SetOfObjects::ConstIterator nodeIt = sob->Begin(); nodeIt != sob->End(); ++nodeIt)
            if (nodeIt->Value()->GetName().find("MVClipper") != nodeIt->Value()->GetName().npos)
                this->GetDataStorage()->Remove(nodeIt->Value());
        CemrgCommonUtils::AddToStorage(mvClipper, "MVClipper", this->GetDataStorage(), false);
        sob = this->GetDataStorage()->GetAll();
        for (mitk::DataStorage::SetOfObjects::ConstIterator nodeIt = sob->Begin(); nodeIt != sob->End(); ++nodeIt) {
            if (nodeIt->Value()->GetName().find("MVClipper") != nodeIt->Value()->GetName().npos) {
                nodeIt->Value()->SetProperty("opacity", mitk::FloatProperty::New(0.4));
                nodeIt->Value()->SetProperty("color", mitk::ColorProperty::New(1.0, 0.0, 0.0));
            }//_if
        }//_for

    }//_if
}

void AtrialFibresView::ClipMV(){
    //Check for selection of points
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) {
        QMessageBox::warning(NULL, "Attention", "Please select the pointsets from the Data Manager to clip the mitral valve!");
        this->GetSite()->GetPage()->ResetPerspective();
        return;
    }//_if

    //Check selection type
    mitk::DataNode::Pointer landMarks = nodes.front();
    mitk::PointSet::Pointer pointSet = dynamic_cast<mitk::PointSet*>(landMarks->GetData());
    if (!pointSet || pointSet->GetSize() != 3) {
        QMessageBox::warning(NULL, "Attention", "Please select landmarks with 3 points from the Data Manager to continue!");
        return;
    }//_if

    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.
    if (!LoadSurfaceChecks()) return;

    //Read in and copy
    // QString path = Path(tagName+".vtk");
    // mitk::Surface::Pointer surface = CemrgCommonUtils::LoadVTKMesh(path.toStdString());
    mitk::Surface::Pointer surface = mitk::IOUtil::Load<mitk::Surface>(StdStringPath(tagName+".vtk"));
    if (surface->GetVtkPolyData() == NULL) {
        QMessageBox::critical(NULL, "Attention", "No mesh was found in the project directory!");
        return;
    }//_if
    mitk::IOUtil::Save(surface, StdStringPath(tagName+"-Original.vtk"));

    /*
     * Producibility Test
     **/
    mitk::IOUtil::Save(pointSet, StdStringPath("prodMVCLandmarks.mps"));
    /*
     * End Test
     **/

    double mvc_C[3];
    double mvc_R = 0;
    mvc_R = CemrgCommonUtils::GetSphereParametersFromLandmarks(pointSet, mvc_C);
    std::cout << "MITRAL VALVE CLIPPER: ";
    std::cout << "Centre = (" << mvc_C[0] << ", " << mvc_C[1] << ", " << mvc_C[2] << ") : Radius = " << mvc_R << '\n';

    surface = CemrgCommonUtils::ClipWithSphere(surface, mvc_C[0], mvc_C[1], mvc_C[2], mvc_R);
    mitk::IOUtil::Save(surface, StdStringPath(tagName+".vtk"));

}


void AtrialFibresView::ClipperPV(){
    MITK_INFO << "TIMELOG|MeshPreprocessing| End";
    MITK_INFO << "TIMELOG|ClipperPV| Start";
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.#
    if (!LoadSurfaceChecks()) return; // Surface was not loaded and user could not find file.

    QString prodPath = Path();

    MITK_INFO << "[ClipperPV] clipping PVs.";

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
        if(automaticPipeline){
            ClipperMV();

            QString correctLabels = prodPath + "prodSeedLabels.txt";
            QString naiveLabels = prodPath + "prodNaiveSeedLabels.txt";

            QMessageBox::information(NULL, "Attention", "Attempting to correct automatic labels to default ones");
            atrium->SetSurfaceLabels(correctLabels, naiveLabels);
            atrium->SaveSurface(path.toStdString());
        }

        QMessageBox::information(NULL, "Attention", "Clipping of PV and MV finished");

    } else{
        QMessageBox::warning(NULL, "Warning", "Radii file not found");
        return;
    }

}

// Visualise
void AtrialFibresView::SimpleVisualisation(){
    MITK_INFO << "[Visualise] ";
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.
    if (!LoadSurfaceChecks()) return;

    //Show the plugin
    this->GetSite()->GetPage()->ResetPerspective();
    AtrialFibresVisualiseView::SetDirectoryFile(directory, tagName+".vtk");
    this->GetSite()->GetPage()->ShowView("org.mitk.views.atrialfibresvisualiseview");
}

// Labelled Mesh to UAC
void AtrialFibresView::SelectLandmarks(){
    MITK_INFO << "TIMELOG|SelectLandmarks| Start";
    MITK_INFO << "[UAC_Landmarks] ";
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.
    if (!LoadSurfaceChecks()) return;

    if(GetUserUacOptionsInputs(false)){
        uac_whichAtrium = uiUac_whichAtrium.at(uiUac_whichAtriumIndex);
        MITK_INFO << ("[UAC_Landmarks] Seleted ["+uac_whichAtrium+"] analysis").toStdString();
    } else{
        MITK_INFO << "User cancelled selection of LA/RA selection";
        return;
    }

    //Show the plugin
    this->GetSite()->GetPage()->ResetPerspective();
    AtrialFibresLandmarksView::SetDirectoryFile(directory, tagName+".vtk", uac_whichAtrium);
    this->GetSite()->GetPage()->ShowView("org.mitk.views.atrialfibreslandmarksview");
}

void AtrialFibresView::CleanMeshQuality(){
    if (!RequestProjectDirectoryFromUser()) return;

    MITK_INFO << "[CleanMeshQuality] Select mesh";
    std::string msg = "Open the mesh file (vtk or CARP). \nFor CARP, select the .elem files";
    QMessageBox::information(NULL, "Open Mesh File", msg.c_str());
    QString meshPath = QFileDialog::getOpenFileName(NULL, "Open Mesh file",
        StdStringPath().c_str(), QmitkIOUtil::GetFileOpenFilterString());

    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);

    QFileInfo fi(meshPath);
    QString inExt = (fi.suffix().contains(".elem")) ? "carp_txt" : fi.suffix();
    bool cleanmesh = true;
    bool userInputsAccepted = GetUserConvertFormatInputs(fi.baseName(), inExt, cleanmesh);

    if(userInputsAccepted){
        MITK_INFO << "[CleanMeshQuality] Cleaning mesh";
        // std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->SetUseDockerContainers(true);

        cmd->DockerCleanMeshQuality(directory, fi.baseName(), uiFormat_outName, 0.2, inExt, uiFormat_outExt);
        cmd->DockerCleanMeshQuality(directory, uiFormat_outName, uiFormat_outName, 0.1, "vtk", "vtk_polydata");
    }
}

void AtrialFibresView::MeshingOptions(){
    MITK_INFO << "TIMELOG|SelectLandmarks| End";
    MITK_INFO << "TIMELOG|MeshImprovement| Start";
    if (!RequestProjectDirectoryFromUser()) return;

    QMessageBox::information(NULL, "Open Mesh File", "Open the mesh file (vtk ONLY)");
    QString meshPath = QFileDialog::getOpenFileName(NULL, "Open Mesh file",
        StdStringPath().c_str(), tr("Mesh (*.vtk)"));
    QFileInfo fi(meshPath);

    QString meshName = fi.baseName();
    // QString prodPath = directory + "/";

    MITK_INFO << "[MeshingOptions] Remeshing";
    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);
    bool userInputsAccepted = GetUserRemeshingInputs();
    if(userInputsAccepted){
        QString refinedPath = cmd->DockerRemeshSurface(directory, meshName, meshName+refinedSuffix, uiRemesh_max, uiRemesh_min, uiRemesh_avrg, uiRemesh_surfcorr);

        if(!cmd->IsOutputSuccessful(refinedPath)){
            QMessageBox::warning(NULL, "Attention", "Surface remeshing output unsuccessful");
            MITK_ERROR << "Surface remeshing output unsuccessful";
            return;
        }

        if(uiRemesh_isscalar){
            QString pathNoExt = Path(meshName);
            QString fieldName = "scar";
            CemrgCommonUtils::VtkCellScalarToFile(pathNoExt+".vtk", pathNoExt+"_elem.dat", fieldName);
            CemrgCommonUtils::VtkPointScalarToFile(pathNoExt+".vtk", pathNoExt+"_pts.dat", fieldName);
            // std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
            cmd->SetUseDockerContainers(true);
            // convert to carp simple
            cmd->DockerConvertMeshFormat(directory, meshName, "vtk", meshName+"-temp", "carp_txt", 1);
            cmd->DockerConvertMeshFormat(directory, meshName+refinedSuffix, "vtk", meshName+refinedSuffix+"-temp", "carp_txt", 1);

            // new functionality: meshtool interpolate elemdata
            cmd->DockerInterpolateCell(directory, meshName+"-temp", meshName+refinedSuffix+"-temp", meshName+"_elem.dat", meshName+refinedSuffix+"_elem.dat");
            cmd->DockerInterpolateCell(directory, meshName+"-temp", meshName+refinedSuffix+"-temp", meshName+"_pts.dat", meshName+refinedSuffix+"_pts.dat");

            std::vector<double> refinedCellScalars = CemrgCommonUtils::ReadScalarField(Path(meshName+refinedSuffix+"_elem.dat"));
            std::vector<double> refinedPointScalars = CemrgCommonUtils::ReadScalarField(Path(meshName+refinedSuffix+"_pts.dat"));
            CemrgCommonUtils::AppendScalarFieldToVtk(refinedPath, fieldName, "POINT", refinedPointScalars, false);
            CemrgCommonUtils::AppendScalarFieldToVtk(refinedPath, fieldName, "CELL", refinedCellScalars, false);

            // cleanup
            QStringList exts = {".elem", ".pts", ".lon", ".fcon"};
            for (int ix = 0; ix < exts.size(); ix++) {
                QFile::remove(Path(meshName+"-temp"+exts.at(ix)));
                QFile::remove(Path(meshName+refinedSuffix+"-temp"+exts.at(ix)));
            }

        }

        atrium->SetSurface(refinedPath);

        if(uiRemesh_extractParts){
            MITK_INFO << "[MeshingOptions] save individual vtks for: RS, RI, LS, LI and LAA";
            atrium->ExtractLabelFromShell(directory, atrium->LSPV(), "LSPV");
            atrium->ExtractLabelFromShell(directory, atrium->LIPV(), "LIPV");
            atrium->ExtractLabelFromShell(directory, atrium->RSPV(), "RSPV");
            atrium->ExtractLabelFromShell(directory, atrium->RIPV(), "RIPV");
            atrium->ExtractLabelFromShell(directory, atrium->LAAP(), "LAAP");
        }

        if(uiRemesh_cleanmesh){
            QString cleanInName = meshName+refinedSuffix;
            QString cleanOutName = "clean-"+cleanInName;

            MITK_INFO(QFile::remove(Path(cleanInName+".fcon"))) << "Removed .fcon file";

            MITK_INFO << "[MeshingOptions] Cleaning up mesh quality";
            // std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
            cmd->SetUseDockerContainers(true);

            cmd->DockerCleanMeshQuality(directory, cleanInName, cleanOutName, 0.2, "vtk", "vtk_polydata");
            cmd->DockerCleanMeshQuality(directory, cleanOutName, cleanOutName, 0.1, "vtk", "vtk_polydata");

            if(!cmd->IsOutputSuccessful(Path(cleanOutName+".vtk"))){
                QMessageBox::warning(NULL, "Attention", "Clean output unsuccessful");
                MITK_ERROR << "Clean output unsuccessful";
                return;
            }
        }

        MITK_INFO << "[MeshingOptions] finished";
        MITK_INFO << "TIMELOG|MeshImprovement| End";
    }
}

bool  AtrialFibresView::UserSelectUacMesh(){
    if(uacMeshName.isEmpty()){
        QMessageBox::information(NULL, "Open Mesh File", "Open the CARP mesh file (pts ONLY)");
        QString meshPath = QFileDialog::getOpenFileName(NULL, "Open the CARP mesh file (pts ONLY)",
            StdStringPath().c_str(), tr("CARP points file (*.pts)"));
        QFileInfo fi(meshPath);

        uacMeshName = fi.baseName();
    }

    int reply = Ask("Check correct filename", ("Is file: " + uacMeshName + " correct?").toStdString());
    return (reply==QMessageBox::Yes);
}


void AtrialFibresView::UacCalculation(){
    MITK_INFO << "TIMELOG|UacCalculationSetup| Start";
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.
    if(!UserSelectUacMesh()) return;

    bool userInputAccepted = GetUserUacOptionsInputs();
    MITK_INFO(userInputAccepted) << "[UacCalculation] User Input accepted";

    if(userInputAccepted){
        SetFibresVariables();

        m_Controls.button_0_1_uacRough->setVisible(true);
        m_Controls.button_0_2_uacRefined->setVisible(true);
    }
    MITK_INFO << "TIMELOG|UacCalculationSetup| End";
}

void AtrialFibresView::UacCalculationRough(){
    MITK_INFO << "TIMELOG|UacCalculation_Stage1| Start";
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.
    if (!UserSelectUacMesh()) return;

    // at least in LA, the uac codes need the refined landmarks for both stages
    bool isRaLabelled = uac_whichAtrium.compare("RA", Qt::CaseInsensitive)==0;

    QString landmarkFilename = (isRaLabelled) ? "prodRaLandmarks" : "prodLaRefinedLandmarks";
    QString landmarkType = (isRaLabelled) ? "RA_LANDMARK" : "LA_REFINED";
    QString path2landmarks = LandmarkFilesCreated(landmarkFilename, landmarkType);

    if(path2landmarks.compare("FILE_NOT_FOUND")==0) return;

    QStringList landmarksFilesList;
    landmarksFilesList << path2landmarks;

    if(isRaLabelled){
        path2landmarks = LandmarkFilesCreated("prodRaRegion", "REGION");
        landmarksFilesList << path2landmarks;
    }

    bool userInputAccepted = GetUserUacOptionsInputs();
    MITK_INFO(userInputAccepted) << "[UacCalculation] User Input accepted";

    if(userInputAccepted){
        QString uacOutput;
        QStringList outputFiles;

        SetFibresVariables("UAC_1_"); // sets: uac_anatomy, uac_surftype, uac_whichAtrium, uac_fibre, uaccmd

        // Cemrg CMD
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->SetUseDockerContainers(true);
        MITK_INFO << "Do Rough UAC code from Docker";

        cmd->SetDockerImageUac();
        MITK_INFO << "TIMELOG|UacCalculation_Stage1| UAC 1 start";
        uacOutput = cmd->DockerUniversalAtrialCoordinates(directory, uaccmd, fibreAtlas, uacMeshName, uiLabels, landmarksFilesList);
        MITK_INFO << ("TIMELOG|UacCalculation_Stage1| UAC 1 end " + uacOutput).toStdString();

        outputFiles << "LSbc1.vtx" << "LSbc2.vtx";
        outputFiles << "PAbc1.vtx" << "PAbc2.vtx";

        if (!IsOutputFileCorrect(directory, outputFiles)){
            MITK_INFO << "TIMELOG|UacCalculation_Stage1| End (FAIL)";
            return;
        }

        MITK_INFO << "Create Laplace Solve files for LR and PA SOLVES";
        QString lr_par, pa_par;
        MITK_INFO << "TIMELOG|UacCalculation_Stage1| openCARP start";
        lr_par = CemrgCommonUtils::OpenCarpParamFileGenerator(directory, "carpf_laplace_LS.par", uacMeshName, "LSbc1", "LSbc2");
        pa_par = CemrgCommonUtils::OpenCarpParamFileGenerator(directory, "carpf_laplace_PA.par", uacMeshName, "PAbc1", "PAbc2");

        MITK_INFO << "Do Laplace Solves using Docker";
        cmd->SetDockerImageOpenCarp();
        QString lrLapSolve, paLapSolve;
        lrLapSolve = cmd->OpenCarpDocker(directory, lr_par, "LR_UAC_N2");
        paLapSolve = cmd->OpenCarpDocker(directory, pa_par, "PA_UAC_N2");
        MITK_INFO << ("TIMELOG|UacCalculation_Stage1| openCARP end " + lrLapSolve + " - " + paLapSolve).toStdString();

        bool uacOutputSuccess = IsOutputFileCorrect(directory, outputFiles);
        MITK_ERROR(!uacOutputSuccess) << ("Problem with " + uaccmd).toStdString();
        std::string msg = "UAC Calculation - Stage 1 ";
        msg += (uacOutputSuccess) ? "successful" : "failed";
        QMessageBox::information(NULL, "Attention", msg.c_str());
        MITK_INFO << "TIMELOG|UacCalculation_Stage1| End";
    }
}

void AtrialFibresView::UacCalculationRefined(){
    MITK_INFO << "TIMELOG|UacCalculation_Stage2| Start";
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.
    if (!UserSelectUacMesh()) return;

    bool isRaLabelled = uac_whichAtrium.compare("RA", Qt::CaseInsensitive)==0;
    QString landmarkFilename = (isRaLabelled) ? "prodRaLandmarks" : "prodLaRefinedLandmarks";
    QString landmarkType = (isRaLabelled) ? "RA_LANDMARK" : "LA_REFINED";
    QString path2landmarks = LandmarkFilesCreated(landmarkFilename, landmarkType);

    if(path2landmarks.compare("FILE_NOT_FOUND")==0){
        return;
    }

    QStringList landmarksFilesList;
    landmarksFilesList << path2landmarks;

    if(isRaLabelled){
        path2landmarks = LandmarkFilesCreated("prodRaRegion", "REGION");
        landmarksFilesList << path2landmarks;
    }

    bool userInputAccepted = GetUserUacOptionsInputs();
    MITK_INFO(userInputAccepted) << "[UacCalculation] User Input accepted";

    if(userInputAccepted){
        QString uacOutput;
        QStringList outputFiles;

        SetFibresVariables("UAC_2A_"); // sets: uac_anatomy, uac_surftype, uac_whichAtrium, uac_fibre, uaccmd

        // Cemrg CMD
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->SetUseDockerContainers(true);
        MITK_INFO << "Do Rough UAC code from Docker";

        outputFiles.clear();
        outputFiles << "AnteriorMesh.elem" << "PosteriorMesh.elem";
        outputFiles << "Ant_Strength_Test_PA1.vtx" << "Ant_Strength_Test_LS1.vtx";
        outputFiles << "Post_Strength_Test_PA1.vtx" << "Post_Strength_Test_LS1.vtx";
        MITK_INFO << "TIMELOG|UacCalculation_Stage2| UAC 2.1 - Start";
        cmd->SetDockerImageUac();
        uacOutput = cmd->DockerUniversalAtrialCoordinates(directory, uaccmd, fibreAtlas, uacMeshName, uiLabels, landmarksFilesList);
        MITK_INFO << ("TIMELOG|UacCalculation_Stage2| UAC 2.1 - End " + uacOutput).toStdString();

        if (!IsOutputFileCorrect(directory, outputFiles)){
            MITK_INFO << "TIMELOG|UacCalculation_Stage2| End (FAILED)";
            return;
        }

        QString lrp_par, udp_par, lra_par, uda_par;
        QString carpf_lr, carpf_ud;
        carpf_lr = "carpf_laplace_single_LR";
        carpf_ud = "carpf_laplace_single_UD";

        MITK_INFO << "TIMELOG|UacCalculation_Stage2| openCARP - Start";
        lrp_par = CemrgCommonUtils::OpenCarpParamFileGenerator(directory, carpf_lr+"_P.par", "PosteriorMesh", "", "Post_Strength_Test_LS1");
        udp_par = CemrgCommonUtils::OpenCarpParamFileGenerator(directory, carpf_ud+"_P.par", "PosteriorMesh", "", "Post_Strength_Test_PA1");
        lra_par = CemrgCommonUtils::OpenCarpParamFileGenerator(directory, carpf_lr+"_A.par", "AnteriorMesh", "", "Ant_Strength_Test_LS1");
        uda_par = CemrgCommonUtils::OpenCarpParamFileGenerator(directory, carpf_ud+"_A.par", "AnteriorMesh", "", "Ant_Strength_Test_PA1");

        cmd->SetDockerImageOpenCarp();

        QString lrpLapSolve, udpLapSolve, lraLapSolve, udaLapSolve;
        lrpLapSolve = cmd->OpenCarpDocker(directory, lrp_par, "LR_Post_UAC");
        udpLapSolve = cmd->OpenCarpDocker(directory, udp_par, "UD_Post_UAC");
        lraLapSolve = cmd->OpenCarpDocker(directory, lra_par, "LR_Ant_UAC");
        udaLapSolve = cmd->OpenCarpDocker(directory, uda_par, "UD_Ant_UAC");
        MITK_INFO << ("TIMELOG|UacCalculation_Stage2| openCARP - End " + lrpLapSolve + "-" + udpLapSolve + "-" + lraLapSolve + "-" + udaLapSolve).toStdString();

        SetFibresVariables("UAC_2B_");

        outputFiles.clear();
        outputFiles << "Labelled_Coords_2D_Rescaling_v3_C.vtk";
        outputFiles << "Labelled_Coords_2D_Rescaling_v3_C.elem";
        outputFiles << "Labelled_Coords_2D_Rescaling_v3_C.pts";
        MITK_INFO << "TIMELOG|UacCalculation_Stage2| UAC 2.2 - Start";
        cmd->SetDockerImageUac();
        uacOutput = cmd->DockerUniversalAtrialCoordinates(directory, uaccmd, fibreAtlas, uacMeshName, uiLabels, QStringList());
        MITK_INFO << "TIMELOG|UacCalculation_Stage2| UAC 2.2 - End";

        bool uacOutputSuccess = IsOutputFileCorrect(directory, outputFiles);
        MITK_ERROR(!uacOutputSuccess) << ("Problem with " + uaccmd).toStdString();
        std::string msg = "UAC Calculation - Stage 2 ";
        msg += (uacOutputSuccess) ? "successful" : "failed";
        QMessageBox::information(NULL, "Attention", msg.c_str());

        MITK_INFO << ("TIMELOG|UacCalculation_Stage2| End " + uacOutput).toStdString();
    }
}

void AtrialFibresView::UacFibreMapping(){
    MITK_INFO << "TIMELOG|UacFibreMapping| Start";
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.
    if (!UserSelectUacMesh()) return;

    QString metadata = Path("prodUacMetadata.txt");
    if(!QFile::exists(metadata)){
        MITK_INFO << "[UacFibreMapping] UAC Metadata not found. Use the UAC Calculation Button";

        std::string msg = "File prodUacMetadata.txt not found.\n";
        msg += "Use the Calculate UAC button.";

        QMessageBox::information(NULL, "File not found", msg.c_str());
        return;
    } else {
        MITK_INFO(GetUserUacOptionsInputs()) << "Reading UAC user inputs";
    }

    QString uacOutput;
    QStringList outputFiles;

    SetFibresVariables("UAC_FibreMapping"); // sets: uac_anatomy, uac_surftype, uac_whichAtrium, uac_fibre, uaccmd

    QMessageBox::information(NULL, "Attention", "Checking for UAC Calculation output");
    outputFiles << "Labelled_Coords_2D_Rescaling_v3_C.vtk";
    outputFiles << "Labelled_Coords_2D_Rescaling_v3_C.elem";
    outputFiles << "Labelled_Coords_2D_Rescaling_v3_C.pts";
    if (!IsOutputFileCorrect(directory, outputFiles)) return;

    std::cout << "[uac_fibreField]" << uac_fibreField.toStdString() << '\n';
    std::cout << "[output]" << uac_fibreFieldOutputName.toStdString() << '\n';

    outputFiles.clear();
    // outputFiles << "Fibre_1.vpts";

    // QStringList cmdargs;
    // cmdargs << (uac_fibreField+".lon");
    // if(uiUac_surftypeIndex==2){
    //     cmdargs << (uac_fibreField+".lon");
    //     uaccmd += "_Bilayer";
    // }
    // cmdargs << uac_fibreFieldOutputName;// output of fibremapping

    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    // cmd->SetDockerImageUac();
    // uacOutput = cmd->DockerUniversalAtrialCoordinates(directory, uaccmd, fibreAtlas, uacMeshName, cmdargs, QStringList(), "Fibre_1.vpts");

    QString _atrium = uiUac_whichAtrium.at(uiUac_whichAtriumIndex).toLower();
    QString _layer = uiUac_surftype.at(uiUac_surftypeIndex).toLower(); // check name
    QString _fibre = uiUac_fibreFile.at(uiUac_fibreFileIndex); // check name
    QString _omsh = "fibres_" + _fibre;

    QString alt_uac = cmd->DockerUacFibreMappingMode(directory, _atrium, _layer, _fibre, uacMeshName, false, _omsh);
    MITK_INFO << alt_uac.toStdString();
    bool uacOutputSuccess = cmd->IsOutputSuccessful(uacOutput);
    MITK_WARN(!uacOutputSuccess) << ("Not found " + uaccmd).toStdString();

    MITK_INFO << "Searching for .vpts files in directory";
    QStringList vptsfiles = GetFilesWithExtension(".vpts");
    std::string msg = "UAC Fibre Mapping \n Following paths found:\n\n";
    if(!vptsfiles.isEmpty()){
        for (int ix = 0; ix < vptsfiles.size(); ix++) {
            msg += (vptsfiles.at(ix).toStdString() + '\n');
        }
    }else{
        msg += "failed";
    }
    QMessageBox::information(NULL, "Fibre Mapping - finished", msg.c_str());


    MITK_INFO << ("Output path of fibres: " + uacOutput).toStdString();
    MITK_INFO << "Clearing auxiliary files";
    QStringList clearFiles;
    clearFiles << "A_Checker_LS.vtx"<< "A_Checker_PA.vtx"<< "Ant_Strength_Test_LS1.vtx";
    clearFiles << "Ant_Strength_Test_PA1.vtx"<< "BorderNodes.vtx"<< "LSbc1.vtx";
    clearFiles << "carpf_laplace_single_LR_A.par"<< "carpf_laplace_single_LR_P.par";
    clearFiles<< "carpf_laplace_single_UD_A.par"<< "carpf_laplace_single_UD_P.par"<< "P_Checker_LS.vtx";
    clearFiles<< "PAbc2.vtx"<< "PAbc1.vtx"<< "NewInferior_PA.vtx"<< "NewInferior_LS.vtx";
    clearFiles<< "LSPV_PA.vtx"<< "LSPV_LS.vtx"<< "LSbc2.vtx"<< "T2.vtk"<< "Test_Split.vtk";
    clearFiles<< "Test_post.vtk"<< "Test_ant.vtk"<< "RSPV_PA.vtx"<< "RSPV_LS.vtx";
    clearFiles<< "Test.vtk"<< "Test_sub.vtk"<< "Test_sub_rspv.vtk"<< "Test_sub_rspv_e.vtk";
    clearFiles << "Test_sub_ripv.vtk"<< "Test_sub_ripv_e.vtk"<< "Test_sub_lspv.vtk";
    clearFiles << "Test_sub_lspv_e.vtk"<< "Test_sub_lipv.vtk"<< "Test_sub_lipv_e.vtk";
    clearFiles << "Test_sub_laa.vtk"<< "Aux_2.elem"<< "Aux_2.lon"<< "Aux_2.pts";
    clearFiles << "carpf_laplace_LS.par"<< "carpf_laplace_PA.par"<< "P_Checker_PA.vtx";
    clearFiles << "Post_Strength_Test_LS1.vtx"<< "Post_Strength_Test_PA1.vtx";

    int delete_files_reply = Ask("Question", "Delete auxiliary and temporary files?");
    if(delete_files_reply==QMessageBox::Yes){
        int count=0;
        for (int idx = 0; idx < clearFiles.size(); idx++) {
            count += (QFile::remove(Path(clearFiles.at(idx)))) ? 1 : 0;
        }

        MITK_INFO(count==clearFiles.size()) << "All aux Files cleared successfully";
    }
    MITK_INFO << "TIMELOG|UacFibreMapping| End";
}

void AtrialFibresView::UacCalculationVerifyLabels(){
    MITK_INFO << "TIMELOG|VerifyLabels| Start";
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.
    if(GetUserUacOptionsInputs(false)){
        uac_whichAtrium = uiUac_whichAtrium.at(uiUac_whichAtriumIndex);
        MITK_INFO << ("[UacCalculationVerifyLabels] Seleted ["+uac_whichAtrium+"] analysis").toStdString();
    } else{
        MITK_INFO << "User cancelled selection of LA/RA selection";
        return;
    }
    if(!GetUserEditLabelsInputs()){
        MITK_INFO << "labels not checked. Stopping";
        return;
    }

    MITK_INFO(LoadSurfaceChecks()) << ("Loaded surface" + tagName).toStdString();
    std::string title, msg;
    mitk::Surface::Pointer surface = mitk::IOUtil::Load<mitk::Surface>(StdStringPath(tagName+".vtk"));
    std::vector<int> incorrectLabels;
    if(atrium->CheckLabelConnectivity(surface, uiLabels, incorrectLabels)){
        title = "Label connectivity verification - failure";
        msg = "Make sure label connectivity is correct.\n\n";
        msg += "Do you want to try automatic fixing? ";

        int reply_auto_fix = Ask(title, msg);
        if(reply_auto_fix == QMessageBox::Yes){
            MITK_INFO << "Solving labelling inconsistencies";
            for (unsigned int ix = 0; ix < incorrectLabels.size(); ix++) {
                atrium->FixSingleLabelConnectivityInSurface(surface, incorrectLabels.at(ix));
            }
            mitk::IOUtil::Save(surface, StdStringPath(tagName+".vtk"));

            title = "Labels fixed";
            msg = "Labels fixed - saved file: \n";
            msg += StdStringPath(tagName+".vtk");
            QMessageBox::information(NULL, title.c_str(), msg.c_str());

        } else{
            title = "Opening Mesh Preprocessing";
            msg = "Use the Fix labelling button to fix labelling connectivity";
            QMessageBox::information(NULL, title.c_str(), msg.c_str());

            //Show the plugin
            this->GetSite()->GetPage()->ResetPerspective();
            AtrialFibresClipperView::SetDirectoryFile(directory, tagName+".vtk", true);
            this->GetSite()->GetPage()->ShowView("org.mitk.views.atrialfibresclipperview");
        }

        return;
    } else{
        title = "Label connectivity verification - success";
        msg = "No connectivity problems found on labels";

        QMessageBox::information(NULL, title.c_str(), msg.c_str());
        MITK_INFO << msg;
    }
    MITK_INFO << "TIMELOG|VerifyLabels| End";
}

bool AtrialFibresView::IsOutputFileCorrect(QString dir, QStringList filenames){
    bool success = true;

    int countfails = 0;
    QString checkOutputMsg = "";

    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    for (int ix = 0; ix < filenames.size(); ix++) {
        bool okSingleTest = cmd->IsOutputSuccessful(dir + "/" + filenames.at(ix));

        if (!okSingleTest){
            MITK_ERROR << ("File(s) not created - " + filenames.at(ix)).toStdString();
            countfails++;
            checkOutputMsg += (countfails==1) ? "File(s) not created: " : "";
            checkOutputMsg += "\n " + filenames.at(ix);
        }
    }

    if (!checkOutputMsg.isEmpty()){
        std::string msg = checkOutputMsg.toStdString();
        QMessageBox::warning(NULL, "Warning", msg.c_str());
        success = false;
    }

    return success;

}

void AtrialFibresView::ConvertFormat(){
    if (!RequestProjectDirectoryFromUser()) return;

    QString meshPath = QFileDialog::getOpenFileName(NULL, "Open Mesh file",
        StdStringPath().c_str(), tr("Mesh (*.vtk)"));

    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);

    QFileInfo fi(meshPath);
    QString inext = (fi.suffix().contains(".elem")) ? "carp_txt" : fi.suffix();
    bool userInputsAccepted = GetUserConvertFormatInputs(fi.baseName(), inext);
    if(userInputsAccepted){
        MITK_INFO << "[ConvertFormat] Creating CARP output in microns";
        cmd->DockerConvertMeshFormat(directory, fi.baseName(), "vtk", uiFormat_outName, uiFormat_outExt, uiFormat_scale);
    }
}

void AtrialFibresView::ScarProjection(){

    if (!RequestProjectDirectoryFromUser()) return;
    QString prodPath = directory + "/";

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
    if(automaticPipeline){
        if(!QFile::copy(cleanPath, segPath)){
            ImageType::Pointer im = atrium->RemoveNoiseFromAutomaticSegmentation(directory);
            mitk::IOUtil::Save(mitk::ImportItkImage(im), cleanPath.toStdString());
            mitk::IOUtil::Save(mitk::ImportItkImage(im), segPath.toStdString());
        }
    } else{
        mitk::Image::Pointer segImageMitk = mitk::IOUtil::Load<mitk::Image>(StdStringPath("PVeinsLabelled.nii"));
        segImageMitk = CemrgCommonUtils::ReturnBinarised(segImageMitk);
        mitk::IOUtil::Save(segImageMitk, segPath.toStdString());
        mitk::IOUtil::Save(segImageMitk, cleanPath.toStdString());
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
        scarShell->GetVtkPolyData()->GetCellData()->GetScalars()->SetName("scar");

        QString scarPath = prodPath + "MaxScar.vtk";
        CemrgCommonUtils::SetCellDataToPointData(scarShell, scarPath, "scar");

        MITK_INFO << "[SCAR_PROJECTION][5] Thresholding";
        double mean = 0.0, stdv = 0.0;
        bool binarise=true;
        ImageType::Pointer segITK = atrium->LoadImage(segPath, binarise);
        mitk::IOUtil::Save(atrium->ImErode(segITK), Path("roi.nii").toStdString());

        atrium->ResampleSegmentationLabelToImage(Path("roi.nii"), lgePath);

        FloatImageType::Pointer roiFloat = FloatImageType::New();
        mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>(StdStringPath("roi.nii")), roiFloat);
        mitk::Image::Pointer roiImage = mitk::ImportItkImage(roiFloat);

        if(scar->CalculateMeanStd(lge, roiImage, mean, stdv)){
            MITK_INFO << "[SCAR_PROJECTION][6] Saving normalised shell and prodThresholds.txt";
            scar->SaveNormalisedScalars(mean, scarShell, (prodPath + "MaxScar_Normalised.vtk"));
            scar->PrintThresholdingResults(directory, uiScar_thresValues, uiScar_thresholdMethod, mean, stdv);
        } else{
            MITK_INFO << "[SCAR_PROJECTION][ERROR] Wrong dimensions between LGE and ROI Image";
        }

        QMessageBox::information(NULL, "Attention", "Scar projection finished!");

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

        if (succesfulAssignment){
            QString now = QDate::currentDate().toString(Qt::ISODate);
            QString logfilename = directory + "/afib_log" + now + ".log";
            std::string logfname_str = logfilename.toStdString();

            mitk::LoggingBackend::SetLogFile(logfname_str.c_str());
            MITK_INFO << ("Changed logfile location to: " + logfilename).toStdString();
        }

    } else {
        MITK_INFO << ("Project directory already set: " + directory).toStdString();
    }//_if


    return succesfulAssignment;
}

bool AtrialFibresView::GetUserConvertFormatInputs(QString inname, QString inext, bool cleanmesh){
    QDialog* inputs = new QDialog(0, Qt::WindowFlags());
    bool userInputAccepted=false;
    m_UIFormat.setupUi(inputs);
    connect(m_UIFormat.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_UIFormat.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    if(cleanmesh){
        m_UIFormat.titleLabel->setText("Clean mesh");
        m_UIFormat.lineEdit_2scale->setVisible(false);
    }

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

bool AtrialFibresView::GetUserUacOptionsInputs(bool enableFullUiOptions){
    QString metadata = Path("prodUacMetadata.txt");
    bool userInputAccepted=false;

    if(uiUac_fibreFile.size() == 0){
        uiUac_fibreFile << "1" << "2" << "3" << "4" << "5" << "6" << "7" << "A" << "L";
        uiUac_whichAtrium << "LA" << "RA";
        uiUac_surftype << "Endo" << "Epi" << "Bilayer";
    }

    if(enableFullUiOptions && QFile::exists(metadata)){
        int reply0 = Ask("Metadata Found", "Load previous analysis metadata found?");
        if(reply0==QMessageBox::Yes){
            std::ifstream fi(metadata.toStdString());
            fi >> uiUac_meshtype_labelled;
            fi >> uiUac_whichAtriumIndex;
            fi >> uiUac_fibreFileIndex;
            fi >> uiUac_surftypeIndex;
            fi.close();

            userInputAccepted=true;
        }
    }

    if(!userInputAccepted){
        QDialog* inputs = new QDialog(0, Qt::WindowFlags());
        m_UIUac.setupUi(inputs);
        connect(m_UIUac.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
        connect(m_UIUac.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));

        // enable or disable parts that might not be used
        m_UIUac.label_3->setVisible(enableFullUiOptions);
        m_UIUac.combo_uac_surftype->setVisible(enableFullUiOptions);
        m_UIUac.label->setVisible(enableFullUiOptions);
        m_UIUac.combo_uac_fibrefile->setVisible(enableFullUiOptions);
        m_UIUac.check_uac_meshtype_labelled->setVisible(enableFullUiOptions);

        int dialogCode = inputs->exec();

        //Act on dialog return code
        if (dialogCode == QDialog::Accepted) {
            userInputAccepted = true;
            uiUac_meshtype_labelled = m_UIUac.check_uac_meshtype_labelled->isChecked();
            uiUac_whichAtriumIndex = m_UIUac.combo_uac_whichAtrium->currentIndex();
            uiUac_fibreFileIndex = m_UIUac.combo_uac_fibrefile->currentIndex();
            uiUac_surftypeIndex = m_UIUac.combo_uac_surftype->currentIndex();

            if (enableFullUiOptions){
                std::ofstream fo(metadata.toStdString());
                fo << uiUac_meshtype_labelled << std::endl;
                fo << uiUac_whichAtriumIndex << std::endl;
                fo << uiUac_fibreFileIndex << std::endl;
                fo << uiUac_surftypeIndex << std::endl;
                fo.close();
            }

        } else if (dialogCode == QDialog::Rejected) {
            inputs->close();
            inputs->deleteLater();
        }//_if
    }

    return userInputAccepted;
}

bool AtrialFibresView::GetUserEditLabelsInputs(){
    bool userInputAccepted=false;

    if(!userInputAccepted){
        QDialog* inputs = new QDialog(0, Qt::WindowFlags());
        m_UIEditLabels.setupUi(inputs);
        connect(m_UIEditLabels.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
        connect(m_UIEditLabels.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
        std::cout << "WHICH ATRIUM: " << uac_whichAtrium.toStdString() << '\n';
        bool isLeftAtrium = (uac_whichAtrium.compare("LA", Qt::CaseInsensitive)==0);

        m_UIEditLabels.lineEdit_LA->setVisible(isLeftAtrium);
        m_UIEditLabels.lineEdit_LAA->setVisible(isLeftAtrium);
        m_UIEditLabels.lineEdit_LSPV->setVisible(isLeftAtrium);
        m_UIEditLabels.lineEdit_LIPV->setVisible(isLeftAtrium);
        m_UIEditLabels.lineEdit_RSPV->setVisible(isLeftAtrium);
        m_UIEditLabels.lineEdit_RIPV->setVisible(isLeftAtrium);
        m_UIEditLabels.lineEdit_RA->setVisible(!isLeftAtrium);
        m_UIEditLabels.lineEdit_RAA->setVisible(!isLeftAtrium);
        m_UIEditLabels.lineEdit_RA_SVC->setVisible(!isLeftAtrium);
        m_UIEditLabels.lineEdit_RA_IVC->setVisible(!isLeftAtrium);
        m_UIEditLabels.lineEdit_RA_CS->setVisible(!isLeftAtrium);

        int dialogCode = inputs->exec();
        //Act on dialog return code
        if (dialogCode == QDialog::Accepted) {
            userInputAccepted = true;
            bool ok1, ok2, ok3, ok4, ok5, ok6;
            if(isLeftAtrium){
                int la, laa, ls, li, rs, ri;

                la = m_UIEditLabels.lineEdit_LA->text().toInt(&ok1);
                laa = m_UIEditLabels.lineEdit_LAA->text().toInt(&ok2);
                ls = m_UIEditLabels.lineEdit_LSPV->text().toInt(&ok3);
                li = m_UIEditLabels.lineEdit_LIPV->text().toInt(&ok4);
                rs = m_UIEditLabels.lineEdit_RSPV->text().toInt(&ok5);
                ri = m_UIEditLabels.lineEdit_RIPV->text().toInt(&ok6);

                if (!ok1) la = 1;
                if (!ok2) laa = 19;
                if (!ok3) ls = 11;
                if (!ok4) li = 13;
                if (!ok5) rs = 15;
                if (!ok6) ri = 17;

                uiLabels.clear();
                uiLabels << QString::number(la);
                uiLabels << QString::number(laa);
                uiLabels << QString::number(ls);
                uiLabels << QString::number(li);
                uiLabels << QString::number(rs);
                uiLabels << QString::number(ri);
            } else{
                int ra, raa, svc, ivc, cs;

                ra = m_UIEditLabels.lineEdit_RA->text().toInt(&ok1);
                svc = m_UIEditLabels.lineEdit_RA_SVC->text().toInt(&ok3);
                ivc = m_UIEditLabels.lineEdit_RA_IVC->text().toInt(&ok4);
                cs = m_UIEditLabels.lineEdit_RA_CS->text().toInt(&ok5);
                raa = m_UIEditLabels.lineEdit_RAA->text().toInt(&ok2);

                if (!ok1) ra = 1;
                if (!ok3) svc = 6;
                if (!ok4) ivc = 7;
                if (!ok5) cs = 5;
                if (!ok2) raa = 2;

                uiLabels.clear();
                uiLabels << QString::number(ra);
                uiLabels << QString::number(svc);
                uiLabels << QString::number(ivc);
                uiLabels << QString::number(cs);
                uiLabels << QString::number(raa);
            }

        } else if (dialogCode == QDialog::Rejected) {
            inputs->close();
            inputs->deleteLater();
        }//_if
    }

    return userInputAccepted;
}

bool AtrialFibresView::GetUserRemeshingInputs(){
    QDialog* inputs = new QDialog(0, Qt::WindowFlags());
    bool userInputAccepted=false;
    m_UIRemesh.setupUi(inputs);
    connect(m_UIRemesh.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_UIRemesh.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    int dialogCode = inputs->exec();

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {
        uiRemesh_isscalar = m_UIRemesh.check_1isscalar->isChecked();
        uiRemesh_extractParts = m_UIRemesh.check_2extractlabels->isChecked();
        uiRemesh_cleanmesh = m_UIRemesh.check_3meshquality->isChecked();

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

bool AtrialFibresView::GetUserAnalysisSelectorInputs(){
    MITK_INFO << "[GetUserAnalysisSelectorInputs]";
    QString metadata = Path("prodMetadata.txt");
    bool userInputAccepted=false;

    if(QFile::exists(metadata)){
        int reply0 = Ask("Metadata Found", "Load previous analysis metadata found?");
        if(reply0==QMessageBox::Yes){
            std::ifstream fi(metadata.toStdString());
            fi >> uiSelector_pipeline;
            fi >> uiSelector_imgauto_skipCemrgNet;
            fi >> uiSelector_imgauto_skipLabel;
            fi >> uiSelector_man_useCemrgNet;
            fi >> uiSelector_img_scar;
            fi.close();

            userInputAccepted=true;

            MITK_INFO(LoadSurfaceChecks()) << ("Loaded surface" + tagName).toStdString();
        }
    }

    if(!userInputAccepted){
        QDialog* inputs = new QDialog(0, Qt::WindowFlags());
        m_UISelector.setupUi(inputs);
        connect(m_UISelector.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
        connect(m_UISelector.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
        int dialogCode = inputs->exec();

        //Act on dialog return code
        if (dialogCode == QDialog::Accepted) {
            uiSelector_img_scar = m_UISelector.check_img_scar->isChecked();

            if(m_UISelector.radioBtn_img_auto->isChecked()){
                uiSelector_pipeline = 0;
                uiSelector_imgauto_skipCemrgNet = m_UISelector.check_img_auto_skipSeg->isChecked();
                uiSelector_imgauto_skipLabel = m_UISelector.check_img_auto_skipLabel->isChecked();
            } else if(m_UISelector.radioBtn_img_man->isChecked()){
                uiSelector_pipeline = 1;
                uiSelector_man_useCemrgNet = m_UISelector.check_img_man_skipSeg->isChecked();
            } else{
                uiSelector_pipeline = 2;
                uiSelector_img_scar = false;
            }

            std::ofstream fo(metadata.toStdString());
            fo << uiSelector_pipeline << std::endl;
            fo << uiSelector_imgauto_skipCemrgNet << std::endl;
            fo << uiSelector_imgauto_skipLabel << std::endl;
            fo << uiSelector_man_useCemrgNet << std::endl;
            fo << uiSelector_img_scar << std::endl;
            fo.close();

            inputs->deleteLater();
            userInputAccepted=true;

        } else if (dialogCode == QDialog::Rejected) {
            inputs->close();
            inputs->deleteLater();
        }//_if
    }
    SetLgeAnalysis(uiSelector_img_scar);

    return userInputAccepted;
}

bool AtrialFibresView::GetUserMeshingInputs(){
    bool userInputAccepted=false;

    if(userHasSetMeshingParams){
        QString msg = "The parameters have been set already, change them?\n";
        msg += "threshold= " + QString::number(uiMesh_th) + '\n';
        msg += "blur= " + QString::number(uiMesh_bl) + '\n';
        msg += "smooth iter= " + QString::number(uiMesh_smth);

        if(Ask("Question", msg.toStdString())==QMessageBox::Yes){
            userHasSetMeshingParams=false;
            return GetUserMeshingInputs();
        } else{
            userInputAccepted=true;
        }
    } else{
        QDialog* inputs = new QDialog(0, Qt::WindowFlags());
        m_UIMeshing.setupUi(inputs);
        connect(m_UIMeshing.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
        connect(m_UIMeshing.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
        int dialogCode = inputs->exec();

        //Act on dialog return code
        if (dialogCode == QDialog::Accepted) {

            bool ok1, ok2, ok3;
            uiMesh_th= m_UIMeshing.lineEdit_1->text().toDouble(&ok1);
            uiMesh_bl= m_UIMeshing.lineEdit_2->text().toDouble(&ok2);
            uiMesh_smth= m_UIMeshing.lineEdit_3->text().toDouble(&ok3);

            //Set default values
            if (!ok1 || !ok2 || !ok3)
            QMessageBox::warning(NULL, "Attention", "Reverting to default parameters!");
            if (!ok1) uiMesh_th=0.5;
            if (!ok2) uiMesh_bl=0.0;
            if (!ok3) uiMesh_smth=10;

            inputs->deleteLater();
            userInputAccepted=true;
            userHasSetMeshingParams=true;

        } else if (dialogCode == QDialog::Rejected) {
            inputs->close();
            inputs->deleteLater();
        }//_if
    }

    return userInputAccepted;
}

bool AtrialFibresView::GetUserScarProjectionInputs(){
    bool userInputAccepted=false;

    QDialog* inputs = new QDialog(0, Qt::WindowFlags());
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
            separated_thresh_list << "0.97" << "1.2" << "1.32";
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
            #if QT_VERSION >= QT_VERSION_CHECK(5, 14, 0)
                separated_thresh_list = thresh_list.split(",", Qt::SkipEmptyParts);
            #else
                separated_thresh_list = thresh_list.split(",", QString::SkipEmptyParts);
            #endif
            int listspaces = separated_thresh_list.removeAll(" ");
            int listduplicates = separated_thresh_list.removeDuplicates();
            std::cout << "Spaces in list: " << listspaces << " Duplicates in list: " << listduplicates << '\n';
            separated_thresh_list.sort();
        }//_if


        bool vOK;
        for(int ix=0; ix<separated_thresh_list.size(); ix++) {
            MITK_INFO << separated_thresh_list.at(ix);
            double tryNumber = separated_thresh_list.at(ix).toDouble(&vOK);
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
    prodPath = directory + "/";
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

    if(b){
        m_Controls.button_auto4_meshpreproc->setText("    Step7: Mesh Preprocessing");
        m_Controls.button_0_landmarks->setText("    Step10: Select Landmarks");
        m_Controls.button_0_calculateUac->setText("    Step11: Calculate UAC");
        m_Controls.button_0_fibreMapUac->setText("    Step12: Fibre Mapping");
    } else{
        m_Controls.button_auto4_meshpreproc->setText("    Step4: Mesh Preprocessing");
    }

}

void AtrialFibresView::SetAutomaticModeButtons(bool b){
    m_Controls.button_auto4_meshpreproc->setVisible(b);
    m_Controls.button_auto5_clipPV->setVisible(b);

    if(b){
        m_Controls.button_0_landmarks->setText("    Step6: Select Landmarks");
        m_Controls.button_0_calculateUac->setText("    Step7: Calculate UAC");
        m_Controls.button_0_fibreMapUac->setText("    Step8: Fibre Mapping");
    }
}

void AtrialFibresView::SetTagNameFromPath(QString path){
    QFileInfo fi(path);
    tagName = fi.baseName();
}

void AtrialFibresView::SetFibresVariables(QString uaccmd_prefix){
    // Variables set in this function: uac_anatomy, uac_surftype, uac_whichAtrium, uac_fibre,
    // uaccmd is set through the uaccmd_prefix input
    //
    MITK_INFO << "[SetFibresVariables] Set variables based on UAC user-defined parameters";

    uac_anatomy = "6"; // might change later
    uac_surftype = (uiUac_surftypeIndex==2) ? "Endo" : uiUac_surftype.at(uiUac_surftypeIndex);
    uac_whichAtrium = uiUac_whichAtrium.at(uiUac_whichAtriumIndex);
    uac_fibre = uiUac_fibreFile.at(uiUac_fibreFileIndex);

    uac_fibreField = "Labelled_" + uac_anatomy + "_" + uac_fibre;
    uac_fibreFieldOutputName = "Fibre_" + uac_fibre;
    if(uiUac_fibreFileIndex==7){ // chosen Avg
        uac_fibreField = "Labelled_" + uac_fibre + "_" + uac_anatomy + "_1";
    }

    fibreAtlas.clear();
    fibreAtlas << "_" + uac_whichAtrium + "_" + uac_surftype;

    if(!uaccmd_prefix.isEmpty()){
        if(uaccmd_prefix.contains("FibreMapping")){
            uaccmd = "UAC_FibreMapping";
            if(uiUac_surftypeIndex==2){ // Check for Bilayer run
                fibreAtlas << ("_" + uac_whichAtrium + "_Epi");
            }
        } else{
            if(uiLabels.isEmpty()){
                uiLabels.clear();
                if(uiUac_meshtype_labelled){
                    QMessageBox::information(NULL, "Attention", "Check the labels are correct");
                    if(!GetUserEditLabelsInputs()){
                        MITK_INFO << "labels not checked. Stopping";
                        return;
                    }
                }
            }
            uaccmd = uaccmd_prefix + uac_whichAtrium;
            uaccmd += (!uiUac_meshtype_labelled) ? "_noPV" : "";
        }
    }

}

bool AtrialFibresView::LoadSurfaceChecks(){
    bool success = true;
    QString prodPath = directory + "/" + tagName + ".vtk";

    if(!QFile::exists(prodPath)){
        int reply1 = Ask("Surface file not found", "Do you have a surface file to load?");
        if(reply1==QMessageBox::Yes){
            UserLoadSurface();
        } else {
            success=false;
        }
    } else{
        std::string msg = ("Load automatically file called: [" + tagName + ".vtk]?").toStdString();
        int reply2 = Ask("Surface file to load", msg);
        if(reply2==QMessageBox::No){
            UserLoadSurface();
        }
    }

    if(success){
        prodPath = directory + "/" + tagName + ".vtk";
        MITK_INFO << ("[LoadSurfaceChecks] Loading surface: " + prodPath).toStdString();
    }

    return success;
}

void AtrialFibresView::UserLoadSurface(){
    QString newpath = QFileDialog::getOpenFileName(NULL, "Select the surface file to load!", directory.toStdString().c_str(), tr("Mesh (*.vtk)"));
    SetTagNameFromPath(newpath);
    CheckLoadedMeshQuality();
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

QStringList AtrialFibresView::GetFilesWithExtension(QString extension){
    QStringList result;
    QDirIterator searchit(directory, QDirIterator::Subdirectories);
    while(searchit.hasNext()) {
        QFileInfo searchfinfo(searchit.next());
        if (searchfinfo.fileName().contains(extension, Qt::CaseSensitive)) {
            result.push_back(searchfinfo.fileName());
        }
    }//_while
    return result;
}

QString AtrialFibresView::UserIncludeLgeAnalysis(QString segPath, ImageType::Pointer segImage){
    // uiSelector_img_scar
    QString prodPath = directory + "/";
    QString resultString = segPath;

    // int replyLAreg = Ask("Register to LGE", "Move to LGE space and prepare for scar projection?");
    if(uiSelector_img_scar){
        QFileInfo fi(segPath);
        QString segRegPath = prodPath + fi.baseName() + "-reg." + fi.suffix();

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
        cmd->ExecuteTransformation(directory, segPath, segRegPath);

        segImage = atrium->LoadImage(segRegPath);

        resultString = segRegPath;
        tagName += "-reg";

        SetLgeAnalysis(true);
    }

    return resultString;
}

void AtrialFibresView::SetLgeAnalysis(bool b){
    analysisOnLge = b;
    // m_Controls.button_z_scar->setEnabled(b);
}

void AtrialFibresView::CheckLoadedMeshQuality(){
    QString prodPath = directory + "/";
    QString meshinput = prodPath + tagName + ".vtk";

    MITK_INFO << "Safety check: converting file to polydata: ";
    std::string msg = "Is [" + (tagName + ".vtk").toStdString() + "] a vtk polydata?";
    int reply = Ask("Question", msg);
    if(reply==QMessageBox::No){
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->SetUseDockerContainers(true);
        meshinput = cmd->DockerConvertMeshFormat(directory, tagName, "vtk", tagName, "vtk_polydata", 1);
    }
    MITK_INFO << "Safety check: finished";

    mitk::Surface::Pointer surface = mitk::IOUtil::Load<mitk::Surface>(meshinput.toStdString());

    vtkSmartPointer<vtkPolyDataConnectivityFilter> cf = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    cf->ScalarConnectivityOff();
    cf->SetInputData(surface->GetVtkPolyData());
    cf->SetExtractionModeToAllRegions();
    cf->Update();

    int numRegions = cf->GetNumberOfExtractedRegions();

    if(numRegions>1){
        resurfaceMesh = true;
        MITK_INFO << ("Number of regions in mesh: " + QString::number(numRegions)).toStdString();
        cf->SetExtractionModeToSpecifiedRegions();
        cf->Modified();
        cf->Update();

        for (int ix = 0; ix < numRegions; ix++) {

            cf->InitializeSpecifiedRegionList();
            cf->AddSpecifiedRegion(ix);
            cf->Modified();
            cf->Update();

            vtkSmartPointer<vtkPolyData> temp = cf->GetOutput();
            surface->SetVtkPolyData(cf->GetOutput());
            QString nameExt = "connectivityTest_"+QString::number(ix)+".vtk";
            mitk::IOUtil::Save(surface, (prodPath+nameExt).toStdString());
        }
    }
}
