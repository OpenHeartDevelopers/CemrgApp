/*=========================================================================BSD 3-Clause License

Copyright (c) 2020, OpenHeartDevelopers
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=========================================================================*/

/*=========================================================================
*
* CemrgApp Plugin - Cemrg Atrial Strain Motion
*
* Cardiac Electromechanics Research Group
* http://www.cemrgapp.com
* info@cemrgapp.com 
*
* This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
* 
=========================================================================*/

#include "kcl_cemrgapp_atrialstrainmotion_Activator.h"
#include "AtrialStrainMotionView.h"
#include "AtrialFibresClipperView.h"

// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>

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

// Qt
#include <QMessageBox>
#include <QPushButton>
#include <QFileDialog>
#include <QInputDialog>
#include <QDirIterator>
#include <QDate>
#include <QFile>

// mitk
#include <QmitkIOUtil.h>
#include <mitkImage.h>
#include <mitkLog.h>
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

//CemrgAppModule
#include <CemrgCommonUtils.h>
#include <CemrgCommandLine.h>
#include <CemrgMeasure.h>
#include <CemrgScar3D.h>
#include "CemrgMultilabelSegmentationUtils.h"
#include <itkLabelStatisticsImageFilter.h>


const std::string AtrialStrainMotionView::VIEW_ID = "org.mitk.views.atrialstrainmotionview";

void AtrialStrainMotionView::SetFocus()
{
  m_Controls.segment_extract->setFocus();
}

void AtrialStrainMotionView::CreateQtPartControl(QWidget *parent)
{
  // create GUI widgets from the Qt Designer's .ui file
  m_Controls.setupUi(parent);
  // 1: Segment and Extract
  connect(m_Controls.segment_extract, SIGNAL(clicked()), this, SLOT(SegmentExtract()));
  m_Controls.segment_extract->setStyleSheet("Text-align:left");
  // 2. Post processing
  // m_Controls.button_man4_2_postproc->setVisible(false);
  connect(m_Controls.button_man4_2_postproc, SIGNAL(clicked()), this, SLOT(SegmentationPostprocessing()));
  m_Controls.button_man4_2_postproc->setStyleSheet("Text-align:left");
  // 3. Identify PVs
  connect(m_Controls.button_man5_idPV, SIGNAL(clicked()), this, SLOT(IdentifyPV()));
  m_Controls.button_man5_idPV->setStyleSheet("Text-align:left");
  // 4. Create Labelled Mesh
  connect(m_Controls.button_man6_labelmesh, SIGNAL(clicked()), this, SLOT(CreateLabelledMesh()));
  m_Controls.button_man6_labelmesh->setStyleSheet("Text-align:left");
  // 5. Mesh Preprocessing
  connect(m_Controls.button_auto4_meshpreproc, SIGNAL(clicked()), this, SLOT(MeshPreprocessing()));
  m_Controls.button_auto4_meshpreproc->setStyleSheet("Text-align:left");
  // 6. Verify Labels
  connect(m_Controls.button_0_3_checklabels, SIGNAL(clicked()), this, SLOT(UacCalculationVerifyLabels()));
  m_Controls.button_0_3_checklabels->setStyleSheet("Text-align:left");
  // 7. Clip PV
  connect(m_Controls.button_man8_clipPV, SIGNAL(clicked()), this, SLOT(ClipperPV()));
  m_Controls.button_man8_clipPV->setStyleSheet("Text-align:left");
  // 8. Mesh Improvement
  connect(m_Controls.btn_mesh_improvement, SIGNAL(clicked()), this, SLOT(MeshImprovement()));
  m_Controls.btn_mesh_improvement->setStyleSheet("Text-align:left");
  // 9. Auto Landmark
  connect(m_Controls.autoLM, SIGNAL(clicked()), this, SLOT(AutoLandMark()));
  m_Controls.autoLM->setStyleSheet("Text-align:left");
  // 10. UAC Stage 1
  connect(m_Controls.uac_stage_1, SIGNAL(clicked()), this, SLOT(UAC_Stage1()));
  m_Controls.uac_stage_1->setStyleSheet("Text-align:left");
  // 11. UAC Stage 2
  connect(m_Controls.uac_stage_2, SIGNAL(clicked()), this, SLOT(UAC_Stage2()));
  m_Controls.uac_stage_2->setStyleSheet("Text-align:left");
  // 12. Create Model
  connect(m_Controls.create_model, SIGNAL(clicked()), this, SLOT(CreateModel()));
  m_Controls.create_model->setStyleSheet("Text-align:left");
  // Crop Image
  connect(m_Controls.crop_images, SIGNAL(clicked()), this, SLOT(CropImage()));
  m_Controls.crop_images->setStyleSheet("Text-align:left");
  // Registration
  connect(m_Controls.registration, SIGNAL(clicked()), this, SLOT(Registration()));
  m_Controls.registration->setStyleSheet("Text-align:left");
  // Deform Mesh
  connect(m_Controls.deform_mesh, SIGNAL(clicked()), this, SLOT(DeformMesh()));
  m_Controls.deform_mesh->setStyleSheet("Text-align:left");
  // Generate Cell Area Strains
  connect(m_Controls.generate_cell_area_strains, SIGNAL(clicked()), this, SLOT(GenerateCellAreaStrains()));
  m_Controls.generate_cell_area_strains->setStyleSheet("Text-align:left");
  // Jacobian Threshold
  connect(m_Controls.jacobian_threshold, SIGNAL(clicked()), this, SLOT(JacobianThreshold()));
  m_Controls.jacobian_threshold->setStyleSheet("Text-align:left");
  // Plot Area Strain
  connect(m_Controls.plot_area_strain, SIGNAL(clicked()), this, SLOT(PlotAreaStrain()));
  m_Controls.plot_area_strain->setStyleSheet("Text-align:left");
  // Calculate Fiber Strain
  connect(m_Controls.calc_fiber_strains, SIGNAL(clicked()), this, SLOT(CalcFiberStrains()));
  m_Controls.calc_fiber_strains->setStyleSheet("Text-align:left");
  // Plot Fiber Strain
  connect(m_Controls.plot_fibers_trains, SIGNAL(clicked()), this, SLOT(PlotFibersTrains()));
  m_Controls.plot_fibers_trains->setStyleSheet("Text-align:left");


  // Set default variables
  atrium = std::unique_ptr<CemrgAtrialTools>(new CemrgAtrialTools());
  tagName = "Labelled";
  uiSelector_pipeline = 0;
  uiSelector_imgauto_skipCemrgNet = false;
  uiSelector_imgauto_skipLabel = false;
  uiSelector_img_scar = false;
  uiSelector_man_useCemrgNet = false;

}

void AtrialStrainMotionView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*source*/,
                                                const QList<mitk::DataNode::Pointer>& /*nodes*/)
{
  // iterate all selected objects, adjust warning visibility
  /*
  foreach (mitk::DataNode::Pointer node, nodes)
  {
    if (node.IsNotNull() && dynamic_cast<mitk::Image *>(node->GetData()))
    {
      m_Controls.labelWarning->setVisible(false);
      m_Controls.segment_extract->setEnabled(true);
      return;
    }
  }

  m_Controls.labelWarning->setVisible(true);
  m_Controls.segment_extract->setEnabled(true);

   */
}

bool AtrialStrainMotionView::GetUserAnalysisSelectorInputs(){
    MITK_INFO << "[GetUserAnalysisSelectorInputs]";
    QString metadata = Path("UAC_CT/prodMetadata.txt");
    bool userInputAccepted=false;

    if(QFile::exists(metadata)){
        std::ifstream fi(metadata.toStdString());
        fi >> uiSelector_pipeline;
        fi >> uiSelector_imgauto_skipCemrgNet;
        fi >> uiSelector_imgauto_skipLabel;
        fi >> uiSelector_man_useCemrgNet;
        fi >> uiSelector_img_scar;
        fi.close();

        userInputAccepted=true;

    }

    if(!userInputAccepted){
        uiSelector_img_scar = 0;
        uiSelector_pipeline = 2;
        uiSelector_img_scar = false;

        std::ofstream fo(metadata.toStdString());
        fo << uiSelector_pipeline << std::endl;
        fo << uiSelector_imgauto_skipCemrgNet << std::endl;
        fo << uiSelector_imgauto_skipLabel << std::endl;
        fo << uiSelector_man_useCemrgNet << std::endl;
        fo << uiSelector_img_scar << std::endl;
        fo.close();

        userInputAccepted=true;

    }
    SetLgeAnalysis(uiSelector_img_scar);

    return userInputAccepted;
}

int AtrialStrainMotionView::Ask(std::string title, std::string msg){
    return QMessageBox::question(NULL, title.c_str(), msg.c_str(), QMessageBox::Yes, QMessageBox::No);
}

bool AtrialStrainMotionView::LoadSurfaceChecks(){
	QString prodPath = directory + "/UAC_CT/" + tagName + ".vtk";
    MITK_INFO << ("[LoadSurfaceChecks] Loading surface: " + prodPath).toStdString();

    // A failed step still leaves a 0-byte mesh behind, so check the size as well as existence.
    QFileInfo meshInfo(prodPath);
    if (!meshInfo.exists() || meshInfo.size() == 0) {
        std::string msg = "The surface mesh is missing or empty:\n  " + prodPath.toStdString() +
            "\n\nAn earlier step did not write it. Look in the log for the command that failed. "
            "Then run that step again.";
        MITK_ERROR << msg;
        QMessageBox::warning(NULL, "Surface mesh not available", msg.c_str());
        return false;
    }

    return UserLoadSurface();
}

void AtrialStrainMotionView::SetLgeAnalysis(bool b){
    analysisOnLge = b;
    // m_Controls.button_z_scar->setEnabled(b);
}

bool AtrialStrainMotionView::CheckLoadedMeshQuality(){
    QString prodPath = directory + "/UAC_CT/";
    QString meshinput = prodPath + tagName + ".vtk";

    mitk::Surface::Pointer surface;
    try {
        surface = mitk::IOUtil::Load<mitk::Surface>(meshinput.toStdString());
    } catch (const std::exception& e) {
        std::string msg = "The application cannot read the surface mesh:\n  " + meshinput.toStdString() +
            "\n\nThe file exists, but it is not a valid VTK surface. The step that wrote it "
            "stopped before it finished. The reader reports:\n  " + e.what() +
            "\n\nRun that step again. Look in the log for the command that failed.";
        MITK_ERROR << msg;
        QMessageBox::critical(NULL, "Could not read the surface mesh", msg.c_str());
        return false;
    }

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

    return true;
}

void AtrialStrainMotionView::SetAutomaticModeButtons(bool b){
    m_Controls.button_auto4_meshpreproc->setVisible(b);
    // m_Controls.button_auto5_clipPV->setVisible(b);

    if(b){
        // m_Controls.button_0_landmarks->setText("    Step6: Select Landmarks");
        // m_Controls.button_0_calculateUac->setText("    Step7: Calculate UAC");
        // m_Controls.button_0_fibreMapUac->setText("    Step8: Fibre Mapping");
    }
}

bool AtrialStrainMotionView::UserLoadSurface(){
    QString newpath = directory + "/UAC_CT/" + tagName + ".vtk";
    SetTagNameFromPath(newpath);
    return CheckLoadedMeshQuality();
}

QString AtrialStrainMotionView::GetFilePath(QString nameSubstring, QString extension){
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

void AtrialStrainMotionView::SetTagNameFromPath(QString path){
    QFileInfo fi(path);
    tagName = fi.baseName();
}

void AtrialStrainMotionView::AutomaticAnalysis(){
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

        cmd->ExecuteSurf(directory, Path("prodClean.nii"), "close", uiMesh_iter, uiMesh_th, uiMesh_bl, uiMesh_smth);
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

QString AtrialStrainMotionView::UserIncludeLgeAnalysis(QString segPath, ImageType::Pointer segImage){
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

bool AtrialStrainMotionView::GetUserMeshingInputs(){
    bool userInputAccepted=false;

    uiMesh_th=0.5;
    uiMesh_bl=0.0;
    uiMesh_smth=10;
    uiMesh_iter=0;

    userInputAccepted=true;
    userHasSetMeshingParams=true;

    return userInputAccepted;
}

bool AtrialStrainMotionView::GetUserUacOptionsInputs(bool enableFullUiOptions){
    QString metadata = Path("UAC_CT/prodUacMetadata.txt");
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

bool AtrialStrainMotionView::GetUserEditLabelsInputs(){
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

/**
 * @brief Describe why an output file is unusable. Return an empty string if the file is usable.
 *
 * CemrgCommandLine::IsOutputSuccessful tests only that a file exists and is not empty.
 * ExecuteCommand also pre-creates each output with ExecuteTouch. A failed tool thus leaves a
 * file that passes these two criteria but may not be usuable. This function reads the file and names the fault.
 */
static QString DescribeOutputProblem(QString fullPath) {
    QFileInfo info(fullPath);

    if (!info.exists())
        return "does not exist";

    if (info.size() == 0)
        return "is empty. The tool started, but it wrote no data";

    QString name = info.fileName();
    if (name.endsWith(".nii", Qt::CaseInsensitive) || name.endsWith(".nii.gz", Qt::CaseInsensitive)) {
        QString message;
        CemrgCommonUtils::NiftiQformStatus status = CemrgCommonUtils::InspectNiftiQform(fullPath, message);
        if (status == CemrgCommonUtils::NiftiQformStatus::NotNifti ||
            status == CemrgCommonUtils::NiftiQformStatus::IoError)
            return "is not a valid NIfTI image (" + message + ")";
        if (status == CemrgCommonUtils::NiftiQformStatus::CannotRepair)
            return "has a NIfTI header the application cannot use (" + message + ")";
    }

    if (name.endsWith(".vtk", Qt::CaseInsensitive)) {
        QFile file(fullPath);
        if (!file.open(QIODevice::ReadOnly))
            return "is not readable. Check the file permissions";
        QByteArray head = file.read(16);
        file.close();
        if (!head.startsWith("# vtk") && !head.startsWith("<?xml") && !head.startsWith("<VTKFile"))
            return "does not start with a VTK header, so it is not a valid VTK file";
    }

    return QString();
}

// The label that the CCTA container gives to the LA body. SegmentExtract changes labels 8, 9
// and 10 to this value. It changes all other labels to 0. A complete LA.nii holds only this value.
static const int LA_LABEL = 4;

/**
 * @brief Describe why the application cannot reuse an output file. Return an empty string if the
 *        application can reuse the file.
 *
 * This test is stricter than DescribeOutputProblem. That function tests a file that the
 * application repairs and then loads. Reuse does not change the file, so this function tests the
 * file in its current state.
 *
 * @param copiedPaths  Files that the pipeline copies from the source data. The pipeline does not
 *                     change their headers. They must only be readable.
 * @param writtenPaths Files that the pipeline writes. The pipeline repairs the qform of each file,
 *                     so a complete run leaves a valid qform. A repairable qform shows an
 *                     incomplete run.
 */
static QString DescribeReusableOutputProblem(const QStringList& copiedPaths, const QStringList& writtenPaths) {
    QStringList allPaths = copiedPaths + writtenPaths;

    for (int ix = 0; ix < allPaths.size(); ix++) {
        QString fullPath = allPaths.at(ix);
        QString name = QFileInfo(fullPath).fileName();

        QString problem = DescribeOutputProblem(fullPath);
        if (!problem.isEmpty())
            return name + " " + problem;

        bool isNifti = name.endsWith(".nii", Qt::CaseInsensitive) ||
                       name.endsWith(".nii.gz", Qt::CaseInsensitive);
        if (isNifti && writtenPaths.contains(fullPath)) {
            QString message;
            if (CemrgCommonUtils::InspectNiftiQform(fullPath, message) !=
                CemrgCommonUtils::NiftiQformStatus::AlreadyValid)
                return name + " does not already declare a usable qform (" + message + ")";
        }
    }

    return QString();
}

/**
 * @brief Check that LA.nii holds only LA_LABEL. Return an empty string if it does.
 *
 * The file checks cannot find this fault. A valid NIfTI file can still hold the wrong image, for
 * example the multilabel segmentation with the wrong name. Reuse gives this file to all later
 * steps, so the application reads the volume one time before it offers reuse.
 */
static QString DescribeLaLabelProblem(QString laPath) {
    QString name = QFileInfo(laPath).fileName();

    mitk::Image::Pointer la;
    try {
        la = mitk::IOUtil::Load<mitk::Image>(laPath.toStdString());
    } catch (const std::exception& e) {
        return name + " is not readable (" + QString(e.what()) + ")";
    }

    std::vector<int> labels;
    std::unique_ptr<CemrgMultilabelSegmentationUtils> cmls(new CemrgMultilabelSegmentationUtils());
    cmls->GetLabels(la, labels); // background is excluded

    QStringList unexpected;
    bool foundLa = false;
    for (size_t ix = 0; ix < labels.size(); ix++) {
        if (labels.at(ix) == LA_LABEL)
            foundLa = true;
        else
            unexpected << QString::number(labels.at(ix));
    }

    if (!foundLa) {
        QString found = unexpected.isEmpty() ? "nothing but background" : ("only " + unexpected.join(", "));
        return name + " contains no LA label (" + QString::number(LA_LABEL) + ") - it holds " + found;
    }

    if (!unexpected.isEmpty())
        return name + " still contains label(s) " + unexpected.join(", ") +
            ". Step 1 removes these labels, so its label arithmetic did not finish";

    return QString();
}

bool AtrialStrainMotionView::IsOutputFileCorrect(QString dir, QStringList filenames){
    QStringList problems;

    for (int ix = 0; ix < filenames.size(); ix++) {
        QString filename = filenames.at(ix);
        QString problem = DescribeOutputProblem(dir + "/" + filename);

        if (!problem.isEmpty()) {
            MITK_ERROR << ("Expected output " + filename + " " + problem).toStdString();
            problems << "  " + filename + "\n      " + problem;
        }
    }

    if (problems.isEmpty())
        return true;

    std::string msg = "The application cannot use these output files:\n\n" +
        problems.join("\n\n").toStdString() +
        "\n\nLook in the log for the command that the step ran. Check that the Docker image "
        "is available. Check that the earlier steps completed.";
    QMessageBox::warning(NULL, "Step did not produce its expected output", msg.c_str());
    return false;
}

bool AtrialStrainMotionView::RequestProjectDirectoryFromUser() {

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

void AtrialStrainMotionView::SegmentExtract() {
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.

    QDir q_directory(directory);
    q_directory.mkdir("UAC_CT") ? MITK_INFO << "UAC_CT created" : MITK_INFO << "UAC_CT exists";
    q_directory.mkdir("UAC_CT_aligned") ? MITK_INFO << "UAC_CT_aligned created" : MITK_INFO << "UAC_CT_aligned exists";
    q_directory.mkdir("tracking") ? MITK_INFO << "tracking folder created" : MITK_INFO << "tracking folder exists";

    // copy the first image and segment it
    QDir nifti_dir = QDir(directory + "/nifti/");
    QStringList nifti_entries = nifti_dir.entryList(QDir::AllEntries | QDir::NoDotAndDotDot);
    if (nifti_entries.isEmpty()) {
        std::string msg = "Step 1 found no images in:\n  " + nifti_dir.absolutePath().toStdString() +
            "\n\nStep 1 expects the cine frames inside a 'nifti' folder in the project directory. "
            "Put the frames there. Then run step 1 again.";
        MITK_ERROR << msg;
        QMessageBox::warning(NULL, "No input images found", msg.c_str());
        return;
    }
    QString input_file_name = nifti_entries.at(0);
    QString input_file_path = directory + "/" + input_file_name;
    MITK_INFO << input_file_path.toStdString();

    // DockerCctaMultilabelSegmentation derives its output name from the input's base name.
    QString pathToSegmentation = directory + "/" + QFileInfo(input_file_path).baseName() + "_label_maps.nii.gz";

    QString la_path = directory + "/LA.nii";
    QString la_msh_path = directory + "/LA_msh.vtk";

    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);

    // Segmenting takes a couple of minutes, so offer to keep a previous result rather than silently
    // repeating it. Only offer that when the whole of the previous run can be trusted: these four
    // files are exactly what reuse would leave untouched, so any doubt about one of them means the
    // step runs in full, with no prompt. The UAC_CT copies and prodMetadata.txt are deliberately
    // absent - both paths rebuild those below, so their state cannot make reuse unsafe.
    QStringList copiedOutputs;
    copiedOutputs << input_file_path;
    QStringList writtenOutputs;
    writtenOutputs << pathToSegmentation << la_path << la_msh_path;

    MITK_INFO << "[SegmentExtract] Checking whether the previous run's output can be reused.";
    QString reuseProblem = DescribeReusableOutputProblem(copiedOutputs, writtenOutputs);
    if (reuseProblem.isEmpty())
        reuseProblem = DescribeLaLabelProblem(la_path);

    bool reuseExistingSegmentation = false;
    if (reuseProblem.isEmpty()) {
        QStringList names;
        names << QFileInfo(input_file_path).fileName()
              << QFileInfo(pathToSegmentation).fileName()
              << QFileInfo(la_path).fileName()
              << QFileInfo(la_msh_path).fileName();

        QStringList descriptions;
        descriptions << "the input image"
                     << "the segmentation"
                     << "the extracted left atrium"
                     << "the extracted surface";

        // The file names come from the project, so the column widths have to be measured rather
        // than hardcoded.
        int nameWidth = 0;
        int descriptionWidth = 0;
        for (int ix = 0; ix < names.size(); ix++) {
            nameWidth = qMax(nameWidth, names.at(ix).length());
            descriptionWidth = qMax(descriptionWidth, descriptions.at(ix).length());
        }

        QString table;
        for (int ix = 0; ix < names.size(); ix++) {
            table += names.at(ix).leftJustified(nameWidth + 3) +
                     descriptions.at(ix).leftJustified(descriptionWidth + 3) + "OK\n";
        }

        // A message box uses a proportional font, so the columns only line up inside <pre>.
        QMessageBox box(QMessageBox::Question, "Re-run the segmentation step?", "");
        box.setTextFormat(Qt::RichText);
        box.setText("Output files for this step already exist. Re-run the segmentation step?<pre>" +
                    table.toHtmlEscaped() + "</pre>");

        QPushButton* rerunButton = box.addButton("Yes, re-run", QMessageBox::YesRole);
        QPushButton* keepButton = box.addButton("No, use existing files", QMessageBox::NoRole);
        box.setDefaultButton(rerunButton);
        box.setEscapeButton(keepButton); // Escape should take the action that destroys nothing
        box.exec();

        reuseExistingSegmentation = (box.clickedButton() == keepButton);
        MITK_INFO << (reuseExistingSegmentation
            ? "[SegmentExtract] Keeping the existing segmentation and surface at the user's request."
            : "[SegmentExtract] Re-running the segmentation at the user's request.");
    } else {
        // Say which file failed, so the log explains why no prompt appeared.
        MITK_INFO << ("[SegmentExtract] Running the step in full - " + reuseProblem).toStdString();
    }

    if (!reuseExistingSegmentation) {
        // QFile::copy will not overwrite, so clear any previous working copy first. Without this
        // the copy fails silently and the pipeline continues against a stale image.
        if (QFile::exists(input_file_path) && !QFile::remove(input_file_path)) {
            std::string msg = "The application cannot replace the working copy:\n  " + input_file_path.toStdString() +
                "\n\nAnother program holds the file open, or the file is read-only. "
                "Close the other program. Then run step 1 again.";
            MITK_ERROR << msg;
            QMessageBox::warning(NULL, "Could not replace the working copy", msg.c_str());
            return;
        }
        if (!QFile::copy(nifti_dir.absolutePath() + "/" + input_file_name, input_file_path)) {
            std::string msg = "The application cannot copy " + input_file_name.toStdString() +
                " into the project directory:\n  " + directory.toStdString() +
                "\n\nCheck that the disk has free space. Check that you can write to the folder.";
            MITK_ERROR << msg;
            QMessageBox::warning(NULL, "Could not copy the input image", msg.c_str());
            return;
        }
        MITK_INFO << ("Copied " + input_file_name).toStdString();

        // extract LA
        QString producedSegmentation = cmd->DockerCctaMultilabelSegmentation(directory, input_file_path, false);
        if (producedSegmentation.isEmpty()) {
            std::string msg = "The multilabel segmentation did not produce an output file.\n\n"
                "Check that Docker runs. Check that the cemrg/ccta image is available. "
                "Look in the log for the command that step 1 ran.";
            MITK_ERROR << msg;
            QMessageBox::warning(NULL, "Segmentation failed", msg.c_str());
            return;
        }
        pathToSegmentation = producedSegmentation;
        MITK_INFO << pathToSegmentation;

        // cmd->DockerAtrialStrainMotion(directory, "fix_vtk_metadata");

        // The segmentation container rebuilds the NIfTI header from the image affine and leaves
        // qform_code at 0, which ITK refuses to read. Repair it before loading. Only files this
        // pipeline created are touched; the originals under nifti/ are left alone.
        QString qformMessage;
        CemrgCommonUtils::NiftiQformStatus qformStatus =
            CemrgCommonUtils::RepairNiftiQform(pathToSegmentation, qformMessage);
        switch (qformStatus) {
            case CemrgCommonUtils::NiftiQformStatus::Repaired:
                MITK_INFO << ("[SegmentExtract] NIfTI header repaired - " + qformMessage).toStdString();
                break;
            case CemrgCommonUtils::NiftiQformStatus::AlreadyValid:
                MITK_INFO << ("[SegmentExtract] NIfTI header needs no repair - " + qformMessage).toStdString();
                break;
            default:
                MITK_WARN << ("[SegmentExtract] NIfTI header could not be repaired (" +
                    CemrgCommonUtils::NiftiQformStatusToString(qformStatus) + ") - " + qformMessage).toStdString();
                break;
        }

        mitk::Image::Pointer segmentation;
        try {
            segmentation = mitk::IOUtil::Load<mitk::Image>(pathToSegmentation.toStdString());
        } catch (const std::exception& e) {
            std::string msg = "The application cannot read the segmentation from the CCTA container:\n  " +
                pathToSegmentation.toStdString() +
                "\n\nAn 'orthonormal direction cosines' error means the file declares no usable "
                "spatial transform. The reader reports:\n  " + e.what() +
                "\n\nLook in the log for the result of the header check.";
            MITK_ERROR << msg;
            QMessageBox::critical(NULL, "Could not read the segmentation", msg.c_str());
            return;
        }

        // create CemrgMultilabelSegmntation (cmls)
        // Keep this in step with DescribeLaLabelProblem, which reads LA.nii back and expects
        // LA_LABEL to be the only label left.
        std::unique_ptr<CemrgMultilabelSegmentationUtils> cmls(new CemrgMultilabelSegmentationUtils());
        segmentation = cmls->ReplaceLabel(segmentation, 8, LA_LABEL);
        segmentation = cmls->ReplaceLabel(segmentation, 9, LA_LABEL);
        segmentation = cmls->ReplaceLabel(segmentation, 10, LA_LABEL);
        for (int i = 1; i <= 7; i++) {
            if (i != LA_LABEL)
                segmentation = cmls->ReplaceLabel(segmentation, i, 0);
        }

        mitk::IOUtil::Save(segmentation, la_path.toStdString());

        // CemrgCommonUtils::AddToStorage(segmentation, "LA", this->GetDataStorage());

        // Extract LA Mesh surface
        // cmd->ExecuteSurf(directory, directory + "/LA.nii", "close", 0, 0.5, 0, 100);
        QString surfaceResult = cmd->ExecuteExtractSurface(directory, la_path, la_msh_path, 0.5, 0);
        // cmd->ExecuteSmoothSurface(directory, la_msh_path, directory + "/LA_msh_smth.vtk", 100);

        // ExecuteCommand pre-creates the output with ExecuteTouch, so a failed run still leaves an
        // empty file behind. Check the size, not just existence, or the empty mesh propagates and
        // surfaces much later as an opaque read error.
        if (surfaceResult.compare("ERROR_IN_PROCESSING") == 0 || QFileInfo(la_msh_path).size() == 0) {
            std::string msg = "Step 1 cannot extract the LA surface from " + QFileInfo(la_path).fileName().toStdString() +
                ".\n\nThis step runs MIRTK's 'extract-surface'. CemrgApp expects that program in:\n  " +
                (QCoreApplication::applicationDirPath() + "/MLib").toStdString() +
                "\n\nCemrgApp does not build the MIRTK binaries. Install them in that folder. "
                "Then check that the programs run. An absent shared library is a common cause.";
            MITK_ERROR << msg;
            QMessageBox::critical(NULL, "Surface extraction failed", msg.c_str());
            return;
        }
    }

    // mitk::Surface::Pointer la_msh = mitk::IOUtil::Load<mitk::Surface>(la_msh_path.toStdString());
    // CemrgCommonUtils::AddToStorage(la_msh, "la_msh", this->GetDataStorage());

    QString la_msh_uac_path = directory + "/UAC_CT" + "/LA_msh.vtk";
    if (QFile::exists(la_msh_uac_path))
        QFile::remove(la_msh_uac_path);
    if (!QFile::copy(la_msh_path, la_msh_uac_path)) {
        std::string msg = "The application cannot copy LA_msh.vtk into the UAC_CT folder:\n  " +
            (directory + "/UAC_CT").toStdString() +
            "\n\nCheck that the disk has free space. Check that you can write to the folder.";
        MITK_ERROR << msg;
        QMessageBox::warning(NULL, "Could not copy the surface mesh", msg.c_str());
        return;
    }
    MITK_INFO << "Copied LA_msh.vtk";

    // Analysis Selector
    bool userInputsAccepted = GetUserAnalysisSelectorInputs();
    uiSelector_pipeline = 2;
    if(userInputsAccepted){
        // uiSelector_pipeline; // =0 (imgAuto), =1 (imgManual), =2 (surf)
        SetAutomaticPipeline(uiSelector_pipeline==0); // Image Pipeline (automatic)
        MITK_INFO << "uiSelector_pipeline " << uiSelector_pipeline;

        MITK_INFO << "[AnalysisChoice] Analysis starting from surface";
        SetAutomaticPipeline(false);

        // Load surface mesh: LA_msh.vtk
		tagName = "LA_msh";
        if(!LoadSurfaceChecks()) return;

        // Create fake segmentation image for labelling
        double origin[3] = {0, 0, 0};
        double spacing[3] = {1, 1, 1};
		MITK_INFO << tagName;
        CemrgCommonUtils::SaveImageFromSurfaceMesh(Path("UAC_CT/" + tagName + ".vtk"), origin, spacing);
        CemrgCommonUtils::SavePadImageWithConstant(Path("UAC_CT/" + tagName + ".nii"));

        mitk::Image::Pointer im = CemrgCommonUtils::ReturnBinarised(mitk::IOUtil::Load<mitk::Image>(StdStringPath("UAC_CT/" + tagName+".nii")));

        mitk::IOUtil::Save(im, StdStringPath("UAC_CT/" + tagName+".nii"));
        CemrgCommonUtils::AddToStorage(im, tagName.toStdString(), this->GetDataStorage());

    }

}

void AtrialStrainMotionView::SegmentationPostprocessing(){
    this->GetSite()->GetPage()->ShowView("org.mitk.views.segmentationutilities");
}


void AtrialStrainMotionView::IdentifyPV(){
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.
    QString path = Path("UAC_CT/") + "segmentation.vtk";

	userHasSetMeshingParams = false;

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
                    mitk::IOUtil::Save(image, StdStringPath("UAC_CT/prodClean.nii"));

                    MITK_INFO << "[IdentifyPV] Create surface file and projecting tags";
                    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
                    cmd->SetUseDockerContainers(true);

                    cmd->ExecuteSurf(Path("UAC_CT"), Path("UAC_CT/prodClean.nii"), "close", uiMesh_iter, uiMesh_th, uiMesh_bl, uiMesh_smth);

                    //Add the mesh to storage
                    std::string meshName = segNode->GetName() + "-Mesh";
                    CemrgCommonUtils::AddToStorage(
                        CemrgCommonUtils::LoadVTKMesh(path.toStdString()), meshName, this->GetDataStorage());
                }
            }
        }
    }

    MITK_INFO(QFile::copy(Path("UAC_CT/segmentation.vtk"), Path("UAC_CT_aligned/segmentation.vtk"))) << "Copied: " << "segmentation.vtk";

    MITK_INFO << "Loading org.mitk.views.atrialfibresclipperview";
    this->GetSite()->GetPage()->ResetPerspective();
    this->GetSite()->GetPage()->ShowView("org.mitk.views.atrialfibresclipperview");
}

void AtrialStrainMotionView::CreateLabelledMesh() {
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

    QString prodPath =  directory + "/UAC_CT/";

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
        mitk::IOUtil::Save(segIm, StdStringPath("UAC_CT/prodClean.nii"));

        MITK_INFO << "[CreateLabelledMesh] Create surface file and projecting tags";
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->SetUseDockerContainers(true);

        cmd->ExecuteSurf(directory + "/UAC_CT", Path("UAC_CT/prodClean.nii"), "close", uiMesh_iter, uiMesh_th, uiMesh_bl, uiMesh_smth);
        atrium->ProjectTagsOnExistingSurface(pveins, directory + "/UAC_CT", tagName+".vtk");

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

void AtrialStrainMotionView::MeshPreprocessing(){
    MITK_INFO << "TIMELOG|MeshPreprocessing| Start";
    MITK_INFO << "[MeshPreprocessing] ";
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.
    if (!LoadSurfaceChecks()) return;

    //Show the plugin
    this->GetSite()->GetPage()->ResetPerspective();
    AtrialFibresClipperView::output_directory(); // directory is empty
    // AtrialFibresClipperView::SetDirectoryFile(directory, tagName+".vtk", true);
    this->GetSite()->GetPage()->ShowView("org.mitk.views.atrialfibresclipperview");
    AtrialFibresClipperView::output_directory(); // directory is empty
}


void AtrialStrainMotionView::ClipperPV(){

	MITK_INFO << "TIMELOG|MeshPreprocessing| End";
    MITK_INFO << "TIMELOG|ClipperPV| Start";
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.#
    if (!LoadSurfaceChecks()) return; // Surface was not loaded and user could not find file.

    QString prodPath = Path("UAC_CT/");

    MITK_INFO << "[ClipperPV] clipping PVs.";

	tagName = "Labelled";
    QString path = prodPath + tagName + ".vtk";

	MITK_INFO << path.toStdString();

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

		// mitk::Surface::Pointer la_msh = mitk::IOUtil::Load<mitk::Surface>(la_msh_path.toStdString());
        CemrgCommonUtils::AddToStorage(surface, "clipped", this->GetDataStorage());

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

void AtrialStrainMotionView::ClipperMV(){

   MITK_INFO << "[ClipperMV] Clipping mitral valve";
   atrium->ClipMitralValveAuto(directory + "/UAC_CT", "prodMVI.nii", tagName+".vtk");
}

void AtrialStrainMotionView::UacCalculationVerifyLabels(){
	tagName = "Labelled";

	MITK_INFO << tagName;
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
	MITK_INFO << tagName;

    if (!LoadSurfaceChecks()) return;
    MITK_INFO << ("Loaded surface " + tagName).toStdString();
    std::string title, msg;
    mitk::Surface::Pointer surface = mitk::IOUtil::Load<mitk::Surface>(StdStringPath("UAC_CT/" + tagName + ".vtk"));
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
            mitk::IOUtil::Save(surface, StdStringPath("UAC_CT/" + tagName + ".vtk"));

            title = "Labels fixed";
            msg = "Labels fixed - saved file: \n";
            msg += StdStringPath("UAC_CT/" + tagName + ".vtk");
            QMessageBox::information(NULL, title.c_str(), msg.c_str());

        } else{
            title = "Opening Mesh Preprocessing";
            msg = "Use the Fix labelling button to fix labelling connectivity";
            QMessageBox::information(NULL, title.c_str(), msg.c_str());

            //Show the plugin
            this->GetSite()->GetPage()->ResetPerspective();
            // AtrialFibresClipperView::SetDirectoryFile(directory, tagName+".vtk", true);
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

void AtrialStrainMotionView::MeshImprovement() {
    if (!RequestProjectDirectoryFromUser()) return; // if the path was chosen incorrectly -> returns.#

    // Mesh Improvement: Resampling
    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);
    QString meshname = "Labelled";
    QString outname = "Labelled-refined";
    double hmax = 0.5;
    double hmin = 0.1;
    double havg = 0.3;
    double surfCorr = 0.95;
    QString refinedPath = cmd->DockerRemeshSurface(Path("UAC_CT/"), meshname, outname, hmax, hmin, havg, surfCorr);

    // Mesh Improvement: Cleaning
    cmd->ExecuteTouch(Path("UAC_CT/") + "clean-Labelled-refined.vtk");
    cmd->DockerCleanMeshQuality(Path("UAC_CT/"), "Labelled-refined", "clean-Labelled-refined", 0.2, "vtk", "vtk_polydata");
    cmd->DockerCleanMeshQuality(Path("UAC_CT/"), "clean-Labelled-refined", "clean-Labelled-refined", 0.1, "vtk", "vtk_polydata");

    // Convert Format
    cmd->DockerConvertMeshFormat(Path("UAC_CT/"), "clean-Labelled-refined", "vtk", "clean-Labelled-refined", "carp_txt", 1000);

}

void AtrialStrainMotionView::AutoLandMark() {
    if (!RequestProjectDirectoryFromUser()) return;

    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);

    QString outputFullPath = Path("UAC_CT/") + "prodRefinedLandmarks.txt";
    cmd->ExecuteTouch(outputFullPath);
    cmd->DockerAtrialStrainMotion(Path(), "autoLM");

    QFileInfo finfo(outputFullPath);
    bool fileExists = finfo.exists();

    MITK_INFO << (fileExists ? "Final decision: Done" : "Final decision: Output file not found.");
}

void AtrialStrainMotionView::UAC_Stage1(){
    if (!RequestProjectDirectoryFromUser()) return;
    // docker run --rm --volume=/home/bzhou6/Desktop/V-0005/UAC_CT:/data cemrg/uac:3.0-beta uac --uac-stage 1 --atrium la --layer endo --fibre l --msh clean-Labelled-refined --tags 1 19 11 13 15 17 --landmarks prodRefinedLandmarks.txt --scale 1000

    MITK_INFO << "TIMELOG|UacCalculation_Stage1| Start";
    // Cemrg CMD
    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);
    MITK_INFO << "Do Rough UAC code from Docker";

    QString uacOutput;
    QStringList outputFiles;

    QStringList tags;
    tags << "1" << "19" << "11" << "13" << "15" << "17";
    QStringList landmarks;
    landmarks << "prodRefinedLandmarks.txt";
    uacOutput = cmd->DockerUacMainMode(directory + "/UAC_CT", "1", "la", "endo", "l", "clean-Labelled-refined", tags, landmarks, false, false, 1000);
    MITK_INFO << ("TIMELOG|UacCalculation_Stage1| UAC 1 end " + uacOutput).toStdString();

    outputFiles << "LSbc1.vtx" << "LSbc2.vtx";
    outputFiles << "PAbc1.vtx" << "PAbc2.vtx";

    if (!IsOutputFileCorrect(directory + "/UAC_CT", outputFiles)){
        MITK_INFO << "TIMELOG|UacCalculation_Stage1| End (FAIL)";
        return;
    }

    MITK_INFO << "Create Laplace Solve files for LR and PA SOLVES";
    QString lr_par, pa_par;
    MITK_INFO << "TIMELOG|UacCalculation_Stage1| openCARP start";
    lr_par = CemrgCommonUtils::OpenCarpParamFileGenerator(directory + "/UAC_CT", "carpf_laplace_LS.par", "clean-Labelled-refined", "LSbc1", "LSbc2");
    pa_par = CemrgCommonUtils::OpenCarpParamFileGenerator(directory + "/UAC_CT", "carpf_laplace_PA.par", "clean-Labelled-refined", "PAbc1", "PAbc2");

    MITK_INFO << "Do Laplace Solves using Docker";
    cmd->SetDockerImageOpenCarp();
    QString lrLapSolve, paLapSolve;
    lrLapSolve = cmd->OpenCarpDocker(directory + "/UAC_CT", lr_par, "LR_UAC_N2");
    paLapSolve = cmd->OpenCarpDocker(directory + "/UAC_CT" , pa_par, "PA_UAC_N2");
    MITK_INFO << ("TIMELOG|UacCalculation_Stage1| openCARP end " + lrLapSolve + " - " + paLapSolve).toStdString();

    bool uacOutputSuccess = IsOutputFileCorrect(directory + "/UAC_CT", outputFiles);
    MITK_ERROR(!uacOutputSuccess) << "problem in UAC_Stage1";
    std::string msg = "UAC Calculation - Stage 1 ";
    msg += (uacOutputSuccess) ? "successful" : "failed";
    QMessageBox::information(NULL, "Attention", msg.c_str());
    MITK_INFO << "TIMELOG|UacCalculation_Stage1| End";
}

void AtrialStrainMotionView::UAC_Stage2(){
    if (!RequestProjectDirectoryFromUser()) return;
    // Cemrg CMD
    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);
    MITK_INFO << "Do Rough UAC code from Docker";

    QStringList outputFiles;
    outputFiles << "AnteriorMesh.elem" << "PosteriorMesh.elem";
    outputFiles << "Ant_Strength_Test_PA1.vtx" << "Ant_Strength_Test_LS1.vtx";
    outputFiles << "Post_Strength_Test_PA1.vtx" << "Post_Strength_Test_LS1.vtx";
    MITK_INFO << "TIMELOG|UacCalculation_Stage2| UAC 2.1 - Start";

    QStringList tags;
    tags << "1" << "19" << "11" << "13" << "15" << "17";
    QStringList landmarks;
    landmarks << "prodRefinedLandmarks.txt";
    QString uacOutput = cmd->DockerUacMainMode(directory + "/UAC_CT", "2a", "la", "endo", "l", "clean-Labelled-refined", tags, landmarks, false, false, 1000);
    MITK_INFO << ("TIMELOG|UacCalculation_Stage2| UAC 2.1 - End " + uacOutput).toStdString();

    if (!IsOutputFileCorrect(directory + "/UAC_CT", outputFiles)){
        MITK_INFO << "TIMELOG|UacCalculation_Stage2| End (FAILED)";
        return;
    }

    QString lrp_par, udp_par, lra_par, uda_par;
    QString carpf_lr = "carpf_laplace_single_LR";
    QString carpf_ud = "carpf_laplace_single_UD";

    MITK_INFO << "TIMELOG|UacCalculation_Stage2| openCARP - Start";
    lrp_par = CemrgCommonUtils::OpenCarpParamFileGenerator(directory + "/UAC_CT", carpf_lr+"_P.par", "PosteriorMesh", "", "Post_Strength_Test_LS1");
    udp_par = CemrgCommonUtils::OpenCarpParamFileGenerator(directory + "/UAC_CT", carpf_ud+"_P.par", "PosteriorMesh", "", "Post_Strength_Test_PA1");
    lra_par = CemrgCommonUtils::OpenCarpParamFileGenerator(directory + "/UAC_CT", carpf_lr+"_A.par", "AnteriorMesh", "", "Ant_Strength_Test_LS1");
    uda_par = CemrgCommonUtils::OpenCarpParamFileGenerator(directory + "/UAC_CT", carpf_ud+"_A.par", "AnteriorMesh", "", "Ant_Strength_Test_PA1");

    cmd->SetDockerImageOpenCarp();
    cmd->ExecuteTouch(Path("UAC_CT/") + "Labelled_Coords_2D_Rescaling_v3_C.vtk");

    QString lrpLapSolve, udpLapSolve, lraLapSolve, udaLapSolve;
    lrpLapSolve = cmd->OpenCarpDocker(directory + "/UAC_CT", lrp_par, "LR_Post_UAC");
    udpLapSolve = cmd->OpenCarpDocker(directory + "/UAC_CT", udp_par, "UD_Post_UAC");
    lraLapSolve = cmd->OpenCarpDocker(directory + "/UAC_CT", lra_par, "LR_Ant_UAC");
    udaLapSolve = cmd->OpenCarpDocker(directory + "/UAC_CT", uda_par, "UD_Ant_UAC");
    MITK_INFO << ("TIMELOG|UacCalculation_Stage2| openCARP - End " + lrpLapSolve + "-" + udpLapSolve + "-" + lraLapSolve + "-" + udaLapSolve).toStdString();


    outputFiles.clear();
    outputFiles << "Labelled_Coords_2D_Rescaling_v3_C.vtk";
    outputFiles << "Labelled_Coords_2D_Rescaling_v3_C.elem";
    outputFiles << "Labelled_Coords_2D_Rescaling_v3_C.pts";
    MITK_INFO << "TIMELOG|UacCalculation_Stage2| UAC 2.2 - Start";
    uacOutput = cmd->DockerUacMainMode(directory + "/UAC_CT", "2b", "la", "endo", "l", "clean-Labelled-refined", tags, landmarks, false, false, 1000);
    MITK_INFO << "TIMELOG|UacCalculation_Stage2| UAC 2.2 - End";

    bool uacOutputSuccess = IsOutputFileCorrect(directory + "/UAC_CT", outputFiles);
    MITK_ERROR(!uacOutputSuccess) << "Problem with UAC_Stage2";
    std::string msg = "UAC Calculation - Stage 2 ";
    msg += (uacOutputSuccess) ? "successful" : "failed";
    QMessageBox::information(NULL, "Attention", msg.c_str());

    MITK_INFO << ("TIMELOG|UacCalculation_Stage2| End " + uacOutput).toStdString();
}

void AtrialStrainMotionView::CreateModel() {
    if (!RequestProjectDirectoryFromUser()) return;

    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);

    // Fibre Map
    // TODO: this might cause issue. Labelled_Endo.lon --> outputLabelled_Endo.lon
    cmd->DockerUacFibreMappingMode(directory + "/UAC_CT", "la", "bilayer", "a", "clean-Labelled-refined", true, "output");

    MITK_INFO(QFile::copy(directory + "/UAC_CT/" + "outputLabelled_Endo.lon", directory + "/UAC_CT/" + "outputLabelled_Endo_avg.lon")) << "Copied Labelled_Endo.lon";
    MITK_INFO(QFile::copy(directory + "/UAC_CT/" + "outputLabelled_Epi.lon", directory + "/UAC_CT/" + "outputLabelled_Epi_avg.lon")) << "Copied Labelled_Epi.lon";

    // generateLabelledMsh
    cmd->ExecuteTouch(Path("UAC_CT/") + "clean-Labelled-refined-aligned.vtk");
    cmd->ExecuteTouch(Path("UAC_CT/") + "clean-Labelled-refined-aligned.vtp");
    cmd->DockerAtrialStrainMotion(Path(), "generateLabelledMsh");

    // appendFibres
    cmd->DockerAtrialStrainMotion(Path(), "appendFibres");

    // rotateSegMsh
    cmd->DockerAtrialStrainMotion(Path(), "rotateSegMsh");

    // alignMesh
    cmd->DockerAtrialStrainMotion(Path(), "alignMesh");
}




void AtrialStrainMotionView::CropImage() {
    MITK_INFO << "crop image";
    // TODO: automate finding x_start, y_start, z_start ...
    if (!RequestProjectDirectoryFromUser()) return;

    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);

    QDir nifti_dir = QDir(directory + "/nifti/");
    // QString input_file_name = nifti_dir.entryList(QDir::AllEntries | QDir::NoDotAndDotDot)[8];
    QString input_file_path = directory + "/" + "dcm-10.nii";
    MITK_INFO << input_file_path.toStdString();

    MITK_INFO(QFile::copy(nifti_dir.absolutePath() + "/" + "dcm-10.nii", input_file_path)) << "Copied " << "dcm-10.nii";
    QString pathToSegmentation = cmd->DockerCctaMultilabelSegmentation(directory, input_file_path, false);
    mitk::Image::Pointer segmentation = mitk::IOUtil::Load<mitk::Image>(pathToSegmentation.toStdString());
    std::unique_ptr<CemrgMultilabelSegmentationUtils> cmls(new CemrgMultilabelSegmentationUtils());
    segmentation = cmls->ReplaceLabel(segmentation, 8, 4);
    segmentation = cmls->ReplaceLabel(segmentation, 9, 4);
    segmentation = cmls->ReplaceLabel(segmentation, 10, 4);

    using LabelStatisticsFilterType = itk::LabelStatisticsImageFilter<ImageType, ImageType>;
    ImageType::Pointer itkImage = ImageType::New();
    mitk::CastToItkImage(segmentation, itkImage);
    LabelStatisticsFilterType::Pointer labelStats = LabelStatisticsFilterType::New();
    labelStats->SetInput(itkImage);
    labelStats->SetLabelInput(itkImage); // Use the segmentation image as label
    labelStats->Update();

    LabelStatisticsFilterType::LabelPixelType labelValue = 4;
    LabelStatisticsFilterType::BoundingBoxType bb = labelStats->GetBoundingBox(labelValue);

     std::vector<int> cogIndx{};
     for (int ix = 1; ix <= 6; ix++) {
         int value = bb[ix - 1];
         if (ix % 2 == 0) {
             int max_value = value + value * 0.3;
             std::cout << max_value << std::endl;
             cogIndx.push_back(max_value);
         }
         else {
             int min_value = value - value * 0.3;
             std::cout << min_value << std::endl;
             cogIndx.push_back(min_value);
         }
    }

    for (auto & c1 : cogIndx)
        std::cout << c1 << " ";
    std::cout << std::endl;
    for (auto & c : bb)
        std::cout << c << " ";
    std::cout << std::endl;




    // for (int frame = 0; frame <= 20; frame ++) {
    //     QString cropped_img_name = "dcm-crop-" + QString::number(frame) + ".nii";
    //    cmd->ExecuteTouch(Path("nifti/") + cropped_img_name);
    // }

    // cmd->DockerAtrialStrainMotion(Path(), "cropImages");
}

void AtrialStrainMotionView::Registration() {
    MITK_INFO << "Registration";
    if (!RequestProjectDirectoryFromUser()) return;
    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);

    cmd->ExecuteTouch(Path("tracking/") + "Final.cfg");
    cmd->ExecuteTouch(Path("tracking/") + "imgTimes.lst");

    cmd->DockerAtrialStrainMotion(Path(), "registration");

    cmd->ExecuteTracking(Path(), Path("tracking/imgTimes.lst"), Path("tracking/Final.cfg"), Path("tracking/tsffd.dof"));
}

void AtrialStrainMotionView::DeformMesh() {
    MITK_INFO << "DeformMesh";
    if (!RequestProjectDirectoryFromUser()) return;
    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);
    for (int frame = 0; frame < 10; frame ++) {
        cmd->ExecuteTransformationOnPoints(directory,
                                           Path("UAC_CT/clean-Labelled-refined-aligned.vtp"),
                                           QString("tracking/cLr-aligned-") + QString::number(frame) + ".vtp",
                                           QString("tracking/tsffd.dof"),
                                           frame * 10);
    }


}

void AtrialStrainMotionView::GenerateCellAreaStrains() {
    if (!RequestProjectDirectoryFromUser()) return;

    MITK_INFO << "GenerateCellAreaStrains";

    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);

    for (int frame = 0; frame < 10; frame ++) {
        QString clr_aligned = "cLr-aligned-" + QString::number(frame) + "-areas.txt";
        cmd->ExecuteTouch(Path("tracking/") + clr_aligned);

        QString area_strains = "area-strains-" + QString::number(frame) + ".csv";
        cmd->ExecuteTouch(Path("tracking/") + area_strains);
    }

    cmd->DockerAtrialStrainMotion(Path(), "generateCellAreaStrains");
}

void AtrialStrainMotionView::JacobianThreshold() {
    MITK_INFO << "JacobianThreshold";

    if (!RequestProjectDirectoryFromUser()) return;

    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);

    for (int frame = 1; frame < 10; frame ++) {
        QString clr_aligned = "cLr-fibres-aligned-" + QString::number(frame) + "-scal_jacob-thr0.2.txt";
        cmd->ExecuteTouch(Path("tracking/") + clr_aligned);
    }


    cmd->DockerAtrialStrainMotion(Path(), "jacobianThreshold");
}

void AtrialStrainMotionView::PlotAreaStrain() {
    if (!RequestProjectDirectoryFromUser()) return;

    MITK_INFO << "PlotAreaStrain";
    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);
    cmd->DockerAtrialStrainMotion(Path(), "plotAreaStrain");
}

void AtrialStrainMotionView::CalcFiberStrains() {
    if (!RequestProjectDirectoryFromUser()) return;

    MITK_INFO << "CalcFiberStrains";
    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);
    cmd->DockerAtrialStrainMotion(Path(), "calcFiberStrains");
}

void AtrialStrainMotionView::PlotFibersTrains() {
    if (!RequestProjectDirectoryFromUser()) return;

    MITK_INFO << "PlotFibersTrains";
    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    cmd->SetUseDockerContainers(true);
    cmd->DockerAtrialStrainMotion(Path(), "plotFiberStrains");
}


