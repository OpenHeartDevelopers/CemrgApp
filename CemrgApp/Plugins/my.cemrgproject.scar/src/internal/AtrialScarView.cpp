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
 * Atrial Scar (AS) Plugin for MITK
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
#include <berryFileEditorInput.h>
#include <berryIWorkbenchPage.h>
#include <berryISelectionProvider.h>

// Qmitk
#include <QmitkIOUtil.h>
#include <mitkIOUtil.h>
#include <mitkProgressBar.h>
#include <mitkLookupTableProperty.h>
#include <mitkVtkScalarModeProperty.h>
#include <mitkIDataStorageService.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkDataStorageEditorInput.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include "mitkPluginActivator.h"
#include "AtrialScarView.h"
#include "AtrialScarClipperView.h"

// VTK
#include <vtkPolyData.h>
#include <vtkLookupTable.h>
#include <vtkSphereSource.h>

// ITK
#include <itkResampleImageFilter.h>
#include <itkBinaryCrossStructuringElement.h>
#include <itkGrayscaleErodeImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>

// MyCemrgLib
#include <CemrgAtriaClipper.h>
#include <CemrgCommandLine.h>
#include <CemrgMeasure.h>
#include <numeric>


const std::string AtrialScarView::VIEW_ID = "my.cemrgproject.views.atrialscarview";

void AtrialScarView::CreateQtPartControl(QWidget *parent) {
    
    //Create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, SIGNAL(clicked()), this, SLOT(LoadDICOM()));
    connect(m_Controls.button_2, SIGNAL(clicked()), this, SLOT(ProcessIMGS()));
    connect(m_Controls.button_3, SIGNAL(clicked()), this, SLOT(SegmentIMGS()));
    connect(m_Controls.button_4, SIGNAL(clicked()), this, SLOT(CreateSurf()));
    connect(m_Controls.button_5, SIGNAL(clicked()), this, SLOT(ScarMap()));
    connect(m_Controls.button_6, SIGNAL(clicked()), this, SLOT(Threshold()));
    connect(m_Controls.button_x, SIGNAL(clicked()), this, SLOT(Registration()));
    connect(m_Controls.button_y, SIGNAL(clicked()), this, SLOT(ClipPVeins()));
    connect(m_Controls.button_z, SIGNAL(clicked()), this, SLOT(ClipperMV()));
    connect(m_Controls.button_s, SIGNAL(clicked()), this, SLOT(Sphericity()));
    connect(m_Controls.button_r, SIGNAL(clicked()), this, SLOT(ResetMain()));

    //Sub-buttons signals
    connect(m_Controls.button_2_1, SIGNAL(clicked()), this, SLOT(ConvertNII()));
    connect(m_Controls.button_x_1, SIGNAL(clicked()), this, SLOT(Register()));
    connect(m_Controls.button_x_2, SIGNAL(clicked()), this, SLOT(Transform()));
    connect(m_Controls.button_z_1, SIGNAL(clicked()), this, SLOT(SelectLandmarks()));
    connect(m_Controls.button_z_2, SIGNAL(clicked()), this, SLOT(ClipMitralValve()));
    
    //Set visibility of buttons
    m_Controls.button_2_1->setVisible(false);
    m_Controls.button_x_1->setVisible(false);
    m_Controls.button_x_2->setVisible(false);
    m_Controls.button_z_1->setVisible(false);
    m_Controls.button_z_2->setVisible(false);

    //Default values
    fileName = ".nii";
}

void AtrialScarView::SetFocus() {

    m_Controls.button_1->setFocus();
}

void AtrialScarView::LoadDICOM() {

    //Use MITK DICOM editor
    QString editor_id = "org.mitk.editors.dicomeditor";
    berry::IEditorInput::Pointer input(new berry::FileEditorInput(QString()));
    this->GetSite()->GetPage()->OpenEditor(input, editor_id);
}

void AtrialScarView::ProcessIMGS() {

    //Toggle visibility of buttons
    if (m_Controls.button_2_1->isVisible())
        m_Controls.button_2_1->setVisible(false);
    else
        m_Controls.button_2_1->setVisible(true);
}

void AtrialScarView::ConvertNII() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 2) {
        QMessageBox::warning(NULL, "Attention", "Please load and select both LGE and CEMRA images from the Data Manager to convert!");
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
        QMessageBox::warning(
                    NULL, "Attention",
                    "Cannot find the type of images automatically. Revert to user order and selections in the data manager: LGE at the top, then CEMRA at the bottom!");
        index.resize(nodes.size());
        std::iota(index.begin(), index.end(), 0);
    }//_if

    //Convert to Nifti
    int ctr = 0;
    QString path, type;
    this->BusyCursorOn();
    mitk::ProgressBar::GetInstance()->AddStepsToDo(index.size());

    foreach (int idx, index) {
        mitk::BaseData::Pointer data = nodes.at(idx)->GetData();
        if (data) {
            //Test if this data item is an image
            mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
            if (image) {

                //Resample image to be iso
                typedef itk::Image<short,3> ImageType;
                typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleImageFilterType;
                typedef itk::LinearInterpolateImageFunction<ImageType, double> LinearInterpolatorType;
                ImageType::Pointer itkImage = ImageType::New();
                mitk::CastToItkImage(image, itkImage);
                ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
                LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();
                resampler->SetInterpolator(interpolator);
                resampler->SetInput(itkImage);
                resampler->SetOutputOrigin(itkImage->GetOrigin());
                ImageType::SizeType input_size = itkImage->GetLargestPossibleRegion().GetSize();
                ImageType::SpacingType input_spacing = itkImage->GetSpacing();
                ImageType::SizeType output_size;
                ImageType::SpacingType output_spacing;
                output_size[0] = input_size[0] * (input_spacing[0] / 1.0);
                output_size[1] = input_size[1] * (input_spacing[1] / 1.0);
                output_size[2] = input_size[2] * (input_spacing[2] / 1.0);
                output_spacing [0] = 1.0;
                output_spacing [1] = 1.0;
                output_spacing [2] = 1.0;
                resampler->SetSize(output_size);
                resampler->SetOutputSpacing(output_spacing);
                resampler->SetOutputDirection(itkImage->GetDirection());
                resampler->UpdateLargestPossibleRegion();
                ImageType::Pointer resampledImage = resampler->GetOutput();
                image = mitk::ImportItkImage(resampledImage)->Clone();

                type = (ctr==0) ? "LGE":"MRA";
                path = directory + mitk::IOUtil::GetDirectorySeparator() + "dcm-" + type + "-" + seriesDscrps.at(idx).c_str() + ".nii";
                mitk::IOUtil::Save(image, path.toStdString());
                this->GetDataStorage()->Remove(nodes.at(idx));

                std::string key = "dicom.series.SeriesDescription";
                mitk::DataStorage::SetOfObjects::Pointer set = mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
                set->Begin().Value()->GetData()->GetPropertyList()->SetStringProperty(key.c_str(), seriesDscrps.at(idx).c_str());
                ctr++;
            } else
                return;
        } else
            return;
        mitk::ProgressBar::GetInstance()->Progress();
    }//for

    nodes.clear();
    this->BusyCursorOff();

    //Load all items
    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
}

void AtrialScarView::SegmentIMGS() {

    int reply = QMessageBox::question(
                NULL, "Question", "Do you have a segmentation to load?", QMessageBox::Yes, QMessageBox::No);

    if (reply == QMessageBox::Yes) {

        QString path = QFileDialog::getOpenFileName(
                    NULL, "Open Segmentation file", mitk::IOUtil::GetProgramPath().c_str(), QmitkIOUtil::GetFileOpenFilterString());
        if (path.isEmpty()) return;
        mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
        mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());

        //Restore image name
        char sep = mitk::IOUtil::GetDirectorySeparator();
        fileName = path.mid(path.lastIndexOf(sep) + 1);

    } else {
        //Show the plugin
        this->GetSite()->GetPage()->ShowView("org.mitk.views.segmentation");
    }//_if
}

/**
 * @brief Relies on the interpolation mesh being added to DataManger
 * @param node Contains segmentation image to be saved
 */
void AtrialScarView::NodeAdded(const mitk::DataNode* node) {

    QString name = QString::fromStdString(node->GetName());
    if (name.endsWith("interpolation")) {

        mitk::DataStorage::SetOfObjects::ConstPointer nodes = this->GetDataStorage()->GetSources(node);
        mitk::DataStorage::SetOfObjects::ConstIterator itSeg = nodes->Begin();
        mitk::DataNode::Pointer segNode = itSeg->Value();

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
        mitk::BaseData::Pointer data = segNode->GetData();
        if (data) {
            //Test if this data item is an image
            mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
            if (image) {

                bool ok;
                QString tmpFileName = fileName;
                fileName = QInputDialog::getText(
                            NULL, tr("Save Segmentation As"), tr("File Name:"), QLineEdit::Normal, fileName, &ok);
                if (ok && !fileName.isEmpty() && fileName.endsWith(".nii")) {
                    segNode->SetName(fileName.left(fileName.length()-4).toStdString());
                    path = directory + mitk::IOUtil::GetDirectorySeparator() + fileName;
                    mitk::IOUtil::Save(image, path.toStdString());
                } else {
                    fileName = tmpFileName;
                    QMessageBox::warning(NULL, "Attention", "Please type a file name with the right extension (i.e. .nii)!");
                    return;
                }//_fileName

            } else {
                QMessageBox::warning(NULL, "Attention", "Segmentation image could not be saved automatically from the Data Manager!");
                return;
            }//_image
        } else
            return;
        this->GetDataStorage()->Remove(node);

    } else if (name.endsWith("PVeinsCroppedImage")) {

        //Listen to image crop change
        fileName = "PVeinsCroppedImage.nii";

    }//_if
}

void AtrialScarView::Registration() {

    //Toggle visibility of buttons
    if (m_Controls.button_x_1->isVisible()) {
        m_Controls.button_x_1->setVisible(false);
        m_Controls.button_x_2->setVisible(false);
    } else {
        m_Controls.button_x_1->setVisible(true);
        m_Controls.button_x_2->setVisible(true);
    }
}

void AtrialScarView::Register() {

    QMessageBox::StandardButton reply;
    reply = QMessageBox::question(
                NULL, "Image Registration",
                "Have you completed steps 1 to 3 before using this feature?", QMessageBox::Yes|QMessageBox::No);
    if (reply == QMessageBox::No)
        return;

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 2) {
        QMessageBox::warning(NULL, "Attention", "Please select both LGE and CEMRA images from the Data Manager to register!");
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

    //Sort the two images
    QString lge, mra;
    mitk::DataNode::Pointer _1st = nodes.at(0);
    mitk::DataNode::Pointer _2nd = nodes.at(1);
    mitk::BaseData::Pointer dat1 = _1st->GetData();
    mitk::BaseData::Pointer dat2 = _2nd->GetData();

    if (dat1 && dat2) {
        //Test if this data items are image
        mitk::Image::Pointer image1 = dynamic_cast<mitk::Image*>(dat1.GetPointer());
        mitk::Image::Pointer image2 = dynamic_cast<mitk::Image*>(dat2.GetPointer());
        if (image1 && image2) {

            if (_1st->GetName().find("LGE") != _1st->GetName().npos && _2nd->GetName().find("MRA") != _2nd->GetName().npos) {

                lge = QString::fromStdString(_1st->GetName());
                mra = QString::fromStdString(_2nd->GetName());

            } else if (_1st->GetName().find("MRA") != _1st->GetName().npos && _2nd->GetName().find("LGE") != _2nd->GetName().npos) {

                lge = QString::fromStdString(_2nd->GetName());
                mra = QString::fromStdString(_1st->GetName());

            } else {

                QMessageBox::warning(NULL, "Attention", "Please select both LGE and CEMRA images from the Data Manager!");
                return;

            }//_if

            //Commandline call
            this->BusyCursorOn();
            mitk::ProgressBar::GetInstance()->AddStepsToDo(1);
            std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
            cmd->ExecuteRegistration(directory, lge, mra);
            QMessageBox::information(NULL, "Attention", "Command Line Operations Finished!");
            this->BusyCursorOff();

        } else {
            QMessageBox::warning(NULL, "Attention", "Please select both LGE and CEMRA images from the Data Manager!");
            return;
        }//_image
    } else
        return;
}

void AtrialScarView::Transform() {

    QMessageBox::StandardButton reply;
    reply = QMessageBox::question(
                NULL, "Image Transformation",
                "Have you completed image registration before using this step?", QMessageBox::Yes|QMessageBox::No);
    if (reply == QMessageBox::No)
        return;

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1) {
        QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation to transform!");
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
    QString path, pathTemp;
    mitk::DataNode::Pointer segNode = nodes.at(0);
    mitk::BaseData::Pointer data = segNode->GetData();
    if (data) {
        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image) {

            //Check seg node name
            if (segNode->GetName().compare(fileName.left(fileName.length()-4).toStdString()) != 0) {
                QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation!");
                return;
            }//_if

            bool ok;
            QString regFileName;
            regFileName = fileName.left(fileName.length()-4) + "-reg.nii";
            regFileName = QInputDialog::getText(NULL, tr("Save Registration As"), tr("File Name:"), QLineEdit::Normal, regFileName, &ok);
            if (ok && !regFileName.isEmpty() && regFileName.endsWith(".nii")) {

                pathTemp = directory + mitk::IOUtil::GetDirectorySeparator() + "temp.nii";
                mitk::IOUtil::Save(image, pathTemp.toStdString());

                //Commandline call
                this->BusyCursorOn();
                mitk::ProgressBar::GetInstance()->AddStepsToDo(1);
                std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
                cmd->ExecuteTransformation(directory, pathTemp.right(8), regFileName);
                QMessageBox::information(NULL, "Attention", "Command Line Operations Finished!");
                this->BusyCursorOff();

                //Load the new segementation
                path = directory + mitk::IOUtil::GetDirectorySeparator() + regFileName;
                mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
                remove(pathTemp.toStdString().c_str());

                //Clear data manager
                mitk::DataStorage::SetOfObjects::ConstPointer sob = this->GetDataStorage()->GetAll();
                for (mitk::DataStorage::SetOfObjects::ConstIterator nodeIt = sob->Begin(); nodeIt != sob->End(); ++nodeIt)
                    if (nodeIt->Value()->GetName().find("MRA") != nodeIt->Value()->GetName().npos)
                        this->GetDataStorage()->Remove(nodeIt->Value());
                this->GetDataStorage()->Remove(segNode);
                fileName = regFileName;

            } else {
                QMessageBox::warning(NULL, "Attention", "Please type a file name with the right extension (i.e. .nii)!");
                return;
            }//_fileName

        } else {
            QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation!");
            return;
        }//_image
    } else
        return;
}

void AtrialScarView::ClipPVeins() {

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
    AtrialScarClipperView::SetDirectoryFile(directory, fileName);
    this->GetSite()->GetPage()->ShowView("my.cemrgproject.views.atrialscarclipperview");
}

void AtrialScarView::CreateSurf() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1) {
        QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation to create a surface!");
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
    QString path, pathTemp;
    mitk::DataNode::Pointer segNode = nodes.at(0);
    mitk::BaseData::Pointer data = segNode->GetData();
    if (data) {
        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image) {

            //Check seg node name
            if (segNode->GetName().compare(fileName.left(fileName.length()-4).toStdString()) != 0) {
                QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation!");
                return;
            }//_if

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
                if(!ok3) blur = 0;
                if(!ok4) smth = 10;
                //_if

                this->BusyCursorOn();
                pathTemp = directory + mitk::IOUtil::GetDirectorySeparator() + "temp.nii";
                mitk::IOUtil::Save(image, pathTemp.toStdString());
                mitk::ProgressBar::GetInstance()->AddStepsToDo(4);
                std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
                path = cmd->ExecuteSurf(directory, pathTemp, iter, th, blur, smth);
                QMessageBox::information(NULL, "Attention", "Command Line Operations Finished!");
                this->BusyCursorOff();

                //Add the mesh to storage
                std::string meshName = segNode->GetName() + "-Mesh";
                AddToStorage(meshName, ReadVTKMesh(path.toStdString()));
                mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
                inputs->deleteLater();
                remove(pathTemp.toStdString().c_str());

            } else if(dialogCode == QDialog::Rejected) {
                inputs->close();
                inputs->deleteLater();
            }//_if

        } else {
            QMessageBox::warning(NULL, "Attention", "Please select the loaded or created segmentation!");
            return;
        }//_image
    } else
        return;
}

void AtrialScarView::ClipperMV() {

    //Toggle visibility of buttons
    if (m_Controls.button_z_1->isVisible()) {
        m_Controls.button_z_1->setVisible(false);
        m_Controls.button_z_2->setVisible(false);
        return;
    } else {
        m_Controls.button_z_1->setVisible(true);
        m_Controls.button_z_2->setVisible(true);
    }//_if
}

void AtrialScarView::SelectLandmarks() {

    if (m_Controls.button_z_1->text() == QString::fromStdString("Select Landmarks")) {

        //Show the plugin
        QMessageBox::information(NULL, "Attention", "Please select 3 points around the mitral valve!");
        this->GetSite()->GetPage()->ShowView("org.mitk.views.pointsetinteraction");
        m_Controls.button_z_1->setText("Display Clipper");

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

        //Reset the button
        m_Controls.button_z_1->setText("Select Landmarks");

        //Retrieve mean and distance of 3 points
        double x_c = 0;
        double y_c = 0;
        double z_c = 0;
        for(int i=0; i<pointSet->GetSize(); i++) {
            x_c = x_c + pointSet->GetPoint(i).GetElement(0);
            y_c = y_c + pointSet->GetPoint(i).GetElement(1);
            z_c = z_c + pointSet->GetPoint(i).GetElement(2);
        }//_for
        x_c /= pointSet->GetSize();
        y_c /= pointSet->GetSize();
        z_c /= pointSet->GetSize();
        double distance[pointSet->GetSize()];
        for(int i=0; i<pointSet->GetSize(); i++) {
            double x_d = pointSet->GetPoint(i).GetElement(0) - x_c;
            double y_d = pointSet->GetPoint(i).GetElement(1) - y_c;
            double z_d = pointSet->GetPoint(i).GetElement(2) - z_c;
            distance[i] = sqrt(pow(x_d,2) + pow(y_d,2) + pow(z_d,2));
        }//_for
        double radius = *std::max_element(distance, distance + pointSet->GetSize());

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
        AddToStorage("MVClipper", mvClipper);
        sob = this->GetDataStorage()->GetAll();
        for (mitk::DataStorage::SetOfObjects::ConstIterator nodeIt = sob->Begin(); nodeIt != sob->End(); ++nodeIt) {
            if (nodeIt->Value()->GetName().find("MVClipper") != nodeIt->Value()->GetName().npos) {
                nodeIt->Value()->SetProperty("opacity", mitk::FloatProperty::New(0.4));
                nodeIt->Value()->SetProperty("color", mitk::ColorProperty::New(1.0, 0.0, 0.0));
            }//_if
        }//_for

    }//_if
}

void AtrialScarView::ClipMitralValve() {

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

    //Read in and copy
    QString path = directory + mitk::IOUtil::GetDirectorySeparator() + "segmentation.vtk";
    mitk::Surface::Pointer surface = ReadVTKMesh(path.toStdString());
    if (surface->GetVtkPolyData() == NULL) {
        QMessageBox::critical(NULL, "Attention", "No mesh was found in the project directory!");
        return;
    }//_if
    QString orgP = path.left(path.length()-4) + "-Original.vtk";
    mitk::IOUtil::Save(mitk::IOUtil::LoadSurface(path.toStdString()), orgP.toStdString());

    /*
     * Producibility Test
     **/
    QString prodPath = directory + mitk::IOUtil::GetDirectorySeparator();
    mitk::IOUtil::Save(pointSet, (prodPath + "prodMVCLandmarks.mps").toStdString());
    /*
     * End Test
     **/

    this->BusyCursorOn();
    std::unique_ptr<CemrgScar3D> scarObj = std::unique_ptr<CemrgScar3D>(new CemrgScar3D());
    surface = scarObj->ClipMesh3D(surface, pointSet);
    this->BusyCursorOff();

    //Check to remove the previous mesh node
    mitk::DataStorage::SetOfObjects::ConstPointer sob = this->GetDataStorage()->GetAll();
    for (mitk::DataStorage::SetOfObjects::ConstIterator nodeIt = sob->Begin(); nodeIt != sob->End(); ++nodeIt) {
        if (nodeIt->Value()->GetName().find("-Mesh") != nodeIt->Value()->GetName().npos)
            this->GetDataStorage()->Remove(nodeIt->Value());
        if (nodeIt->Value()->GetName().find("MVClipper") != nodeIt->Value()->GetName().npos)
            this->GetDataStorage()->Remove(nodeIt->Value());
    }//_for
    AddToStorage("MVClipped-Mesh", surface);

    //Reverse coordination of surface for writing MIRTK style
    mitk::Surface::Pointer surfCloned = surface->Clone();
    vtkSmartPointer<vtkPolyData> pd = surfCloned->GetVtkPolyData();
    for (int i=0; i<pd->GetNumberOfPoints(); i++) {
        double* point = pd->GetPoint(i);
        point[0] = -point[0];
        point[1] = -point[1];
        pd->GetPoints()->SetPoint(i, point);
    }//_for
    mitk::IOUtil::Save(surfCloned, path.toStdString());
}

void AtrialScarView::ScarMap() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1) {
        QMessageBox::warning(NULL, "Attention", "Please select the LGE image from the Data Manager to calculate the scar map!");
        return;
    }
    std::size_t found = nodes.at(0)->GetName().find("LGE");
    bool check = found != std::string::npos ? false : true;
    if (check) {
        QMessageBox::warning(NULL, "Attention", "Please select the LGE image from the Data Manager to calculate the scar map!");
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

    //Check for mesh in the project directory
    try {
        QString path = directory + mitk::IOUtil::GetDirectorySeparator() + "segmentation.vtk";
        mitk::IOUtil::LoadSurface(path.toStdString());
    } catch (...) {
        QMessageBox::critical(NULL, "Attention", "No mesh was found in the project directory!");
        return;
    }//_try

    //Find the selected node
    mitk::DataNode::Pointer imgNode = nodes.at(0);
    mitk::BaseData::Pointer data = imgNode->GetData();
    if (data) {

        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image.IsNotNull()) {

            scar = std::unique_ptr<CemrgScar3D>(new CemrgScar3D());
            if (scar) {

                //Ask for user input to set the parameters
                QDialog* inputs = new QDialog(0,0);
                m_UIScar.setupUi(inputs);
                connect(m_UIScar.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
                connect(m_UIScar.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
                int dialogCode = inputs->exec();

                //Act on dialog return code
                if (dialogCode == QDialog::Accepted) {

                    bool ok1, ok2;
                    int minStep = m_UIScar.lineEdit_1->text().toInt(&ok1);
                    int maxStep = m_UIScar.lineEdit_2->text().toInt(&ok2);
                    int methodType = m_UIScar.radioButton_1->isChecked() ? 2 : 1;
                    QString meType = m_UIScar.radioButton_1->isChecked() ? "Max" : "Mean";

                    //Set default values
                    if (!ok1 || !ok2)
                        QMessageBox::warning(NULL, "Attention", "Reverting to default parameters!");
                    if(!ok1) minStep =-1;
                    if(!ok2) maxStep = 3;
                    //_if

                    /*
                     * Producibility Test
                     **/
                    QString prodPath = directory + mitk::IOUtil::GetDirectorySeparator();
                    ofstream prodFile1;
                    prodFile1.open((prodPath + "prodScarMapInputs.txt").toStdString());
                    prodFile1 << minStep << "\n";
                    prodFile1 << maxStep << "\n";
                    prodFile1 << meType.toStdString() << "\n";
                    prodFile1.close();
                    /*
                     * End Test
                     **/

                    this->BusyCursorOn();
                    mitk::ProgressBar::GetInstance()->AddStepsToDo(1);

                    scar->SetMinStep(minStep);
                    scar->SetMaxStep(maxStep);
                    scar->SetMethodType(methodType);
                    mitk::Image::Pointer scarSegImg;
                    try {
                        QString path = directory + mitk::IOUtil::GetDirectorySeparator() + fileName;
                        scarSegImg = mitk::IOUtil::LoadImage(path.toStdString());
                    } catch (...) {
                        QMessageBox::critical(NULL, "Attention", "The loaded or created segmentation was not found!");
                        mitk::ProgressBar::GetInstance()->Progress();
                        this->BusyCursorOff();
                        return;
                    }//_try

                    //Resample to fit LGE
                    typedef itk::Image<short, 3> ImageType;
                    itk::ResampleImageFilter<ImageType, ImageType>::Pointer resampleFilter;
                    ImageType::Pointer scarSegITK = ImageType::New();
                    CastToItkImage(scarSegImg, scarSegITK);
                    ImageType::Pointer lgeITK = ImageType::New();
                    CastToItkImage(image, lgeITK);
                    resampleFilter = itk::ResampleImageFilter<ImageType, ImageType >::New();
                    resampleFilter->SetInput(scarSegITK);
                    resampleFilter->SetReferenceImage(lgeITK);
                    resampleFilter->SetUseReferenceImage(true);
                    resampleFilter->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<ImageType>::New());
                    resampleFilter->SetDefaultPixelValue(0);
                    resampleFilter->UpdateLargestPossibleRegion();
                    scarSegImg = mitk::ImportItkImage(resampleFilter->GetOutput());

                    //Update datamanger
                    std::string nodeSegImgName = fileName.left(fileName.length()-4).toStdString();
                    mitk::DataStorage::SetOfObjects::ConstPointer sob = this->GetDataStorage()->GetAll();
                    for (mitk::DataStorage::SetOfObjects::ConstIterator nodeIt = sob->Begin(); nodeIt != sob->End(); ++nodeIt)
                        if (nodeIt->Value()->GetName().find(nodeSegImgName) != nodeIt->Value()->GetName().npos)
                            this->GetDataStorage()->Remove(nodeIt->Value());
                    QString savePath = directory + mitk::IOUtil::GetDirectorySeparator() + fileName;
                    mitk::IOUtil::Save(scarSegImg, savePath.toStdString());

                    //Projection
                    scar->SetScarSegImage(scarSegImg);
                    mitk::Surface::Pointer shell = scar->Scar3D(directory.toStdString(), image);
                    mitk::DataNode::Pointer node = AddToStorage((meType + "Scar3D").toStdString(), shell);
                    mitk::ProgressBar::GetInstance()->Progress();

                    //Check to remove the previous mesh node
                    sob = this->GetDataStorage()->GetAll();
                    for (mitk::DataStorage::SetOfObjects::ConstIterator nodeIt = sob->Begin(); nodeIt != sob->End(); ++nodeIt)
                        if (nodeIt->Value()->GetName().find("-Mesh") != nodeIt->Value()->GetName().npos)
                            this->GetDataStorage()->Remove(nodeIt->Value());

                    mitk::LookupTable::Pointer scarLut = mitk::LookupTable::New();
                    vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
                    lut->SetTableRange(0, scar->GetMaxScalar());
                    lut->SetHueRange(0.33, 0.0);
                    lut->Build();
                    scarLut->SetVtkLookupTable(lut);
                    node->SetProperty("LookupTable", mitk::LookupTableProperty::New(scarLut));
                    node->SetProperty("scalar visibility", mitk::BoolProperty::New(true));
                    node->SetProperty("scalar mode", mitk::VtkScalarModeProperty::New(2));
                    node->SetFloatProperty("ScalarsRangeMinimum", 0);
                    node->SetFloatProperty("ScalarsRangeMaximum", scar->GetMaxScalar());

                    //Save the vtk mesh
                    QString name(imgNode->GetName().c_str());
                    name = name.right(name.length() - name.lastIndexOf("-") - 1);
                    QString path = directory + mitk::IOUtil::GetDirectorySeparator() + name + "-" + meType + "Scar.vtk";
                    mitk::IOUtil::Save(shell, path.toStdString());

                    this->BusyCursorOff();
                    inputs->deleteLater();

                } else if(dialogCode == QDialog::Rejected) {
                    inputs->close();
                    inputs->deleteLater();
                }//_if
            }//_scar
        } else
            return;
    } else
        return;
}

void AtrialScarView::Threshold() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1) {
        QMessageBox::warning(NULL, "Attention", "Please select the LGE image from the Data Manager to quantify the scar!");
        return;
    }
    std::size_t found = nodes.at(0)->GetName().find("LGE");
    bool check = found != std::string::npos ? false : true;
    if (check) {
        QMessageBox::warning(NULL, "Attention", "Please select the LGE image from the Data Manager to quantify the scar!");
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
    double mean = 0.0, stdv = 0.0;
    mitk::DataNode::Pointer imgNode = nodes.at(0);
    mitk::BaseData::Pointer data = imgNode->GetData();
    if (data) {

        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image.IsNotNull()) {

            //Convert images to right type
            mitk::Image::Pointer roi;
            mitk::Image::Pointer lgeImage = mitk::Image::New();
            itk::Image<float,3>::Pointer itkImage = itk::Image<float,3>::New();
            mitk::CastToItkImage(image, itkImage);
            mitk::CastToMitkImage(itkImage, lgeImage);
            try {
                QString path = directory + mitk::IOUtil::GetDirectorySeparator() + fileName;
                roi = mitk::IOUtil::LoadImage(path.toStdString());
            } catch (...) {
                QMessageBox::critical(NULL, "Attention", "The loaded or created segmentation was not found!");
                return;
            }//_try
            mitk::Image::Pointer roiImage = mitk::Image::New();
            itk::Image<float,3>::Pointer roiItkImage = itk::Image<float,3>::New();
            mitk::CastToItkImage(roi, roiItkImage);
            if (scar) {

                //Erosion of bloodpool
                typedef itk::Image<float, 3> ImageType;
                typedef itk::BinaryCrossStructuringElement<ImageType::PixelType, 3> CrossType;
                typedef itk::GrayscaleErodeImageFilter<ImageType, ImageType, CrossType> ErosionFilterType;
                bool ok;
                int vxls = QInputDialog::getInt(NULL, tr("Bloodpool Erosion"), tr("Voxels:"), 3, 1, 20, 1, &ok);
                if (!ok) {
                    QMessageBox::warning(NULL, "Attention", "Enter a correct value!");
                    return;
                }//_if
                CrossType binaryCross;
                binaryCross.SetRadius(vxls);
                binaryCross.CreateStructuringElement();
                ErosionFilterType::Pointer erosionFilter = ErosionFilterType::New();
                erosionFilter->SetInput(roiItkImage);
                erosionFilter->SetKernel(binaryCross);
                erosionFilter->UpdateLargestPossibleRegion();
                roiImage = mitk::ImportItkImage(erosionFilter->GetOutput())->Clone();
                AddToStorage("Eroded ROI", roiImage);

                //Calculate mean, std of ROI
                bool success = scar->CalculateMeanStd(lgeImage, roiImage, mean, stdv);
                if (!success)
                    return;

            } else {
                QMessageBox::warning(NULL, "Attention", "The scar map from the previous step has not been generated!");
                return;
            }//_if_scar
        }//_image
    }//_data

    //Ask for user input to set the parameters
    QDialog* inputs = new QDialog(0,0);
    m_UISQuant.setupUi(inputs);
    connect(m_UISQuant.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_UISQuant.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    int dialogCode = inputs->exec();

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {

        bool ok;
        double value = m_UISQuant.lineEdit->text().toDouble(&ok);
        int methodType = m_UISQuant.radioButton_1->isChecked() ? 1 : 2;

        //Set default values
        if (!ok) {
            QMessageBox::warning(NULL, "Attention", "Please enter a valid value!");
            return;
        }//_ok

        this->BusyCursorOn();
        mitk::ProgressBar::GetInstance()->AddStepsToDo(1);

        double percentage = 0;
        double thresh = (methodType == 1) ? mean*value : mean+value*stdv;

        /*
         * Producibility Test
         **/
        QString prodPath = directory + mitk::IOUtil::GetDirectorySeparator();
        ofstream prodFile1;
        prodFile1.open((prodPath + "prodThresholds.txt").toStdString());
        prodFile1 << value << "\n";
        prodFile1 << methodType << "\n";
        prodFile1 << mean << "\n";
        prodFile1 << stdv << "\n";
        prodFile1 << thresh << "\n";
        prodFile1.close();
        /*
         * End Test
         **/

        if (scar) percentage = scar->Thresholding(thresh);
        std::ostringstream os;
        os << std::fixed << std::setprecision(2) << percentage;
        QString message = "The percentage scar is " + QString::fromStdString(os.str()) + "% of total segmented volume.";
        QMessageBox::information(NULL, "Scar Quantification", message);

        mitk::ProgressBar::GetInstance()->Progress();
        this->BusyCursorOff();
        inputs->deleteLater();

    } else if(dialogCode == QDialog::Rejected) {
        inputs->close();
        inputs->deleteLater();
    }//_if
}

void AtrialScarView::Sphericity() {

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

    //Read in the mesh
    QString path = directory + mitk::IOUtil::GetDirectorySeparator() + "segmentation.vtk";
    mitk::Surface::Pointer surface = ReadVTKMesh(path.toStdString());

    double result = 0;
    std::unique_ptr<CemrgMeasure> Sphericity = std::unique_ptr<CemrgMeasure>(new CemrgMeasure());
    if (surface->GetVtkPolyData() != NULL)
        result = Sphericity->GetSphericity(surface->GetVtkPolyData());
    else {
        QMessageBox::warning(NULL, "Attention", "No mesh was found in the project directory!");
        return;
    }//_if

    QString message = "The sphericity is " + QString::number(result) + "%.";
    QMessageBox::information(NULL, "Sphericity", message);
}

void AtrialScarView::ResetMain() {

    Reset(true);
}

void AtrialScarView::Reset(bool allItems) {

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

        //Remove everything or keep images
        mitk::DataStorage::SetOfObjects::ConstPointer nodesToRemove = dataStorage->GetAll();
        if (allItems) {
            dataStorage->Remove(nodesToRemove);
        } else {
            mitk::DataStorage::SetOfObjects::ConstIterator iter = nodesToRemove->Begin();
            mitk::DataStorage::SetOfObjects::ConstIterator iterEnd = nodesToRemove->End();
            for (; iter != iterEnd; ++iter) {
                QString name = QString::fromStdString(iter->Value()->GetName());
                if (name.startsWith("dcm-"))
                    continue;
                dataStorage->Remove(iter->Value());
            }//_for
        }//_if

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
    fileName = ".nii";
    directory.clear();
    m_Controls.button_z_1->setText("Select Landmarks");
    this->GetSite()->GetPage()->ResetPerspective();
}

mitk::DataNode::Pointer AtrialScarView::AddToStorage(std::string nodeName, mitk::BaseData* data) {

    mitk::DataNode::Pointer node = mitk::DataNode::New();
    if (!data) return node;
    node->SetData(data);
    node->SetName(nodeName);
    this->GetDataStorage()->Add(node);
    return node;
}

mitk::Surface::Pointer AtrialScarView::ReadVTKMesh(std::string meshPath) {

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
