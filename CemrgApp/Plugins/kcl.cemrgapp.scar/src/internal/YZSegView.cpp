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
 * Scar Quantification
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * YZ
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>
#include <berryFileEditorInput.h>
#include <berryIWorkbenchPage.h>

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
#include <mitkImage.h>
#include "kcl_cemrgapp_scar_Activator.h"
#include "YZSegView.h"

// VTK
#include <vtkPolyData.h>
#include <vtkLookupTable.h>

// ITK
#include <itkResampleImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImage.h>

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>

const std::string YZSegView::VIEW_ID = "org.mitk.views.scaryzseg";

void YZSegView::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, SIGNAL(clicked()), this, SLOT(LoadDICOM()));
    connect(m_Controls.button_2, SIGNAL(clicked()), this, SLOT(ProcessIMGS()));
    connect(m_Controls.button_3, SIGNAL(clicked()), this, SLOT(SegmentIMGS()));
    connect(m_Controls.button_4, SIGNAL(clicked()), this, SLOT(ScarSeg()));
    connect(m_Controls.button_2_1, SIGNAL(clicked()), this, SLOT(ConvertNII()));
    connect(m_Controls.button_2_2, SIGNAL(clicked()), this, SLOT(ResampIMGS()));
    connect(m_Controls.button_3_1, SIGNAL(clicked()), this, SLOT(SaveSEG()));
    connect(m_Controls.button_4_1, SIGNAL(clicked()), this, SLOT(ScarSeg_FWHM()));
    connect(m_Controls.button_4_2, SIGNAL(clicked()), this, SLOT(ScarSeg_SD()));
    connect(m_Controls.button_4_2_1, SIGNAL(clicked()), this, SLOT(ScarSeg_4SD()));
    connect(m_Controls.button_4_2_2, SIGNAL(clicked()), this, SLOT(ScarSeg_6SD()));
    connect(m_Controls.button_4_2_3, SIGNAL(clicked()), this, SLOT(ScarSeg_CustomisedSD()));
    connect(m_Controls.button_4_3, SIGNAL(clicked()), this, SLOT(ScarSeg_save()));

    //Set visibility of buttons
    m_Controls.label_1->setVisible(false);
    m_Controls.label_2->setVisible(false);
    m_Controls.label_3->setVisible(false);
    m_Controls.label_4->setVisible(false);
    m_Controls.label_5->setVisible(false);
    m_Controls.label_6->setVisible(false);
    m_Controls.button_2_1->setVisible(false);
    m_Controls.button_2_2->setVisible(false);
    m_Controls.button_3_1->setVisible(false);
    m_Controls.button_4_1->setVisible(false);
    m_Controls.button_4_2->setVisible(false);
    m_Controls.button_4_2_1->setVisible(false);
    m_Controls.button_4_2_2->setVisible(false);
    m_Controls.button_4_2_3->setVisible(false);
    m_Controls.button_4_3->setVisible(false);

    //Default values
    fileName = ".nii";
}

void YZSegView::SetFocus() {
    m_Controls.button_1->setFocus();
}

void YZSegView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*source*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

void YZSegView::LoadDICOM() {
    //Use MITK DICOM editor
    QString editor_id = "org.mitk.editors.dicomeditor";
    berry::IEditorInput::Pointer input(new berry::FileEditorInput(QString()));
    this->GetSite()->GetPage()->OpenEditor(input, editor_id);
}

void YZSegView::ProcessIMGS() {
    //Toggle visibility of buttons
    if (m_Controls.button_2_1->isVisible()) {
        m_Controls.button_2_1->setVisible(false);
        m_Controls.button_2_2->setVisible(false);
    } else {
        m_Controls.button_2_1->setVisible(true);
        m_Controls.button_2_2->setVisible(true);
    }
}

void YZSegView::ConvertNII() {
    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1) {
        QMessageBox::warning(
            NULL, "Attention",
            "Please load and select an LGE image from the Data Manager before starting this step!");
        return;
    }//_if

    //Ask the user for a dir to store data
    if (directory.isEmpty()) {
        directory = QFileDialog::getExistingDirectory(
            NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
            QFileDialog::ShowDirsOnly | QFileDialog::DontUseNativeDialog);
        if (directory.isEmpty() || directory.simplified().contains(" ")) {
            QMessageBox::warning(NULL, "Attention", "Please select a project directory with no spaces in the path!");
            directory = QString();
            return;
        }//_if
    }

    QString path;
    mitk::DataNode::Pointer imNode = nodes.at(0);
    mitk::BaseData::Pointer data = imNode->GetData();
    if (data) {

        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image) {

            bool ok;
            QString tmpFileName = fileName;
            fileName = QInputDialog::getText(NULL, tr("Save As"), tr("File Name:"), QLineEdit::Normal, fileName, &ok);
            if (ok && !fileName.isEmpty() && fileName.endsWith(".nii")) {

                path = directory + "/" + fileName;
                mitk::IOUtil::Save(image, path.toStdString());
                mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
                this->GetDataStorage()->Remove(imNode);
                mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());

            } else {

                fileName = tmpFileName;
                QMessageBox::warning(NULL, "Attention", "Please type a file name with the right extension (i.e. .nii)!");
                return;
            }//_fileName

        } else {
            QMessageBox::warning(NULL, "Attention", "Please select a DICOM image from the Data Manager!");
            return;
        }//_image
    } else
        return;
}

void YZSegView::ResampIMGS() {

    //Show the plugin
    QMessageBox::critical(
        NULL, "Attention",
        "Make sure you close this view before closing the application! Failure to do so will cause serious problems!");
    this->GetSite()->GetPage()->ShowView("org.mitk.views.basicimageprocessing");
}

void YZSegView::SegmentIMGS() {

    //Toggle visibility of buttons
    if (m_Controls.button_3_1->isVisible()) {
        m_Controls.button_3_1->setVisible(false);
        return;
    } else {
        m_Controls.button_3_1->setVisible(true);
    }

    int reply = QMessageBox::question(
        NULL, "Question", "Do you have a segmentation to load?", QMessageBox::Yes, QMessageBox::No);

    if (reply == QMessageBox::Yes) {

        QString path = QFileDialog::getOpenFileName(
            NULL, "Open Segmentation file", mitk::IOUtil::GetProgramPath().c_str(), QmitkIOUtil::GetFileOpenFilterString());
        if (path.isEmpty()) return;
        mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
        mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());

    } else {
        //Show the plugin
        this->GetSite()->GetPage()->ShowView("org.mitk.views.segmentation");
    }//_if
}

void YZSegView::SaveSEG() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1) {
        QMessageBox::warning(NULL, "Attention", "Please select a segmentation from the Data Manager to save!");
        return;
    }

    //Ask the user for a dir to store data
    if (directory.isEmpty()) {
        directory = QFileDialog::getExistingDirectory(
            NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
            QFileDialog::ShowDirsOnly | QFileDialog::DontUseNativeDialog);
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

            bool ok;
            QString tmpFileName = fileName;
            fileName = QInputDialog::getText(NULL, tr("Save As"), tr("File Name:"), QLineEdit::Normal, fileName, &ok);
            if (ok && !fileName.isEmpty() && fileName.endsWith(".nii")) {

                path = directory + "/" + fileName;
                mitk::IOUtil::Save(image, path.toStdString());

            } else {
                fileName = tmpFileName;
                QMessageBox::warning(NULL, "Attention", "Please type a file name with the right extension (i.e. .nii)!");
                return;
            }//_fileName

        } else {
            QMessageBox::warning(NULL, "Attention", "Please select a segmentation image from the Data Manager!");
            return;
        }//_image
    } else
        return;
}

void YZSegView::ScarSeg() {

    //Set buttons visibility
    if (m_Controls.button_4_1->isVisible()) {
        m_Controls.button_4_1->setVisible(false);
        m_Controls.button_4_2->setVisible(false);
        m_Controls.button_4_3->setVisible(false);
        m_Controls.label_1->setVisible(false);
        m_Controls.label_2->setVisible(false);
        m_Controls.label_3->setVisible(false);
    } else {
        m_Controls.button_4_1->setVisible(true);
        m_Controls.button_4_2->setVisible(true);
        m_Controls.button_4_3->setVisible(true);
        m_Controls.label_1->setVisible(true);
        m_Controls.label_2->setVisible(true);
        m_Controls.label_3->setVisible(true);
    }//_if
}

void YZSegView::ScarSeg_FWHM() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) {
        QMessageBox::warning(NULL, "Attention", "Please select an image from the Data Manager to add landmarks!");
        return;
    }

    mitk::DataNode::Pointer imgNode = nodes.at(0);
    mitk::DataNode::Pointer segNode = nodes.at(1);
    mitk::BaseData::Pointer imgdata = imgNode->GetData();
    mitk::BaseData::Pointer segdata = segNode->GetData();

    //if the selected data are images
    if (imgdata && segdata) {

        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(imgdata.GetPointer());
        mitk::Image::Pointer seg = dynamic_cast<mitk::Image*>(segdata.GetPointer());

        if (image && seg) {

            typedef itk::Image<short, 3> ImageType;
            ImageType::Pointer itkImage = ImageType::New();
            mitk::CastToItkImage(image, itkImage);

            ImageType::Pointer itkSeg = ImageType::New();
            mitk::CastToItkImage(seg, itkSeg);

            //Create an empty image LV of the same size as the input images
            ImageType::Pointer LV = ImageType::New();
            ImageType::IndexType start;
            start[0] = 0;
            start[1] = 0;
            start[2] = 0;

            ImageType::SizeType size = itkImage->GetRequestedRegion().GetSize();
            ImageType::RegionType region;
            region.SetSize(size);
            region.SetIndex(start);

            LV->SetRegions(region);
            LV->Allocate();

            //Segment the LV from the input image
            ConstIteratorType inIm(itkImage, itkImage->GetRequestedRegion());
            ConstIteratorType inSeg(itkSeg, itkSeg->GetRequestedRegion());
            IteratorType outLV(LV, LV->GetRequestedRegion());

            //Loop over the image
            for (inIm.GoToBegin(), inSeg.GoToBegin(), outLV.GoToBegin(); !inIm.IsAtEnd(); ++inIm, ++inSeg, ++outLV) {

                //if we are on the left ventricular wall
                if (inSeg.Get() != 0) {
                    outLV.Set(inIm.Get());
                } else {
                    //Set all the values not on the LV wall to be 0
                    outLV.Set(0);
                }//_if
            }

            //Find the maximum infart value within LV
            double max_inf = FWMH_getPeak_SORT(LV);

            //Apply the threshold
            thresholdImage(itkImage, LV, 0.5 * max_inf); //, max_inf, 1); // based on Schmidt et al. that SICore > 0.5 x Peak-infarct
        }//_if_image

    } else {

        QMessageBox::warning(NULL, "Attention", "Please select an image from the Data Manager to add landmarks!");
        return;
    }//_if_data
}

void YZSegView::ScarSeg_SD() {

    //Set visibility
    m_Controls.label_1->setVisible(false);
    m_Controls.label_2->setVisible(false);
    m_Controls.label_3->setVisible(false);
    m_Controls.label_4->setVisible(true);
    m_Controls.label_5->setVisible(true);
    m_Controls.label_6->setVisible(true);
    m_Controls.button_4_2_1->setVisible(true);
    m_Controls.button_4_2_2->setVisible(true);
    m_Controls.button_4_2_3->setVisible(true);

    //Select a small region of remote left ventricular myocardium
    int reply = QMessageBox::question(
        NULL, "Question", "Do you have a remote myocardium segmentation to load?", QMessageBox::Yes, QMessageBox::No);

    if (reply == QMessageBox::Yes) {

        QString path = QFileDialog::getOpenFileName(
            NULL, "Open remote myocardium Segmentation file", mitk::IOUtil::GetProgramPath().c_str(), QmitkIOUtil::GetFileOpenFilterString());
        if (path.isEmpty()) return;
        mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
        mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());

    } else {

        //Show the plugin
        this->GetSite()->GetPage()->ShowView("org.mitk.views.segmentation");
    }//_if
}

void YZSegView::ScarSeg_4SD() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 3) {
        QMessageBox::warning(
            NULL, "Attention",
            "Please select the LGE image, LV segmentation and the remote myocardium in this order from the Data Manager!");
        return;
    }

    mitk::DataNode::Pointer imgNode = nodes.at(0);
    mitk::DataNode::Pointer segNode = nodes.at(1);
    mitk::DataNode::Pointer myoNode = nodes.at(2);
    mitk::BaseData::Pointer imgdata = imgNode->GetData();
    mitk::BaseData::Pointer segdata = segNode->GetData();
    mitk::BaseData::Pointer myodata = myoNode->GetData();

    //if the selected data are images
    if (imgdata && segdata && myodata) {

        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(imgdata.GetPointer());
        mitk::Image::Pointer seg = dynamic_cast<mitk::Image*>(segdata.GetPointer());
        mitk::Image::Pointer myo = dynamic_cast<mitk::Image*>(myodata.GetPointer());

        if (image && seg && myo) {

            ImageType::Pointer itkImage = ImageType::New();
            mitk::CastToItkImage(image, itkImage);

            ImageType::Pointer itkSeg = ImageType::New();
            mitk::CastToItkImage(seg, itkSeg);

            ImageType::Pointer itkMyo = ImageType::New();
            mitk::CastToItkImage(myo, itkMyo);

            //Create an empty image LV of the same size as the input images
            ImageType::Pointer LV = ImageType::New();
            ImageType::IndexType start;
            start[0] = 0;
            start[1] = 0;
            start[2] = 0;

            ImageType::SizeType size = itkImage->GetRequestedRegion().GetSize();
            ImageType::RegionType region;
            region.SetSize(size);
            region.SetIndex(start);

            LV->SetRegions(region);
            LV->Allocate();

            //Segment the LV from the input image
            ConstIteratorType inIm(itkImage, itkImage->GetRequestedRegion());
            ConstIteratorType inSeg(itkSeg, itkSeg->GetRequestedRegion());
            IteratorType outLV(LV, LV->GetRequestedRegion());

            //Loop over the image
            for (inIm.GoToBegin(), inSeg.GoToBegin(), outLV.GoToBegin(); !inIm.IsAtEnd(); ++inIm, ++inSeg, ++outLV) {

                //if we are on the left ventricular wall
                if (inSeg.Get() != 0) {
                    outLV.Set(inIm.Get());
                } else {
                    //Set all the values not on the LV wall to be 0
                    outLV.Set(0);
                }
            }//_for

            double mean = GetStats(itkImage, itkMyo, 1);
            double std = GetStats(itkImage, itkMyo, 2);
            double max = GetStats(itkImage, LV, 3);

            if (max > mean + 4 * std)
                thresholdImage(itkImage, LV, mean + 4 * std);
            else {
                QMessageBox::warning(
                    NULL, "Attention!",
                    "4 Standard Deviation away from the mean is out of range! Please try again with the customised SD!");
                return;
            }
        }//_if_images

    } else {
        QMessageBox::warning(NULL, "Attention", "Please select images from the Data Manager to segment!");
        return;
    }//_if_data
}

void YZSegView::ScarSeg_6SD() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty() || nodes.size() != 3) {
        QMessageBox::warning(
            NULL, "Attention",
            "Please select the LGE image, LV segmentation and the remote myocardium in this order from the Data Manager!");
        return;
    }

    mitk::DataNode::Pointer imgNode = nodes.at(0);
    mitk::DataNode::Pointer segNode = nodes.at(1);
    mitk::DataNode::Pointer myoNode = nodes.at(2);
    mitk::BaseData::Pointer imgdata = imgNode->GetData();
    mitk::BaseData::Pointer segdata = segNode->GetData();
    mitk::BaseData::Pointer myodata = myoNode->GetData();

    //if the selected data are images
    if (imgdata && segdata && myodata) {

        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(imgdata.GetPointer());
        mitk::Image::Pointer seg = dynamic_cast<mitk::Image*>(segdata.GetPointer());
        mitk::Image::Pointer myo = dynamic_cast<mitk::Image*>(myodata.GetPointer());

        if (image && seg && myo) {

            ImageType::Pointer itkImage = ImageType::New();
            mitk::CastToItkImage(image, itkImage);

            ImageType::Pointer itkSeg = ImageType::New();
            mitk::CastToItkImage(seg, itkSeg);

            ImageType::Pointer itkMyo = ImageType::New();
            mitk::CastToItkImage(myo, itkMyo);

            //Create an empty image LV of the same size as the input images
            ImageType::Pointer LV = ImageType::New();
            ImageType::IndexType start;
            start[0] = 0;
            start[1] = 0;
            start[2] = 0;

            ImageType::SizeType size = itkImage->GetRequestedRegion().GetSize();
            ImageType::RegionType region;
            region.SetSize(size);
            region.SetIndex(start);

            LV->SetRegions(region);
            LV->Allocate();

            //Segment the LV from the input image
            ConstIteratorType inIm(itkImage, itkImage->GetRequestedRegion());
            ConstIteratorType inSeg(itkSeg, itkSeg->GetRequestedRegion());
            IteratorType outLV(LV, LV->GetRequestedRegion());

            //Loop over the image
            for (inIm.GoToBegin(), inSeg.GoToBegin(), outLV.GoToBegin(); !inIm.IsAtEnd(); ++inIm, ++inSeg, ++outLV) {
                //if we are on the left ventricular wall
                if (inSeg.Get() != 0) {
                    outLV.Set(inIm.Get());
                } else {
                    //Set all the values not on the LV wall to be 0
                    outLV.Set(0);
                }
            }//_for

            double mean = GetStats(itkImage, itkMyo, 1);
            double std = GetStats(itkImage, itkMyo, 2);
            double max = GetStats(itkImage, LV, 3);

            if (max > mean + 6 * std)
                thresholdImage(itkImage, LV, mean + 6 * std);
            else {
                QMessageBox::warning(
                    NULL, "Attention!",
                    "6 Standard Deviation away from the mean is out of range! Please try again with the customised SD!");
                return;
            }
        }//_if_image

    } else {
        QMessageBox::warning(NULL, "Attention", "Please select an image from the Data Manager to add landmarks!");
        return;
    }//_if_data
}

void YZSegView::ScarSeg_CustomisedSD() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty() || nodes.size() != 3) {
        QMessageBox::warning(
            NULL, "Attention",
            "Please select the LGE image, LV segmentation and the remote myocardium in this order from the Data Manager!");
        return;
    }

    mitk::DataNode::Pointer imgNode = nodes.at(0);
    mitk::DataNode::Pointer segNode = nodes.at(1);
    mitk::DataNode::Pointer myoNode = nodes.at(2);
    mitk::BaseData::Pointer imgdata = imgNode->GetData();
    mitk::BaseData::Pointer segdata = segNode->GetData();
    mitk::BaseData::Pointer myodata = myoNode->GetData();

    //if the selected data are images
    if (imgdata && segdata && myodata) {

        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(imgdata.GetPointer());
        mitk::Image::Pointer seg = dynamic_cast<mitk::Image*>(segdata.GetPointer());
        mitk::Image::Pointer myo = dynamic_cast<mitk::Image*>(myodata.GetPointer());

        if (image && seg && myo) {

            ImageType::Pointer itkImage = ImageType::New();
            mitk::CastToItkImage(image, itkImage);

            ImageType::Pointer itkSeg = ImageType::New();
            mitk::CastToItkImage(seg, itkSeg);

            ImageType::Pointer itkMyo = ImageType::New();
            mitk::CastToItkImage(myo, itkMyo);

            //Create an empty image LV of the same size as the input images
            ImageType::Pointer LV = ImageType::New();
            ImageType::IndexType start;
            start[0] = 0;
            start[1] = 0;
            start[2] = 0;

            ImageType::SizeType size = itkImage->GetRequestedRegion().GetSize();
            ImageType::RegionType region;
            region.SetSize(size);
            region.SetIndex(start);

            LV->SetRegions(region);
            LV->Allocate();

            //Segment the LV from the input image
            ConstIteratorType inIm(itkImage, itkImage->GetRequestedRegion());
            ConstIteratorType inSeg(itkSeg, itkSeg->GetRequestedRegion());
            IteratorType outLV(LV, LV->GetRequestedRegion());

            //Loop over the image
            for (inIm.GoToBegin(), inSeg.GoToBegin(), outLV.GoToBegin(); !inIm.IsAtEnd(); ++inIm, ++inSeg, ++outLV) {

                //if we are on the left ventricular wall
                if (inSeg.Get() != 0) {
                    outLV.Set(inIm.Get());
                } else {
                    //Set all the values not on the LV wall to be 0
                    outLV.Set(0);
                }
            }//_for

            double mean = GetStats(itkImage, itkMyo, 1);
            double std = GetStats(itkImage, itkMyo, 2);
            double max = GetStats(itkImage, LV, 3);

            bool ok;
            QString msg = "Please enter a whole number to set the threshold ";
            int num = QInputDialog::getInt(NULL, tr("Scar Segmentation: SD method: "), msg, 0, 0, 10, 1, &ok);
            if (ok) {
                if (max > mean + num * std)
                    thresholdImage(itkImage, LV, mean + num * std);
                else {
                    QMessageBox::warning(
                        NULL, "Attention!",
                        QString::number(num) + " Standard Deviation away from the mean is out of range! Please try again with the customised SD!");
                    return;
                }
            }
        }//_if_images

    } else {
        QMessageBox::warning(NULL, "Attention", "Please select an image from the Data Manager!");
        return;
    }//_if_data
}

void YZSegView::ScarSeg_save() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1) {
        QMessageBox::warning(NULL, "Attention", "Please select a segmentation from the Data Manager to save!");
        return;
    }

    //Ask the user for a dir to store data
    if (directory.isEmpty()) {
        directory = QFileDialog::getExistingDirectory(
            NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(),
            QFileDialog::ShowDirsOnly | QFileDialog::DontUseNativeDialog);
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

            bool ok;
            QString tmpFileName = fileName;
            fileName = QInputDialog::getText(NULL, tr("Save As"), tr("File Name:"), QLineEdit::Normal, fileName, &ok);
            if (ok && !fileName.isEmpty() && fileName.endsWith(".nii")) {

                path = directory + "/" + fileName;
                mitk::IOUtil::Save(image, path.toStdString());

            } else {
                fileName = tmpFileName;
                QMessageBox::warning(NULL, "Attention", "Please type a file name with the right extension (i.e. .nii)!");
                return;
            }//_fileName

        } else {
            QMessageBox::warning(NULL, "Attention", "Please select a segmentation image from the Data Manager!");
            return;
        }//_image
    } else
        return;
}

/**
 * @brief This function returns the average top 5% of the maximum intensity value in the left ventricular myocardium
 * @param x
 * @return
 */
double YZSegView::FWMH_getPeak_SORT(ImageType::Pointer x) {

    //Copy all the intensity values to a temporary 1D array arr, which will be sorted later
    int m = 0;
    ConstIteratorType in(x, x->GetRequestedRegion());

    //Find the number of non-zero elements in the image
    for (in.GoToBegin(); !in.IsAtEnd(); ++in) {
        if (in.Get() != 0)
            m++;
    }

    //Assign the non-zero elements to the new array arr
    int n = 0;
    int * arr = new int[m];
    for (in.GoToBegin(); !in.IsAtEnd(); ++in) {
        if (in.Get() != 0) {
            arr[n] = in.Get();
            n++;
        }
    }

    //Sort the temporary array arr in the ascending order
    insertion_sort(arr, m);

    //Find the average value of the top 5% of the intensity values
    double max = 0;
    int start = m * 0.998;
    for (int i = start - 1; i < m; i++) {
        max += arr[i];
    }
    delete[] arr;
    max = max / (m - start);
    return max;
}

void YZSegView::thresholdImage(ImageType::Pointer img, ImageType::Pointer mask, double threshold) {

    //Create two iterators for the image and the mask
    IteratorType in_img(img, img->GetRequestedRegion());
    ConstIteratorType in_mask(mask, mask->GetRequestedRegion());

    //Loop over the image
    for (in_img.GoToBegin(), in_mask.GoToBegin(); !in_img.IsAtEnd(); ++in_img, ++in_mask) {

        //if inside the mask
        if (in_img.Get() >= threshold && in_mask.Get() != 0) {
            in_img.Set(256);
        } else {
            //Set all the voxels outside the mask to be 0
            in_img.Set(0);
        }
    }//_for
}

void YZSegView::insertion_sort(int array[], int l) {

    for (int i = 0; i < l; i++) {
        int j = i;
        while (j > 0 && array[j] < array[j - 1]) {
            int temp = array[j];
            array[j] = array[j - 1];
            array[j - 1] = temp;
            j--;
        }
    }//_for
}

double YZSegView::GetStats(ImageType::Pointer im, ImageType::Pointer mask, int a) {

    int n = 0;
    double sum = 0;
    double diff = 0;
    double max = -1;

    //Loop over the selected region of myocardium segmentation
    IteratorType in_im(im, im->GetRequestedRegion());
    ConstIteratorType in_mask(mask, mask->GetRequestedRegion());

    for (in_im.GoToBegin(), in_mask.GoToBegin(); !in_im.IsAtEnd(); ++in_im, ++in_mask) {
        if (in_mask.Get() != 0) {
            if (max <= in_im.Get())
                max = in_im.Get();
            sum += in_im.Get();
            n++;
        }
    }//_for

    double mean = sum / n;
    for (in_im.GoToBegin(), in_mask.GoToBegin(); !in_im.IsAtEnd(); ++in_im, ++in_mask) {
        if (in_mask.Get() != 0)
            diff += (in_im.Get() - mean) * (in_im.Get() - mean);
    }

    double sd = sqrt(diff / n);
    if (a == 1)
        return mean;
    else if (a == 2)
        return sd;
    else if (a == 3)
        return max;
    else {
        QMessageBox::warning(NULL, "Attention", "error in function GetStats");
        return -1;
    }
}
