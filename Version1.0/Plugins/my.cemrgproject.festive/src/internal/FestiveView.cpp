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
 * Fetal Motion Corrected Slice to Volume Registration (Festive) Plugin
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
#include <berryIWorkbenchWindow.h>
#include <berryIWorkbenchPage.h>
#include <berryFileEditorInput.h>

// Qmitk
#include <QmitkIOUtil.h>
#include <mitkIOUtil.h>
#include <mitkProgressBar.h>
#include <mitkIDataStorageService.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>
#include <mitkDataStorageEditorInput.h>
#include <mitkDataStorage.h>
#include "mitkPluginActivator.h"
#include "FestiveView.h"

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QInputDialog>

// MyCemrgLib
#include <CemrgCommandLine.h>


const std::string FestiveView::VIEW_ID = "my.cemrgproject.views.festiveview";

void FestiveView::SetFocus() {

    m_Controls.button_1->setFocus();
}

void FestiveView::CreateQtPartControl( QWidget *parent ) {

    //create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_1, SIGNAL(clicked()), this, SLOT(LoadDICOM()));
    connect(m_Controls.button_2, SIGNAL(clicked()), this, SLOT(ProcessIMGS()));
    connect(m_Controls.button_3, SIGNAL(clicked()), this, SLOT(MaskIMGS()));
    connect(m_Controls.button_4, SIGNAL(clicked()), this, SLOT(EstablishConnection()));
    connect(m_Controls.button_5, SIGNAL(clicked()), this, SLOT(CopyServer()));
    connect(m_Controls.button_6, SIGNAL(clicked()), this, SLOT(Reconstruction()));
    connect(m_Controls.button_7, SIGNAL(clicked()), this, SLOT(Download()));
    connect(m_Controls.button_r, SIGNAL(clicked()), this, SLOT(Reset()));

    //Set visibility of buttons
    m_Controls.button_2_1->setVisible(false);
    m_Controls.button_3_1->setVisible(false);
    connect(m_Controls.button_2_1, SIGNAL(clicked()), this, SLOT(ConvertNII()));
    connect(m_Controls.button_3_1, SIGNAL(clicked()), this, SLOT(SaveMASK()));
}

void FestiveView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*source*/,
                                     const QList<mitk::DataNode::Pointer>& /*nodes*/) {

}

void FestiveView::LoadDICOM() {

    //Use MITK DICOM editor
    QString editor_id = "org.mitk.editors.dicomeditor";
    berry::IEditorInput::Pointer input(new berry::FileEditorInput(QString()));
    this->GetSite()->GetPage()->OpenEditor(input, editor_id);
}

void FestiveView::ProcessIMGS() {

    //Toggle visibility of buttons
    if (m_Controls.button_2_1->isVisible()) {
        m_Controls.button_2_1->setVisible(false);
    } else {
        m_Controls.button_2_1->setVisible(true);
    }
}

void FestiveView::ConvertNII() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() < 3) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please load and select at least 3 images from the Data Manager before starting this step!");
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

    //Convert to Nifti
    int ctr = 0;
    QString path;
    this->BusyCursorOn();
    mitk::ProgressBar::GetInstance()->AddStepsToDo(nodes.size());
    foreach (mitk::DataNode::Pointer node, nodes) {
        mitk::BaseData::Pointer data = node->GetData();
        mitk::ProgressBar::GetInstance()->Progress();
        if (data) {
            //Test if this data item is an image
            mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
            if (image) {
                path = directory + mitk::IOUtil::GetDirectorySeparator() + "dcm-" + QString::number(ctr++) + ".nii";
                mitk::IOUtil::Save(image, path.toStdString());
                this->GetDataStorage()->Remove(node);
                mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
            } else
                return;
        } else
            return;
    }//_for

    nodes.clear();
    this->BusyCursorOff();
    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());
}

void FestiveView::MaskIMGS() {

    //Toggle visibility of buttons
    if (m_Controls.button_3_1->isVisible()) {
        m_Controls.button_3_1->setVisible(false);
        return;
    } else {
        m_Controls.button_3_1->setVisible(true);
    }

    int reply = QMessageBox::question(
                NULL, "Question", "Do you have a mask to load?", QMessageBox::Yes, QMessageBox::No);

    if (reply == QMessageBox::Yes) {

        QString path = QFileDialog::getOpenFileName(
                    NULL, "Open the mask file", mitk::IOUtil::GetProgramPath().c_str(), QmitkIOUtil::GetFileOpenFilterString());
        if (path.isEmpty()) return;
        mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage());
        mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());

    } else {

        //Show the plugin
        this->GetSite()->GetPage()->ShowView("org.mitk.views.segmentation");
    }//_if
}

void FestiveView::SaveMASK() {

    //Check for selection of images
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please select a mask from the Data Manager to save!");
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
    mitk::DataNode::Pointer mskNode = nodes.at(0);
    mitk::BaseData::Pointer data = mskNode->GetData();
    if (data) {
        //Test if this data item is an image
        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
        if (image) {

            bool ok;
            QString fileName = "Mask.nii";
            fileName = QInputDialog::getText(NULL, tr("Save As"), tr("File Name:"), QLineEdit::Normal, fileName, &ok);
            if (ok && !fileName.isEmpty() && fileName.endsWith(".nii")) {

                path = directory + mitk::IOUtil::GetDirectorySeparator() + fileName;
                mitk::IOUtil::Save(image, path.toStdString());

            } else {
                QMessageBox::warning(NULL, "Attention", "Please type a file name with the right extension (i.e. .nii)!");
                return;
            }//_fileName

        } else {
            QMessageBox::warning(NULL, "Attention", "Please select a mask image from the Data Manager!");
            return;
        }//_image
    } else
        return;
}

void FestiveView::EstablishConnection() {

    //Prompt for username and password
    QDialog* inputs = new QDialog(0,0);
    m_UILogin.setupUi(inputs);
    connect(m_UILogin.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_UILogin.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    int dialogCode = inputs->exec();

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {

        userID = m_UILogin.lineEdit_1->text();
        server = m_UILogin.lineEdit_2->text();

        //Set default values
        if (userID.isEmpty()) {
            QMessageBox::warning(NULL, "Attention", "Please enter a valid username!");
            return;
        }//_if
        if (server.isEmpty()) {
            QMessageBox::warning(NULL, "Attention", "Reverting to default server!");
            server = "gpubeastie03-pc";
        }//_if

        this->BusyCursorOn();
        mitk::ProgressBar::GetInstance()->AddStepsToDo(1);

        cmd = std::unique_ptr<CemrgCommandLine>(new CemrgCommandLine());
        if (cmd) {
            if (cmd->ConnectToServer(userID,server)) {
                QMessageBox::information(NULL, "Attention", "Connection has been established!");
                mitk::ProgressBar::GetInstance()->Progress();
            } else {
                QMessageBox::critical(NULL, "Attention", "Connection has not been established!");
                mitk::ProgressBar::GetInstance()->Progress();
            }
        }//_if_cmd

        //Clean up
        mitk::ProgressBar::GetInstance()->Progress();
        this->BusyCursorOff();
        inputs->deleteLater();

    } else if(dialogCode == QDialog::Rejected) {
        inputs->close();
        inputs->deleteLater();
    }//_if
}

void FestiveView::CopyServer() {

    //Check for selection of data
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please select all your items from the Data Manager for transfering to the server!");
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

    /*** UNIX based commands ***/
    char dirSep = mitk::IOUtil::GetDirectorySeparator();
    std::system(("rm -rf " + (directory + dirSep + "Transfer").toStdString()).c_str());
    std::system(("mkdir "  + (directory + dirSep + "Transfer").toStdString()).c_str());
    imgsList.clear();

    QString path;
    this->BusyCursorOn();
    mitk::ProgressBar::GetInstance()->AddStepsToDo(nodes.size());
    foreach (mitk::DataNode::Pointer node, nodes) {
        mitk::BaseData::Pointer data = node->GetData();
        mitk::ProgressBar::GetInstance()->Progress();
        if (data) {
            //Test if this data item is an image
            mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(data.GetPointer());
            if (image) {
                if (QString::fromStdString(node->GetName()).contains("Mask", Qt::CaseInsensitive)) node->SetName("Mask");
                imgsList << QString::fromStdString(node->GetName() + ".nii.gz");
                path = directory + dirSep + "Transfer" + dirSep + node->GetName().c_str() + ".nii.gz";
                mitk::IOUtil::Save(image, path.toStdString());
            } else
                return;
        } else
            return;
    }//_for

    //Transfer to server
    if (cmd) {
        mitk::ProgressBar::GetInstance()->AddStepsToDo(1);
        std::unique_ptr<CemrgCommandLine> scpCMD(new CemrgCommandLine());
        if (scpCMD->TransferTFServer(directory, "Transfer", userID, server, false)) {
            QMessageBox::information(NULL, "Attention", "Transfer was successful!");
            mitk::ProgressBar::GetInstance()->Progress();
        } else {
            QMessageBox::critical(NULL, "Attention", "Transfer was not successful!");
            mitk::ProgressBar::GetInstance()->Progress();
        }//_if
    } else {
        QMessageBox::warning(NULL, "Attention", "Server Connection has not been established!");
    }//_if_cmd

    nodes.clear();
    this->BusyCursorOff();
}

void FestiveView::Reconstruction() {

    //Check for selection of data
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.size() != 1) {
        QMessageBox::warning(
                    NULL, "Attention",
                    "Please select the target image from the Data Manager before running the reconstruction step!");
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

    //Prompt for parameters
    QDialog* inputs = new QDialog(0,0);
    m_UIReconst.setupUi(inputs);
    connect(m_UIReconst.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_UIReconst.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    int dialogCode = inputs->exec();

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {

        bool ok1, ok2, ok3;
        double reso = m_UIReconst.lineEdit_1->text().toDouble(&ok1);
        double delt = m_UIReconst.lineEdit_2->text().toDouble(&ok2);
        int package = m_UIReconst.lineEdit_3->text().toInt(&ok3);
        QString out = m_UIReconst.lineEdit_4->text();

        //Set default values
        if (!ok1 || !ok2 || !ok3 || out.isEmpty() || !out.endsWith(".nii"))
            QMessageBox::warning(NULL, "Attention", "Reverting to default parameters!");
        if (!ok1) reso = 0.75;
        if (!ok2) delt = 150.0;
        if (!ok3) package = 4;
        outName = (out.isEmpty() || !out.endsWith(".nii")) ? "output.nii" : out;
        //_if

        this->BusyCursorOn();
        mitk::ProgressBar::GetInstance()->AddStepsToDo(1);

        if (cmd) {
            std::unique_ptr<CemrgCommandLine> conCMD(new CemrgCommandLine());
            QString targetImage = QString::fromStdString(nodes.at(0)->GetName() + ".nii.gz");
            conCMD->GPUReconstruction(userID, server, imgsList, targetImage, reso, delt, package, outName);
            QMessageBox::information(NULL, "Attention", "Reconstruction Completed!");
            mitk::ProgressBar::GetInstance()->Progress();
        } else {
            QMessageBox::warning(NULL, "Attention", "Server Connection has not been established!");
            mitk::ProgressBar::GetInstance()->Progress();
        }//_if_cmd

        //Clean up
        this->BusyCursorOff();
        inputs->deleteLater();

    } else if(dialogCode == QDialog::Rejected) {
        inputs->close();
        inputs->deleteLater();
    }//_if
}

void FestiveView::Download() {

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

    //File to download
    bool ok;
    outName = QInputDialog::getText(NULL, tr("To Download"), tr("File Name:"), QLineEdit::Normal, outName, &ok);
    if (!ok || outName.isEmpty() || !outName.endsWith(".nii"))
        QMessageBox::warning(NULL, "Attention", "Reverting to default file name!");
    outName = (outName.isEmpty() || !outName.endsWith(".nii")) ? "output.nii" : outName;

    //Transfer from server
    this->BusyCursorOn();
    if (cmd) {
        mitk::ProgressBar::GetInstance()->AddStepsToDo(1);
        std::unique_ptr<CemrgCommandLine> scpCMD(new CemrgCommandLine());
        if (scpCMD->TransferTFServer(directory, outName, userID, server, true)) {

            QString path = directory + mitk::IOUtil::GetDirectorySeparator() + outName;
            if (mitk::IOUtil::Load(path.toStdString(), *this->GetDataStorage()).IsNull()) {
                QMessageBox::critical(NULL, "Attention", "Transfer was not successful!");
                mitk::ProgressBar::GetInstance()->Progress();
                return;
            }//_if
            QMessageBox::information(NULL, "Attention", "Transfer was successful!");
            mitk::ProgressBar::GetInstance()->Progress();

        } else {
            QMessageBox::critical(NULL, "Attention", "Transfer was not successful!");
            mitk::ProgressBar::GetInstance()->Progress();
        }//_if
    } else {
        QMessageBox::warning(NULL, "Attention", "Server Connection has not been established!");
    }//_if_cmd
    this->BusyCursorOff();
}

void FestiveView::Reset() {

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
    directory.clear();
    cmd.reset();
}
