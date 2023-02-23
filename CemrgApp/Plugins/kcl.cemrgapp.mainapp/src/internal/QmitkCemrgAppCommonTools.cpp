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
 * CemrgApp Common Tools
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
#include <berryIWorkbench.h>
#include <berryIWorkbenchWindow.h>
#include <berryIWorkbenchPage.h>

// Qmitk
#include <mitkCoreObjectFactory.h>
#include <mitkProgressBar.h>
#include <QmitkIOUtil.h>

// CemrgAppModule
#include <CemrgCommonUtils.h>
#include <CemrgCommandLine.h>
#include "QmitkCemrgAppCommonTools.h"

// Qt
#include <QMessageBox>
#include <QStringList>
#include <QFileDialog>
#include <QSignalMapper>
#include <QDirIterator>
#include <QFileInfo>

const std::string QmitkCemrgAppCommonTools::VIEW_ID = "org.mitk.views.cemrgappcommontools";

void QmitkCemrgAppCommonTools::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.btn_loadmesh, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::LoadMesh);
    connect(m_Controls.btn_convert2carto, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::ConvertToCarto);
    connect(m_Controls.btn_vtk2cart, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::ConvertCarpToVtk);
    connect(m_Controls.btn_padimage, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::PadImageEdgesWithConstant);
    connect(m_Controls.btn_binariseimage, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::BinariseImage);
    connect(m_Controls.btn_resamplereorient, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::ResampleReorientConvert);
    connect(m_Controls.button_mirtk, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::MirtkOptsSelection);
    connect(m_Controls.button_mirtk_reg, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::MirtkOptsRegister);
    connect(m_Controls.button_mirtk_tx, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::MirtkOptsTransform);
    connect(m_Controls.button_mirtk_invreg, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::MirtkOptsInvRegister);

    m_Controls.button_mirtk_reg->setVisible(false);
    m_Controls.button_mirtk_tx->setVisible(false);
    m_Controls.button_mirtk_invreg->setVisible(false);
}

void QmitkCemrgAppCommonTools::SetFocus() {
}

void QmitkCemrgAppCommonTools::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*source*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

void QmitkCemrgAppCommonTools::LoadMesh() {
    QString path = QFileDialog::getOpenFileName(NULL, "Open Mesh Data File", QmitkIOUtil::GetFileOpenFilterString());
    CemrgCommonUtils::AddToStorage(CemrgCommonUtils::LoadVTKMesh(path.toStdString()), "Mesh", this->GetDataStorage());
}

void QmitkCemrgAppCommonTools::ConvertToCarto() {

    QString path = QFileDialog::getOpenFileName(NULL, "Open Mesh Data File", QmitkIOUtil::GetFileOpenFilterString());
    if (path.isEmpty() || !path.endsWith(".vtk")) {
        QMessageBox::warning(NULL, "Attention", "Select Correct Input File!");
        return;
    }

    //Find thresholds file
    double fileMeanBP = 0.0;
    double fileStdvBP = 0.0;
    bool threshFileExist = false;
    QFileInfo fullPathInfo(path);
    if (fullPathInfo.dir().exists("prodThresholds.txt")) {
        //Threshold file
        QString tPath = fullPathInfo.absolutePath() + "/prodThresholds.txt";
        ifstream prodFileRead(tPath.toStdString());
        if (prodFileRead.is_open()) {

            std::string line;
            std::vector<std::string> lines;
            while (getline(prodFileRead, line))
                lines.push_back(line);
            fileMeanBP = QString::fromStdString(lines.at(2)).toDouble();
            fileStdvBP = QString::fromStdString(lines.at(3)).toDouble();
            prodFileRead.close();
            threshFileExist = true;

        }//_if
    }//_if

    //Ask for user input to set the parameters
    QDialog* inputs = new QDialog(0, 0);
    m_CartoUIThresholding.setupUi(inputs);
    connect(m_CartoUIThresholding.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_CartoUIThresholding.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    connect(m_CartoUIThresholding.radioButton_1, SIGNAL(toggled(bool)), this, SLOT(ConvertToCartoUITextUpdate()));
    connect(m_CartoUIThresholding.comboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(ConvertToCartoUIUpdate()));
    m_CartoUIThresholding.lineEdit_2->setPlaceholderText(QString::number(fileMeanBP));
    m_CartoUIThresholding.lineEdit_3->setPlaceholderText(QString::number(fileStdvBP));
    int dialogCode = inputs->exec();

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {

        //Methods
        int methodType = m_CartoUIThresholding.radioButton_1->isChecked() ? 1 : 2;
        bool discreteScheme = m_CartoUIThresholding.comboBox->currentIndex() == 0 ? true : false;

        //Thresholds
        bool ok0;
        std::vector<double> thresholds;
        QRegExp separator("(\\ |\\,|\\;|\\:|\\t)");
        QStringList thresholdsInput = m_CartoUIThresholding.lineEdit_1->text().trimmed().split(separator);
        for (QString item : thresholdsInput)
            if (item.isEmpty())
                thresholdsInput.removeOne(item);
        if (discreteScheme && thresholdsInput.size() == 0) {
            QMessageBox::warning(NULL, "Attention", "Reverting to default threshold values!");
            thresholds.push_back(methodType == 1 ? 1.20 : 3.0);
            thresholds.push_back(methodType == 1 ? 1.32 : 4.0);
        } else if (discreteScheme && thresholdsInput.count() > 2) {
            QMessageBox::warning(NULL, "Attention", "Parsing thresholds failed!\nReverting to default values.");
            thresholds.push_back(methodType == 1 ? 1.20 : 3.0);
            thresholds.push_back(methodType == 1 ? 1.32 : 4.0);
        } else {
            for (QString item : thresholdsInput) {
                double thresh = item.toDouble(&ok0);
                if (!ok0) {
                    QMessageBox::warning(NULL, "Attention", "Parsing thresholds failed!\nReverting to default values.");
                    thresholds.clear();
                    thresholds.push_back(methodType == 1 ? 1.20 : 3.0);
                    thresholds.push_back(methodType == 1 ? 1.32 : 4.0);
                    break;
                } else
                    thresholds.push_back(thresh);
            }//_for
        }//_if

        //BP values
        bool ok1, ok2;
        double meanBP = m_CartoUIThresholding.lineEdit_2->text().toDouble(&ok1);
        double stdvBP = m_CartoUIThresholding.lineEdit_3->text().toDouble(&ok2);
        if (discreteScheme && (!ok1 || !ok2)) {
            if (threshFileExist) {
                meanBP = fileMeanBP;
                stdvBP = fileStdvBP;
                QMessageBox::information(NULL, "Attention", "Blood pool intensity values from the file was successfully restored.");
            } else {
                QMessageBox::warning(NULL, "Attention", "Parsing blood pool intensity values failed! Try again please.");
                return;
            }//_if
        }//_if

        //Conversion
        bool result = CemrgCommonUtils::ConvertToCarto(path.toStdString(), thresholds, meanBP, stdvBP, methodType, discreteScheme);
        if (result)
            QMessageBox::information(NULL, "Attention", "Conversion Completed!");
        else
            QMessageBox::information(NULL, "Attention", "Conversion Failed!");
        inputs->deleteLater();

    } else if (dialogCode == QDialog::Rejected) {
        inputs->close();
        inputs->deleteLater();
    }//_if
}

void QmitkCemrgAppCommonTools::ConvertToCartoUIUpdate() {

    if (m_CartoUIThresholding.comboBox->currentIndex() == 1) {
        m_CartoUIThresholding.lineEdit_1->setEnabled(false);
        m_CartoUIThresholding.lineEdit_2->setEnabled(false);
        m_CartoUIThresholding.lineEdit_3->setEnabled(false);
        m_CartoUIThresholding.radioButton_1->setEnabled(false);
        m_CartoUIThresholding.radioButton_2->setEnabled(false);
    } else {
        m_CartoUIThresholding.lineEdit_1->setEnabled(true);
        m_CartoUIThresholding.lineEdit_2->setEnabled(true);
        m_CartoUIThresholding.lineEdit_3->setEnabled(true);
        m_CartoUIThresholding.radioButton_1->setEnabled(true);
        m_CartoUIThresholding.radioButton_2->setEnabled(true);
    }//_if
}

void QmitkCemrgAppCommonTools::ConvertToCartoUITextUpdate() {

    if (m_CartoUIThresholding.radioButton_1->isChecked())
        m_CartoUIThresholding.lineEdit_1->setPlaceholderText("1.2; 1.32");
    else
        m_CartoUIThresholding.lineEdit_1->setPlaceholderText("3; 4");
}

void QmitkCemrgAppCommonTools::ConvertCarpToVtk() {
    QString pathElem = QFileDialog::getOpenFileName(NULL, "Open Mesh .elem File");
    if (pathElem.isEmpty() || !pathElem.endsWith(".elem")) {
        QMessageBox::warning(NULL, "Attention", "Select Correct Input (.elem) File!");
        return;
    }
    QFileInfo fi(pathElem);
    QString dir = fi.absolutePath();
    QString vtkPath = dir + "/" + fi.baseName() + ".vtk";

    QString pathPts = QFileDialog::getOpenFileName(NULL, "Open Mesh .pts File", dir.toStdString().c_str());
    if (pathPts.isEmpty() || !pathPts.endsWith(".pts")) {
        QMessageBox::warning(NULL, "Attention", "Select Correct Input (.pts) File!");
        return;
    }

    int regionScalarsReply = QMessageBox::question(NULL, "Question", "Include region as (cell) scalar field?", QMessageBox::Yes, QMessageBox::No);
    CemrgCommonUtils::CarpToVtk(pathElem, pathPts, vtkPath, (regionScalarsReply == QMessageBox::Yes));

    int appendScalarFieldReply = QMessageBox::question(NULL, "Question", "Append a scalar field from a file?", QMessageBox::Yes, QMessageBox::No);
    if (appendScalarFieldReply == QMessageBox::Yes) {
        QString typeData = "";
        int nElem = CemrgCommonUtils::GetTotalFromCarpFile(pathElem);
        int nPts = CemrgCommonUtils::GetTotalFromCarpFile(pathPts);
        int countFields = 0;

        while (appendScalarFieldReply == QMessageBox::Yes) {
            QString path = QFileDialog::getOpenFileName(NULL, "Open Scalar field (.dat) file", dir.toStdString().c_str());
            QFileInfo fi2(path);
            std::vector<double> field = CemrgCommonUtils::ReadScalarField(path);

            int nField = field.size();
            MITK_INFO << ("FieldSize: " + QString::number(nField)).toStdString();
            if (nField == nElem) {
                typeData = "CELL";
            } else if (nField == nPts) {
                typeData = "POINT";
            } else {
                MITK_INFO << "Inconsistent file size";
                break;
            }
            CemrgCommonUtils::AppendScalarFieldToVtk(vtkPath, fi2.baseName(), typeData, field, (countFields == 0));
            countFields++;
            appendScalarFieldReply = QMessageBox::question(NULL, "Question",
                "Append another scalar field from a file?", QMessageBox::Yes, QMessageBox::No);
        }
    }
}

void QmitkCemrgAppCommonTools::PadImageEdgesWithConstant(){
    QString pathToImage = "";
    pathToImage = QFileDialog::getOpenFileName(NULL, "Open image file");
    if (pathToImage.isEmpty()) {
        QMessageBox::warning(NULL, "Attention", "Select Correct Input (.nii) File!");
        return;
    }

    //Ask for user input to set the parameters
    QDialog* inputs = new QDialog(0,0);
    m_ImagePadding.setupUi(inputs);
    connect(m_ImagePadding.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_ImagePadding.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    int dialogCode = inputs->exec();
    if (dialogCode == QDialog::Accepted) {
        bool ok1, ok2;
        int paddingSize = m_ImagePadding.lineEdit_2->text().toInt(&ok1);
        int constantForPadding = m_ImagePadding.lineEdit_3->text().toDouble(&ok2);
        QString outputName = m_ImagePadding.lineEdit_3->text();
        QString outputPath = pathToImage;

        if(!ok1){
            paddingSize = 2;
        }
        if(!ok2){
            constantForPadding = 0;
        }
        if(!outputName.isEmpty()){
            QFileInfo fi(pathToImage);
            outputPath = fi.absolutePath() + "/" + outputName + fi.suffix();
        }

        CemrgCommonUtils::SavePadImageWithConstant(pathToImage, outputPath, paddingSize, constantForPadding);

        QMessageBox::information(NULL, "Attention", "Operation finished. File created");
    }

}

void QmitkCemrgAppCommonTools::BinariseImage(){
    QString pathToImage = "";
    pathToImage = QFileDialog::getOpenFileName(NULL, "Open image file");
    if (pathToImage.isEmpty()) {
        QMessageBox::warning(NULL, "Attention", "Select Correct Input (.nii) File!");
        return;
    }

    QFileInfo fi(pathToImage);
    QString outPath = fi.absolutePath() + "/" + fi.baseName() + "-bin." + fi.suffix();

    mitk::Image::Pointer im = mitk::IOUtil::Load<mitk::Image>(pathToImage.toStdString());
    mitk::Image::Pointer outIm = CemrgCommonUtils::ReturnBinarised(im);

    mitk::IOUtil::Save(outIm, outPath.toStdString());
}

void QmitkCemrgAppCommonTools::ResampleReorientConvert(){
    QString pathToImage = "";
    pathToImage = QFileDialog::getOpenFileName(NULL, "Open image file");
    if (pathToImage.isEmpty()) {
        QMessageBox::warning(NULL, "Attention", "Incorrect input!");
        return;
    }

    std::string title, msg;
    title = "Choose Image Type";
    msg = "Is this a binary image (i.e a segmentation)?";
    int replyImBinary = QMessageBox::question(NULL, title.c_str(), msg.c_str(), QMessageBox::Yes, QMessageBox::No);

    bool resamplebool=true, reorientbool=true;
    bool isBinary=(replyImBinary==QMessageBox::Yes);
    QFileInfo fi(pathToImage);
    QString pathToOutput=fi.absolutePath() + "/" + fi.baseName() + ".nii";
    bool success = CemrgCommonUtils::ImageConvertFormat(pathToImage, pathToOutput, resamplebool, reorientbool, isBinary);

    if(success){
        title = "Attention";
        msg = "Image resampled, reoriented and converted to NIFTI";
        QMessageBox::information(NULL, title.c_str(), msg.c_str());
    }
}

void QmitkCemrgAppCommonTools::MirtkOptsSelection(){
    if (m_Controls.button_mirtk_reg->isVisible()){
        m_Controls.button_mirtk_reg->setVisible(false);
        m_Controls.button_mirtk_tx->setVisible(false);
        m_Controls.button_mirtk_invreg->setVisible(false);
    } else {
        m_Controls.button_mirtk_reg->setVisible(true);
        m_Controls.button_mirtk_tx->setVisible(true);
        m_Controls.button_mirtk_invreg->setVisible(true);
    }
}

void QmitkCemrgAppCommonTools::MirtkOptsRegister(){
    QString directory = QFileDialog::getExistingDirectory(NULL, "Open directory with files",
        mitk::IOUtil::GetProgramPath().c_str(), QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);

    QDialog* inputs = new QDialog(0,0);
    QSignalMapper* signalMapper = new QSignalMapper(this);

    m_MirtkUIOptions.setupUi(inputs);
    connect(m_MirtkUIOptions.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_MirtkUIOptions.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    connect(m_MirtkUIOptions.button_browse1, SIGNAL(clicked()), signalMapper, SLOT(map()));
    connect(m_MirtkUIOptions.button_browse2, SIGNAL(clicked()), signalMapper, SLOT(map()));
    signalMapper->setMapping(m_MirtkUIOptions.button_browse1, "1,0,"+directory);
    signalMapper->setMapping(m_MirtkUIOptions.button_browse2, "2,0,"+directory);
    connect(signalMapper, SIGNAL(mapped(QString)), this, SLOT(MirtkOptsBrowse(const QString&)));

    m_MirtkUIOptions.check_tx_points->setVisible(false);
    m_MirtkUIOptions.check_multiple_tx->setVisible(false);
    QString msgInput1, msgInput2, msgOutput;
    msgInput1 = "Select input 1 filename (moving image)";
    msgInput2 = "Select input 2 filename (fixed image)";
    msgOutput = "Output name for DOF file (no extension, default = modelname)";
    m_MirtkUIOptions.lineEdit_input1->setPlaceholderText(msgInput1);
    m_MirtkUIOptions.lineEdit_input2->setPlaceholderText(msgInput2);
    m_MirtkUIOptions.lineEdit_output->setPlaceholderText(msgOutput);
    int dialogCode = inputs->exec();

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {
        QDir home(directory);
        QString moving = home.relativeFilePath(m_MirtkUIOptions.lineEdit_input1->text());
        QString fixed = home.relativeFilePath(m_MirtkUIOptions.lineEdit_input2->text());
        QString model = m_MirtkUIOptions.combo_reg_model->currentText();
        QString outputname = m_MirtkUIOptions.lineEdit_output->text();

        if(moving.isEmpty() || fixed.isEmpty()){
            QMessageBox::information(NULL, "Attention", "Images not selected correctly");
            return;
        }

        if(outputname.isEmpty()){
            int indexOfPlus = model.indexOf("+");
            outputname = model;
            outputname.replace(indexOfPlus, 1, "_");
            MITK_INFO << ("No output name specified, using : " + outputname).toStdString();
        }

        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->ExecuteRegistration(directory, fixed, moving, outputname+".dof", model);
    }
}

void QmitkCemrgAppCommonTools::MirtkOptsTransform(){
    QString directory = QFileDialog::getExistingDirectory(NULL, "Open directory with files",
        mitk::IOUtil::GetProgramPath().c_str(), QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);

    QDialog* inputs = new QDialog(0,0);
    QSignalMapper* signalMapper = new QSignalMapper(this);

    m_MirtkUIOptions.setupUi(inputs);
    connect(m_MirtkUIOptions.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_MirtkUIOptions.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    connect(m_MirtkUIOptions.button_browse1, SIGNAL(clicked()), signalMapper, SLOT(map()));
    connect(m_MirtkUIOptions.button_browse2, SIGNAL(clicked()), signalMapper, SLOT(map()));
    signalMapper->setMapping(m_MirtkUIOptions.button_browse1, "1,1,"+directory);
    signalMapper->setMapping(m_MirtkUIOptions.button_browse2, "2,1,"+directory);
    connect(signalMapper, SIGNAL(mapped(QString)), this, SLOT(MirtkOptsBrowse(const QString&)));

    m_MirtkUIOptions.check_tx_points->setVisible(true);
    m_MirtkUIOptions.check_multiple_tx->setVisible(true);
    m_MirtkUIOptions.label_reg_model->setVisible(false);
    m_MirtkUIOptions.combo_reg_model->setVisible(false);
    QString msgInput1, msgInput2, msgOutput;
    msgInput1 = "Select input 1 filename (image or point set)";
    msgInput2 = "Select input 2 filename (Registration file DOF)";
    msgOutput = "Output image or point set (no extension, default = transformation)";
    m_MirtkUIOptions.lineEdit_input1->setPlaceholderText(msgInput1);
    m_MirtkUIOptions.lineEdit_input2->setPlaceholderText(msgInput2);
    m_MirtkUIOptions.lineEdit_output->setPlaceholderText(msgOutput);
    int dialogCode = inputs->exec();

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {
        QDir home(directory);
        QString objectToTransform = home.relativeFilePath(m_MirtkUIOptions.lineEdit_input1->text());
        QString dofFile = home.relativeFilePath(m_MirtkUIOptions.lineEdit_input2->text());
        bool transformPoints = m_MirtkUIOptions.check_tx_points->isChecked();
        bool multipleTransformations = m_MirtkUIOptions.check_multiple_tx->isChecked();
        QString outputname = m_MirtkUIOptions.lineEdit_output->text();

        MITK_INFO << ("objectToTransform = " + objectToTransform).toStdString();

        if(objectToTransform.isEmpty() || dofFile.isEmpty()){
            QMessageBox::warning(NULL, "Attention", "Inputs not selected correctly");
            return;
        }

        if(outputname.isEmpty()){
            QFileInfo fi1(objectToTransform);
            QFileInfo fi2(dofFile);

            outputname = "tx_" + fi1.baseName() + "-by-" + fi2.baseName();
            MITK_INFO << ("No output name specified, using : " + outputname).toStdString();
        }

        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        if(transformPoints){
            QStringList objects, outputFiles;
            if(multipleTransformations){
                QFileInfo fi(directory+"/"+objectToTransform);
                QString base = fi.baseName().left(fi.baseName().size()-1);
                QDir objectdirect(fi.absolutePath());
                MITK_INFO << ("Search dir: " + fi.absolutePath()).toStdString();
                QStringList dirEntryList = objectdirect.entryList();

                for (int jx = 0; jx < dirEntryList.size(); jx++) {
                    QString thisFile = home.relativeFilePath(fi.absolutePath()+"/"+dirEntryList.at(jx));

                    if(thisFile.contains(base) && thisFile.contains(".vtk")){
                        QFileInfo thisFi(thisFile);

                        objects.push_back(thisFile);
                        outputFiles.push_back("out_"+thisFi.baseName()+".vtk");
                    }
                }
            } else{
                objects.push_back(objectToTransform);
                outputFiles.push_back(outputname+".vtk");
            }

            MITK_INFO << objects.size();
            for (int ix = 0; ix < objects.size(); ix++) {
                MITK_INFO << ("Input: " + objects.at(ix)).toStdString();
                MITK_INFO << ("Output: " + outputFiles.at(ix)).toStdString();
                cmd->ExecuteTransformationOnPoints(directory, objects.at(ix), outputFiles.at(ix), dofFile);
            }

        } else {
            cmd->ExecuteTransformation(directory, objectToTransform, outputname+".nii", dofFile);
        }
    }
}

void QmitkCemrgAppCommonTools::MirtkOptsInvRegister(){
    QString directory = QFileDialog::getExistingDirectory(NULL, "Open directory with files",
        mitk::IOUtil::GetProgramPath().c_str(), QFileDialog::ShowDirsOnly|QFileDialog::DontUseNativeDialog);

    QDialog* inputs = new QDialog(0,0);
    QSignalMapper* signalMapper = new QSignalMapper(this);

    m_MirtkUIOptions.setupUi(inputs);
    connect(m_MirtkUIOptions.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
    connect(m_MirtkUIOptions.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));
    connect(m_MirtkUIOptions.button_browse1, SIGNAL(clicked()), signalMapper, SLOT(map()));
    signalMapper->setMapping(m_MirtkUIOptions.button_browse1, "1,2,"+directory);
    connect(signalMapper, SIGNAL(mapped(QString)), this, SLOT(MirtkOptsBrowse(const QString&)));

    m_MirtkUIOptions.check_tx_points->setVisible(false);
    m_MirtkUIOptions.check_multiple_tx->setVisible(false);
    m_MirtkUIOptions.lineEdit_input2->setVisible(false);
    m_MirtkUIOptions.button_browse2->setVisible(false);
    m_MirtkUIOptions.label_reg_model->setVisible(false);
    m_MirtkUIOptions.combo_reg_model->setVisible(false);

    QString msgInput1, msgOutput;
    msgInput1 = "Select input 1 filename (dof file)";
    msgOutput = "Output name for DOF file (no extension, default = inverse_inputname)";
    m_MirtkUIOptions.lineEdit_input1->setPlaceholderText(msgInput1);
    m_MirtkUIOptions.lineEdit_output->setPlaceholderText(msgOutput);
    int dialogCode = inputs->exec();

    //Act on dialog return code
    if (dialogCode == QDialog::Accepted) {
        QString inputDof = m_MirtkUIOptions.lineEdit_input1->text();
        QString outputname = m_MirtkUIOptions.lineEdit_output->text();

        if(inputDof.isEmpty()){
            QMessageBox::information(NULL, "Attention", "Images not selected correctly");
            return;
        }

        if(outputname.isEmpty()){
            QFileInfo fi(inputDof);
            outputname = "inverse_" + fi.baseName() + ".dof";
            MITK_INFO << ("No output name specified, using : " + outputname).toStdString();
        }
        outputname += (!outputname.contains(".dof")) ? ".dof" : "";

        QDir home(directory);
        std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
        cmd->SetUseDockerContainersOn();
        cmd->SetDockerImage("biomedia/mirtk:v1.1.0");
        QStringList arguments = cmd->GetDockerArguments(directory, "invert-dof");
        arguments << home.relativeFilePath(inputDof);
        arguments << home.relativeFilePath(outputname);

        QString outputFilePath = directory + "/" + outputname;

        bool successful = cmd->ExecuteCommand("docker", arguments, outputFilePath);
        MITK_INFO(successful) << "File created successfully";

    }
}

void QmitkCemrgAppCommonTools::MirtkOptsBrowse(const QString& buttDir){
    MITK_INFO << buttDir.toStdString();

    QStringList list = buttDir.split(",");
    QString buttonIdx = list.at(0);
    QString operationIdx = list.at(1);
    QString dir = list.at(2);
    QString titlelabel, input1, input2 = "";

    std::string msg1, msg2 = "";
    switch (operationIdx.toInt()) {
        case 0: // Registration
            titlelabel = "Select inputs to register images";
            msg1 = "Open moving (source) image for registration";
            msg2 = "Open fixed (target) image for registration";
            break;
        case 1: // Transformation
            titlelabel = "Select inputs to transform image or point set";
            msg1 = "Open image or point set for transformation";
            msg2 = "Open DOF file for transformation";
            break;
        case 2: // Inverse registration
            titlelabel = "Select inputs to inverse registration file";
            msg1 = "Open DOF file to calculate inverse";
            msg2 = "";
            break;
    }

    switch(buttonIdx.toInt()){
        case 1:
            QMessageBox::information(NULL, "Attention", msg1.c_str());
            input1 = QFileDialog::getOpenFileName(NULL, msg1.c_str(), dir, QmitkIOUtil::GetFileOpenFilterString());
            m_MirtkUIOptions.lineEdit_input1->setText(input1);
            break;
        case 2:
            QMessageBox::information(NULL, "Attention", msg2.c_str());
            input2 = QFileDialog::getOpenFileName(NULL, msg2.c_str(), dir, QmitkIOUtil::GetFileOpenFilterString());
            m_MirtkUIOptions.lineEdit_input2->setText(input2);
            break;
    }
    m_MirtkUIOptions.titleLabel->setText(titlelabel);
}
