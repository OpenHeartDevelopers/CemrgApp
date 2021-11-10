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
#include "QmitkCemrgAppCommonTools.h"

// Qt
#include <QMessageBox>
#include <QFileDialog>
#include <QSignalMapper>

const std::string QmitkCemrgAppCommonTools::VIEW_ID = "org.mitk.views.cemrgappcommontools";

void QmitkCemrgAppCommonTools::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.btn_loadmesh, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::LoadMesh);
    connect(m_Controls.btn_convert2carto, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::ConvertToCarto);
    connect(m_Controls.btn_vtk2cart, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::ConvertCarpToVtk);
    connect(m_Controls.btn_roiboxes, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::RoiControls);
    connect(m_Controls.btn_roiboxes_pts, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::RoiControlsSelectPoints);
    connect(m_Controls.btn_roiboxes_create, &QPushButton::clicked, this, &QmitkCemrgAppCommonTools::RoiControlsExtract);

    m_Controls.btn_roiboxes_pts->setVisible(false);
    m_Controls.btn_roiboxes_create->setVisible(false);
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

void QmitkCemrgAppCommonTools::RoiControls(){
    // loads image, opens pointsetinteraction, activates other buttons.
    QString path2image = QFileDialog::getOpenFileName(NULL, "Open image");
    if (path2image.isEmpty()){
        return;
    }

    mitk::IOUtil::Load(path2image.toStdString(), *this->GetDataStorage());
    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->GetDataStorage());

    int reply_found_points = QMessageBox::question(NULL, "Question", "Do you have points previously saved?");

    if (reply_found_points==QMessageBox::No){
        QMessageBox::information(NULL, "Attention", "Select positions of ROI using the PointSet Interaction Tool");
        this->GetSite()->GetPage()->ShowView("org.mitk.views.pointsetinteraction");
    } else if(reply_found_points==QMessageBox::Yes){
        QString path2roiboxes = QFileDialog::getOpenFileName(NULL, "Open file with box parameters");
        QFileInfo qfi(path2roiboxes);

        QString path = qfi.absolutePath();

        ifstream fi(path2roiboxes.toStdString());
        int numPoints;
        double xL, yL, zL, xC, yC, zC;

        fi >> numPoints;
        fi >> xL >> yL >> zL;
        std::cout << xL <<  ", " << yL << ", " << zL << '\n';
        std::cout << "Points: " << '\n';
        for (int ix = 0; ix < numPoints; ix++) {
            fi >> xC >> yC >> zC;
            std::cout << xC <<  ", " << yC << ", " << zC << '\n';

            std::vector<double> centre = {xC, yC, zC};
            std::vector<double> sides = {xL, yL, zL};
            mitk::Surface::Pointer cube = CemrgCommonUtils::CreateCube(centre, sides);

            mitk::DataStorage::SetOfObjects::ConstPointer sob = this->GetDataStorage()->GetAll();

            std::string cubename = ("cube_"+QString::number(ix)).toStdString();
            CemrgCommonUtils::AddToStorage(cube, cubename, this->GetDataStorage(), false);

            sob = this->GetDataStorage()->GetAll();
            for (mitk::DataStorage::SetOfObjects::ConstIterator nodeIt = sob->Begin(); nodeIt != sob->End(); ++nodeIt) {
                if (nodeIt->Value()->GetName().find("cube_") != nodeIt->Value()->GetName().npos) {
                    nodeIt->Value()->SetProperty("opacity", mitk::FloatProperty::New(0.4));
                    nodeIt->Value()->SetProperty("color", mitk::ColorProperty::New(0.0, 0.5, 0.5));
                }//_if
            }//_for
            m_Controls.btn_roiboxes_pts->setEnabled(false);
        }
    } else{
        return;
    }

    if(!m_Controls.btn_roiboxes_pts->isVisible()){
        m_Controls.btn_roiboxes_pts->setVisible(true);
        m_Controls.btn_roiboxes_create->setVisible(true);
    } else{
        m_Controls.btn_roiboxes_pts->setVisible(false);
        m_Controls.btn_roiboxes_create->setVisible(false);
    }
}

void QmitkCemrgAppCommonTools::RoiControlsSelectPoints(){
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) {
        QMessageBox::warning(NULL, "Attention", "Please select the PointSet from the Data Manager!");
        return;
    }//_if

    //Check selection type
    mitk::DataNode::Pointer landMarks = nodes.front();
    mitk::PointSet::Pointer pointSet = dynamic_cast<mitk::PointSet*>(landMarks->GetData());
    if (!pointSet) {
        QMessageBox::warning(NULL, "Attention", "PointSet error. Try selecting the points again!");
        return;
    }//_if

    double xL=1.0, yL=1.0, zL=1.0, xC=1.0, yC=1.0, zC=1.0;
    QString roiBoxesFilename = "prodRoiParameters.txt";

    if(m_Controls.btn_roiboxes_pts->text() == QString::fromStdString("Display ROI Boxes")){

        //Ask for user input to set the parameters
        QDialog* inputs = new QDialog(0, 0);
        m_RoiControls.setupUi(inputs);
        connect(m_RoiControls.buttonBox, SIGNAL(accepted()), inputs, SLOT(accept()));
        connect(m_RoiControls.buttonBox, SIGNAL(rejected()), inputs, SLOT(reject()));

        int dialogCode = inputs->exec();

        if (dialogCode == QDialog::Accepted) {
            std::cout << "x spinbox: " << m_RoiControls.spinbox_x->value() << ", ";
            std::cout << "y spinbox: " << m_RoiControls.spinbox_y->value() << ", ";
            std::cout << "z spinbox: " << m_RoiControls.spinbox_z->value() << '\n';
            xL = m_RoiControls.spinbox_x->value();
            yL = m_RoiControls.spinbox_y->value();
            zL = m_RoiControls.radio_set2d->isChecked() ? 1.0 : m_RoiControls.spinbox_z->value();

            if(!m_RoiControls.lineedit_name->text().isEmpty()){
                roiBoxesFilename = m_RoiControls.lineedit_name->text();
                roiBoxesFilename += (!roiBoxesFilename.contains(".txt")) ? ".txt" : "";
            }

            inputs->deleteLater();
        } else if (dialogCode == QDialog::Rejected) {
            inputs->close();
            inputs->deleteLater();

            return;
        }//_if

        QString dir = QFileDialog::getExistingDirectory(NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(), QFileDialog::ShowDirsOnly | QFileDialog::DontUseNativeDialog);
        QString prodFile = dir + "/" + roiBoxesFilename;
        ofstream fo(prodFile.toStdString());
        fo << pointSet->GetSize() << '\n'; // number of points
        fo << xL << " " << yL << " " << zL << '\n'; // sizes

        for (int ix = 0; ix < pointSet->GetSize(); ix++) {
            xC = pointSet->GetPoint(ix).GetElement(0);
            yC = pointSet->GetPoint(ix).GetElement(1);
            zC = pointSet->GetPoint(ix).GetElement(2);

            fo << std::setprecision(6) << xC << " ";
            fo << std::setprecision(6) << yC << " ";
            fo << std::setprecision(6) << zC << '\n';

            std::vector<double> centre = {xC, yC, zC};
            std::vector<double> sides = {xL, yL, zL};
            mitk::Surface::Pointer cube = CemrgCommonUtils::CreateCube(centre, sides);

            mitk::DataStorage::SetOfObjects::ConstPointer sob = this->GetDataStorage()->GetAll();

            std::string cubename = ("cube_"+QString::number(ix)).toStdString();
            CemrgCommonUtils::AddToStorage(cube, cubename, this->GetDataStorage(), false);

            sob = this->GetDataStorage()->GetAll();
            for (mitk::DataStorage::SetOfObjects::ConstIterator nodeIt = sob->Begin(); nodeIt != sob->End(); ++nodeIt) {
                if (nodeIt->Value()->GetName().find("cube_") != nodeIt->Value()->GetName().npos) {
                    nodeIt->Value()->SetProperty("opacity", mitk::FloatProperty::New(0.4));
                    nodeIt->Value()->SetProperty("color", mitk::ColorProperty::New(0.0, 0.5, 0.5));
                }//_if
            }//_for
        }//_for each point selected by user

        m_Controls.btn_roiboxes_pts->setText("Clear ROI Boxes");
    } else{
        std::string title = "Question";
        std::string msg = "Clear current boxes?";
        int reply = QMessageBox::question(NULL, title.c_str(), msg.c_str());
        if(reply==QMessageBox::Yes){
            // remove everything
            mitk::DataStorage::SetOfObjects::ConstPointer sob = this->GetDataStorage()->GetAll();

            for (mitk::DataStorage::SetOfObjects::ConstIterator nodeIt = sob->Begin(); nodeIt != sob->End(); ++nodeIt){
                if (nodeIt->Value()->GetName().find("cube_") != nodeIt->Value()->GetName().npos){
                    this->GetDataStorage()->Remove(nodeIt->Value());
                }
            }
        }
        m_Controls.btn_roiboxes_pts->setText("Display ROI Boxes");
    }
}

void QmitkCemrgAppCommonTools::RoiControlsExtract(){
    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    if (nodes.empty()) {
        QMessageBox::warning(NULL, "Attention",
            "Please select an image from the Data Manager to perform cropping!");
        return;
    }//_if

    QString dir = QFileDialog::getExistingDirectory(NULL, "Open Project Directory", mitk::IOUtil::GetProgramPath().c_str(), QFileDialog::ShowDirsOnly | QFileDialog::DontUseNativeDialog);
    QString prodFile = dir + "/roi_";

    // get things form node
    mitk::DataNode::Pointer imageNode = nodes.at(0);
    mitk::Image::Pointer imageToCut;
    mitk::BaseData::Pointer data = imageNode->GetData();

    if(data){
        imageToCut = dynamic_cast<mitk::Image*>(data.GetPointer());
        if(imageToCut){

            CemrgCommonUtils::SetImageNode(imageNode);
            CemrgCommonUtils::SetImageToCut(imageToCut);

            int ix = 0;
            mitk::DataStorage::SetOfObjects::ConstPointer sob = this->GetDataStorage()->GetAll();
            for (mitk::DataStorage::SetOfObjects::ConstIterator nodeIt = sob->Begin(); nodeIt != sob->End(); ++nodeIt) {
                if (nodeIt->Value()->GetName().find("cube_") != nodeIt->Value()->GetName().npos) {
                    MITK_INFO << nodeIt->Value()->GetName();
                    mitk::BoundingObject::Pointer cuttingCube;
                    cuttingCube = dynamic_cast<mitk::BoundingObject*>(nodeIt->Value()->GetData());
                    cuttingCube->FitGeometry(imageToCut->GetGeometry());

                    CemrgCommonUtils::SetCuttingCube(cuttingCube);

                    mitk::Image::Pointer outputImage = CemrgCommonUtils::CropImage();
                    mitk::IOUtil::Save(outputImage, (prodFile+QString::number(ix)+".nii").toStdString());
                    ix++;
                }//_if
            }//_for
        } else return;
    } else return;
}
