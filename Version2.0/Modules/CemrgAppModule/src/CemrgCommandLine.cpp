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
 * Commandline Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

//Qmitk
#include <mitkIOUtil.h>
#include <mitkProgressBar.h>

//Qt
#include <QFileDialog>
#include <QFileInfo>
#include <QFile>
#include <QCoreApplication>
#include <QDebug>
#include <QDir>
#include <QMessageBox>

//Generic
#include <thread>
#include <chrono>
#include <sys/stat.h>
#include "CemrgCommandLine.h"

CemrgCommandLine::CemrgCommandLine() {

    isUI = true;
    _useDockerContainers = true;
    _dockerimage = "biomedia/mirtk:v1.1.0";

    //Setup panel
    panel = new QTextEdit(0,0);
    QPalette palette = panel->palette();
    palette.setColor(QPalette::Base, Qt::black);
    palette.setColor(QPalette::Text, Qt::red);
    panel->setPalette(palette);
    panel->setReadOnly(true);

    //Setup dialog
    layout = new QVBoxLayout();
    dial = new QDialog(0,0);
    dial->setFixedSize(640, 480);
    dial->setLayout(layout);
    dial->layout()->addWidget(panel);
    dial->show();

    //Setup the process
    process = std::unique_ptr<QProcess>(new QProcess(this));
    process->setProcessChannelMode(QProcess::MergedChannels);
    connect(process.get(), SIGNAL(readyReadStandardOutput()), this, SLOT(UpdateStdText()));
    connect(process.get(), SIGNAL(readyReadStandardError()), this, SLOT(UpdateErrText()));
    connect(process.get(), SIGNAL(finished(int)), this, SLOT(FinishedAlert()));
}

CemrgCommandLine::CemrgCommandLine(bool cmd) {

    isUI = cmd;
    _useDockerContainers = true;
    _dockerimage = "biomedia/mirtk:v1.1.0";

    if(cmd) {
        //Setup panel
        panel = new QTextEdit(0,0);
        QPalette palette = panel->palette();
        palette.setColor(QPalette::Base, Qt::black);
        palette.setColor(QPalette::Text, Qt::red);
        panel->setPalette(palette);
        panel->setReadOnly(true);

        //Setup dialog
        layout = new QVBoxLayout();
        dial = new QDialog(0,0);
        dial->setFixedSize(640, 480);
        dial->setLayout(layout);
        dial->layout()->addWidget(panel);
        dial->show();
    }

    //Setup the process
    process = std::unique_ptr<QProcess>(new QProcess(this));
    process->setProcessChannelMode(QProcess::MergedChannels);
    connect(process.get(), SIGNAL(readyReadStandardOutput()), this, SLOT(UpdateStdText()));
    connect(process.get(), SIGNAL(readyReadStandardError()), this, SLOT(UpdateErrText()));
    connect(process.get(), SIGNAL(finished(int)), this, SLOT(FinishedAlert()));
}

CemrgCommandLine::CemrgCommandLine(std::string dockerimage) {

    isUI = true;
    _useDockerContainers = true;
    _dockerimage = QString::fromStdString(dockerimage);

    //Setup panel
    panel = new QTextEdit(0,0);
    QPalette palette = panel->palette();
    palette.setColor(QPalette::Base, Qt::black);
    palette.setColor(QPalette::Text, Qt::red);
    panel->setPalette(palette);
    panel->setReadOnly(true);

    //Setup dialog
    layout = new QVBoxLayout();
    dial = new QDialog(0,0);
    dial->setFixedSize(640, 480);
    dial->setLayout(layout);
    dial->layout()->addWidget(panel);
    dial->show();

    //Setup the process
    process = std::unique_ptr<QProcess>(new QProcess(this));
    process->setProcessChannelMode(QProcess::MergedChannels);
    connect(process.get(), SIGNAL(readyReadStandardOutput()), this, SLOT(UpdateStdText()));
    connect(process.get(), SIGNAL(readyReadStandardError()), this, SLOT(UpdateErrText()));
    connect(process.get(), SIGNAL(finished(int)), this, SLOT(FinishedAlert()));
}

CemrgCommandLine::~CemrgCommandLine() {

    process->close();
    dial->deleteLater();
    panel->deleteLater();
    layout->deleteLater();
}

/***************************************************************************
 ************************* BUILDING MESH UTILITIES *************************
 ***************************************************************************/

QString CemrgCommandLine::ExecuteSurf(QString dir, QString segPath, int iter, float th, int blur, int smth) {
     MITK_INFO << "[ATTENTION] SURFACE CREATION: Close -> Surface -> Smooth";
     QString closeOutputPath, surfOutputPath;
     QString outAbsolutepath = "ERROR_IN_PROCESSING";

     closeOutputPath = ExecuteMorphologicalOperation("close", dir, segPath, "segmentation.s.nii", iter);
     if (QString::compare(closeOutputPath, "ERROR_IN_PROCESSING")!=0){
         surfOutputPath = ExecuteExtractSurface(dir, closeOutputPath, "segmentation.vtk", th, blur);
         if (QString::compare(surfOutputPath, "ERROR_IN_PROCESSING")!=0){
             outAbsolutepath = ExecuteSmoothSurface(dir, surfOutputPath, surfOutputPath, smth);
             remove((dir + mitk::IOUtil::GetDirectorySeparator() + "segmentation.s.nii").toStdString().c_str());
         }
     }

     return outAbsolutepath;
 }

 QString CemrgCommandLine::ExecuteCreateCGALMesh(QString dir, QString outputName, QString paramsFullPath, QString segmentationName) {
     MITK_INFO << "[ATTENTION] Attempting meshtools3d libraries.";

     QString revertDockerImage = "";
     if (!_dockerimage.contains("meshtools3d", Qt::CaseInsensitive)){
         MITK_INFO << "Changing docker image name to meshtools3d.";
         revertDockerImage = getDockerImage(); //get current docker image
         setDockerImage(QString("alonsojasl/meshtools3d:v1.0"));
     }

     QString executablePath = "";
     QString executableName;
     QDir meshtools3dhome(dir);

     QString outAbsolutepath, outputDirectory, segmentationDirectory;
     QStringList arguments;

     segmentationDirectory = dir + mitk::IOUtil::GetDirectorySeparator();
     outputDirectory = segmentationDirectory + "CGALMeshDir";
     outAbsolutepath = outputDirectory + mitk::IOUtil::GetDirectorySeparator() + outputName;

     if(_useDockerContainers){
         MITK_INFO << "Using docker containers.";
         #if defined(__APPLE__)
         executablePath = "/usr/local/bin/";
         #endif

         executableName = executablePath+"docker";
         arguments = getDockerArguments(meshtools3dhome.absolutePath());

         arguments << "-f" << meshtools3dhome.relativeFilePath(paramsFullPath);
         arguments << "-seg_dir" << meshtools3dhome.relativeFilePath(segmentationDirectory);;
         arguments << "-seg_name" << segmentationName;
         arguments << "-out_dir" << meshtools3dhome.relativeFilePath(outputDirectory);
         arguments << "-out_name" << outputName;

     } else{
         MITK_INFO << "Using static MIRTK libraries.";
         #if defined(__APPLE__)
         executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
             mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
             mitk::IOUtil::GetDirectorySeparator() + QString("M3DLib");
         #endif

         executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "M3DLib";
         executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + "meshtools3d";
         QDir apathd(executablePath);
         if (apathd.exists()) {
             process->setWorkingDirectory(executablePath);

             arguments << "-f" << paramsFullPath;
             arguments << "-seg_dir" << segmentationDirectory;;
             arguments << "-seg_name" << segmentationName;
             arguments << "-out_dir" << outputDirectory;
             arguments << "-out_name" << outputName;

         } else {
             QMessageBox::warning(NULL, "Please check the LOG", "MESHTOOLS3D libraries not found");
             MITK_WARN << "MESHTOOLS3D libraries not found. Please make sure the M3DLib folder is inside the directory:\n\t"+
             mitk::IOUtil::GetProgramPath();
         }//_if
     }

     //Setup EnVariable - in windows TBB_NUM_THREADS should be set in the system environment variables
 #ifndef _WIN32
     QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
     env.insert("TBB_NUM_THREADS","12");
     process->setProcessEnvironment(env);
 #endif
     bool successful = ExecuteCommand(executableName, arguments, outAbsolutepath);

     if (!revertDockerImage.isEmpty()){ // revert to original docker image (in case the object is used later).
         setDockerImage(revertDockerImage);
     }

     if(!successful) {
         if (_useDockerContainers){
             MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
             setUseDockerContainersOff();
             return ExecuteCreateCGALMesh(dir, outputName, paramsFullPath, segmentationName);
         } else{
             MITK_WARN << "Local MESHTOOLS3D libraries did not produce a good outcome.";
             return "ERROR_IN_PROCESSING";
         }
     } else {
         return outAbsolutepath;
     }
 }

/***************************************************************************
 **************************** TRACKING UTILITIES ***************************
 ***************************************************************************/

 void CemrgCommandLine::ExecuteTracking(QString dir, QString imgTimes, QString param, QString output) {
     MITK_INFO << "[ATTENTION] Attempting Registration.";

     QString executablePath = "";
     QString executableName;
     QString commandName = "register";
     QDir mirtkhome(dir);

     QString imgTimesFilePath, outAbsolutepath;
     QString prodPath = dir + mitk::IOUtil::GetDirectorySeparator();
     QStringList arguments;

     imgTimesFilePath = imgTimes.contains(dir, Qt::CaseSensitive) ? imgTimes : prodPath + imgTimes;
     outAbsolutepath = output.contains(dir, Qt::CaseSensitive) ? output : prodPath + output;

     if(!outAbsolutepath.contains(".dof", Qt::CaseSensitive)) outAbsolutepath += ".dof";

     MITK_INFO << ("[...] IMAGE FILES PATH: " + imgTimesFilePath).toStdString();
//     MITK_INFO << ("[...] FIXED (target): " + fixedfullpath).toStdString();
     MITK_INFO << ("[...] OUTPUT DOF: " + outAbsolutepath).toStdString();

     if(_useDockerContainers){
         MITK_INFO << "Using docker containers.";
         #if defined(__APPLE__)
         executablePath = "/usr/local/bin/";
         #endif

         executableName = executablePath+"docker";
         arguments = getDockerArguments(mirtkhome.absolutePath(), commandName);

         arguments << "-images" << mirtkhome.relativeFilePath(imgTimesFilePath);
         if (!param.isEmpty()) arguments << "-parin" << mirtkhome.relativeFilePath(param);
         arguments << "-dofout" << mirtkhome.relativeFilePath(outAbsolutepath);
         arguments << "-threads" << "12";
         arguments << "-verbose" << "3";

     } else{
         MITK_INFO << "Using static MIRTK libraries.";
         #if defined(__APPLE__)
         executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
             mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
             mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
         #endif

         executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
         executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + commandName;
         QDir apathd(executablePath);
         if (apathd.exists()) {
             process->setWorkingDirectory(executablePath);

             arguments << "-images" << imgTimesFilePath;
             if (!param.isEmpty()) arguments << "-parin" << param;
             arguments << "-dofout" << outAbsolutepath;
             arguments << "-threads" << "12";
             arguments << "-verbose" << "3";

         } else {
             QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
             MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
             mitk::IOUtil::GetProgramPath();
         }//_if
     }

     bool successful = ExecuteCommand(executableName, arguments, outAbsolutepath);

     if(!successful) {
         if (_useDockerContainers){
             MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
             setUseDockerContainersOff();
             ExecuteTracking(dir, imgTimes, param, output);
         } else{
             MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
         }
     }
 }

 void CemrgCommandLine::ExecuteApplying(QString dir, QString inputMesh, double iniTime, QString dofin, int noFrames, int smooth) {
     int fctTime = 10;
     noFrames *= smooth;
     if (smooth == 2){
         fctTime = 5;
     } else if (smooth == 5){
         fctTime = 2;
     }

     QString output = dir + mitk::IOUtil::GetDirectorySeparator() + "transformed-";
     QString thisOutput;
     for (int i=0; i<noFrames; i++) {
         thisOutput = output+QString::number(i)+".vtk";
         ExecuteTransformationOnPoints(dir, inputMesh, thisOutput, dofin, iniTime);
     }
 }

void CemrgCommandLine::ExecuteRegistration(QString dir, QString fixed, QString moving, QString txname, QString modelname) {
    MITK_INFO << "[ATTENTION] Attempting Registration.";
    //lge : fixed   ||   mra : moving
    QString executablePath = "";
    QString executableName;
    QString commandName = "register";
    QDir mirtkhome(dir);

    QString fixedfullpath, movingfullpath, outAbsolutepath;
    QString prodPath = dir + mitk::IOUtil::GetDirectorySeparator();
    QStringList arguments;

    fixedfullpath = fixed.contains(dir, Qt::CaseSensitive) ? fixed : prodPath + fixed;
    movingfullpath = moving.contains(dir, Qt::CaseSensitive) ? moving : prodPath + moving;
    outAbsolutepath = txname.contains(dir, Qt::CaseSensitive) ? txname : prodPath + txname;

    if(!fixedfullpath.contains(".nii", Qt::CaseSensitive)) fixedfullpath += ".nii";
    if(!movingfullpath.contains(".nii", Qt::CaseSensitive)) movingfullpath += ".nii";
    if(!outAbsolutepath.contains(".dof", Qt::CaseSensitive)) movingfullpath += ".dof";

    MITK_INFO << ("[...] MOVING (source): " + movingfullpath).toStdString();
    MITK_INFO << ("[...] FIXED (target): " + fixedfullpath).toStdString();
    MITK_INFO << ("[...] OUTPUT DOF: " + outAbsolutepath).toStdString();

    if(_useDockerContainers){
        MITK_INFO << "Using docker containers.";
        #if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
        #endif

        executableName = executablePath+"docker";
        arguments = getDockerArguments(mirtkhome.absolutePath(), commandName);

        arguments << mirtkhome.relativeFilePath(movingfullpath); // input1
        arguments << mirtkhome.relativeFilePath(fixedfullpath); // input2
        arguments << "-dofout" << mirtkhome.relativeFilePath(outAbsolutepath);
        arguments << "-model" << modelname;
        arguments << "-verbose" << "3";
        arguments << "-color";

    } else{
        MITK_INFO << "Using static MIRTK libraries.";
        #if defined(__APPLE__)
        executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
            mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
            mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
        #endif

        executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
        executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + commandName;
        QDir apathd(executablePath);
        if (apathd.exists()) {
            process->setWorkingDirectory(executablePath);

            arguments << movingfullpath;
            arguments << fixedfullpath;
            arguments << "-dofout" << outAbsolutepath;
            arguments << "-model" << modelname;
            arguments << "-verbose" << "3";

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
            mitk::IOUtil::GetProgramPath();
        }//_if
    }

    MITK_INFO << ("Performing a " + modelname + " registration").toStdString();
    bool successful = ExecuteCommand(executableName, arguments, outAbsolutepath);

    if(!successful) {
        if (_useDockerContainers){
            MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
            setUseDockerContainersOff();
            ExecuteRegistration(dir, fixed, moving, txname, modelname);
        } else{
            MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
        }
    }
}

void CemrgCommandLine::ExecuteTransformation(QString dir, QString imgname, QString regname, QString txfullpath) {
    MITK_INFO << "[ATTENTION] Attempting Image Transformation.";

    QString executablePath = "";
    QString executableName;
    QString commandName = "transform-image";
    QDir mirtkhome(dir);

    QString dofpath, imgNamefullpath, outAbsolutepath;
    QString prodPath = dir + mitk::IOUtil::GetDirectorySeparator();
    QStringList arguments;

    imgNamefullpath = imgname.contains(dir, Qt::CaseSensitive) ? imgname : prodPath + imgname;
    outAbsolutepath = regname.contains(dir, Qt::CaseSensitive) ? regname : prodPath + regname;
    dofpath = txfullpath.contains(dir, Qt::CaseSensitive) ? txfullpath : prodPath + txfullpath;

    MITK_INFO << ("[...] INPUT IMAGE: " + imgNamefullpath).toStdString();
    MITK_INFO << ("[...] INPUT DOF: " + dofpath).toStdString();
    MITK_INFO << ("[...] OUTPUT IMAGE: " + outAbsolutepath).toStdString();

    if(_useDockerContainers){
        MITK_INFO << "Using docker containers.";
        #if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
        #endif

        executableName = executablePath+"docker";
        arguments = getDockerArguments(mirtkhome.absolutePath(), commandName);

        arguments << mirtkhome.relativeFilePath(imgNamefullpath); //input
        arguments << mirtkhome.relativeFilePath(outAbsolutepath); //output
        arguments << "-dof" << mirtkhome.relativeFilePath(dofpath);
        arguments << "-verbose" << "3";

    } else{
        MITK_INFO << "Using static MIRTK libraries.";
        #if defined(__APPLE__)
        executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
            mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
            mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
        #endif

        executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
        executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + commandName;
        QDir apathd(executablePath);
        if (apathd.exists()) {
            process->setWorkingDirectory(executablePath);

            arguments << imgNamefullpath; //input
            arguments << outAbsolutepath; //output
            arguments << "-dof" << dofpath;
            arguments << "-verbose" << "3";

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
            mitk::IOUtil::GetProgramPath();
        }//_if
    }

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutepath);

    if(!successful) {
        if (_useDockerContainers){
            MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
            setUseDockerContainersOff();
            ExecuteTransformation(dir, imgname, regname, txfullpath);
        } else{
            MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
        }
    }
}

void CemrgCommandLine::ExecuteResamplingOnNifti(
        QString niifullpath, QString outputtniifullpath, int isovalue) {
    MITK_INFO << "[ATTENTION] Attempting Image Transformation.";

    QString executablePath = "";
    QString executableName;
    QString commandName = "resample-image";
    QFileInfo inputnii(niifullpath);
    QString dir = inputnii.absolutePath();
    QDir mirtkhome(dir);

    QString dofpath, imgNamefullpath, outAbsolutepath;
    QString prodPath = dir + mitk::IOUtil::GetDirectorySeparator();
    QStringList arguments;

    imgNamefullpath = niifullpath;
    outAbsolutepath = outputtniifullpath;

    MITK_INFO << ("[...] INPUT IMAGE: " + imgNamefullpath).toStdString();
    MITK_INFO << ("[...] OUTPUT IMAGE: " + outAbsolutepath).toStdString();

    if(_useDockerContainers){
        MITK_INFO << "Using docker containers.";
        #if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
        #endif

        executableName = executablePath+"docker";
        arguments = getDockerArguments(mirtkhome.absolutePath(), commandName);

        arguments << mirtkhome.relativeFilePath(imgNamefullpath); //input
        arguments << mirtkhome.relativeFilePath(outAbsolutepath); //output
        arguments << "-isotropic" << QString::number(isovalue);
        arguments << "-interp" << "CSpline";
        arguments << "-verbose" << "3";

    } else{
        MITK_INFO << "Using static MIRTK libraries.";
        #if defined(__APPLE__)
        executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") + mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") + mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
        #endif

        executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
        executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + commandName;
        QDir apathd(executablePath);
        if (apathd.exists()) {
            process->setWorkingDirectory(executablePath);

            arguments << imgNamefullpath; //input
            arguments << outAbsolutepath; //output
            arguments << "-isotropic" << QString::number(isovalue);
            arguments << "-interp" << "CSpline";
            arguments << "-verbose" << "3";

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
            mitk::IOUtil::GetProgramPath();
        }//_if
    }

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutepath);

    if(!successful) {
        if (_useDockerContainers){
            MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
            setUseDockerContainersOff();
            ExecuteResamplingOnNifti(niifullpath, outputtniifullpath,isovalue);
        } else{
            MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
        }
    }
}

void CemrgCommandLine::ExecuteTransformationOnPoints(QString dir, QString meshfullpath, QString outputtmeshfullpath, QString txfullpath, double applyingIniTime) {
    MITK_INFO << "[ATTENTION] Attempting Pointset transformation.";
    QString executablePath = "";
    QString executableName;
    QString commandName = "transform-points";
    QDir mirtkhome(dir);

    QString dofpath, inputMeshFullPath, outAbsolutepath;
    QString prodPath = dir + mitk::IOUtil::GetDirectorySeparator();
    QStringList arguments;

    inputMeshFullPath = meshfullpath.contains(dir, Qt::CaseSensitive) ? meshfullpath : prodPath + meshfullpath;
    outAbsolutepath = outputtmeshfullpath.contains(dir, Qt::CaseSensitive) ? outputtmeshfullpath : prodPath + outputtmeshfullpath;
    dofpath = txfullpath.contains(dir, Qt::CaseSensitive) ? txfullpath : prodPath + txfullpath;

    MITK_INFO << ("[...] INPUT IMAGE: " + inputMeshFullPath).toStdString();
    MITK_INFO << ("[...] INPUT DOF: " + dofpath).toStdString();
    MITK_INFO << ("[...] OUTPUT IMAGE: " + outAbsolutepath).toStdString();

    if(_useDockerContainers){
        MITK_INFO << "Using docker containers.";
        #if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
        #endif

        executableName = executablePath+"docker";
        arguments = getDockerArguments(mirtkhome.absolutePath(), commandName);

        arguments << mirtkhome.relativeFilePath(inputMeshFullPath); // input
        arguments << mirtkhome.relativeFilePath(outAbsolutepath); // output
        arguments << "-dofin" << mirtkhome.relativeFilePath(dofpath);
        arguments << "-ascii";
        if (applyingIniTime != -100){
            // -100 is the default value indicating ExecuteApplying is not being called.
            arguments << "-St";
            arguments << QString::number(applyingIniTime);
        }
        arguments << "-verbose" << "3";

    } else{
        MITK_INFO << "Using static MIRTK libraries.";
        #if defined(__APPLE__)
        executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
            mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
            mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
        #endif

        executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
        executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + commandName;
        QDir apathd(executablePath);
        if (apathd.exists()) {
            process->setWorkingDirectory(executablePath);

            arguments << inputMeshFullPath; // input
            arguments << outAbsolutepath; // output
            arguments << "-dofin" << dofpath;
            arguments << "-ascii";
            if (applyingIniTime != -100){
                // -100 is the default value indicating ExecuteApplying is not being called.
                arguments << "-St";
                arguments << QString::number(applyingIniTime);
            }
            arguments << "-verbose" << "3";

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
            mitk::IOUtil::GetProgramPath();
        }//_if
    }

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutepath);

    if(!successful) {
        if (_useDockerContainers){
            MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
            setUseDockerContainersOff();
            ExecuteTransformationOnPoints(dir, meshfullpath, outputtmeshfullpath, txfullpath);
        } else{
            MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
        }
    }
}

QString CemrgCommandLine::ExecuteExpandSurf(QString dir, QString segPath, int iter, float th, int blur, int smth){
    MITK_INFO << "[ATTENTION] EXPANDED SURFACE CREATION: Dilate -> Surface -> Smooth";
    QString dilateOutputPath, surfOutputPath;
    QString outAbsolutepath = "ERROR_IN_PROCESSING";

    dilateOutputPath = ExecuteMorphologicalOperation("dilate", dir, segPath, "segmentation.s.nii", iter);
    if (QString::compare(dilateOutputPath, "ERROR_IN_PROCESSING")!=0){
        surfOutputPath = ExecuteExtractSurface(dir, dilateOutputPath, "segmentation.vtk", th, blur);
        if (QString::compare(surfOutputPath, "ERROR_IN_PROCESSING")!=0){
            outAbsolutepath = ExecuteSmoothSurface(dir, surfOutputPath, surfOutputPath, smth);
            remove((dir + mitk::IOUtil::GetDirectorySeparator() + "segmentation.s.nii").toStdString().c_str());
        }
    }

    return outAbsolutepath;
}

void CemrgCommandLine::ExecuteSimpleTranslation(QString dir, QString sourceMeshP, QString targetMeshP, QString txName, bool transformThePoints){
    MITK_INFO << "[ATTENTION] Attempting INIT-DOF.";
    QString aPath, executableName, commandName, sourceMeshPath, targetMeshPath, outAbsolutepath, prodPath;
    QStringList arguments;
    QDir mirtkhome(dir);

    aPath = "";
    commandName = "init-dof"; //simple translation

    sourceMeshPath = sourceMeshP.contains(dir, Qt::CaseSensitive) ? sourceMeshP : prodPath + sourceMeshP;
    targetMeshPath = targetMeshP.contains(dir, Qt::CaseSensitive) ? targetMeshP : prodPath + targetMeshP;
    outAbsolutepath = txName.contains(dir, Qt::CaseSensitive) ? txName : prodPath + txName;

    if(_useDockerContainers){
        MITK_INFO << "Using docker containers.";
        #if defined(__APPLE__)
        aPath = "/usr/local/bin/";
        #endif

        executableName = aPath+"docker";
        arguments = getDockerArguments(mirtkhome.absolutePath(),  commandName);

        arguments << mirtkhome.relativeFilePath(outAbsolutepath);
        arguments << "-translations" << "-norotations" << "-noscaling" << "-noshearing";
        if(transformThePoints) {
            arguments << "-displacements";
            arguments << mirtkhome.relativeFilePath(sourceMeshPath);
            arguments << mirtkhome.relativeFilePath(targetMeshPath);
        }
        else{
            arguments << "-source" << mirtkhome.relativeFilePath(sourceMeshPath);
            arguments << "-target" << mirtkhome.relativeFilePath(targetMeshPath);
        }
        arguments << "-verbose" << "3";
    }else{
        MITK_INFO << "Using static MIRTK libraries.";
        aPath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
        executableName = aPath + mitk::IOUtil::GetDirectorySeparator() + commandName;
        #if defined(__APPLE__)
        aPath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
        #endif
        QDir apathd(aPath);
        if (apathd.exists()) {
            process->setWorkingDirectory(aPath);

            arguments << outAbsolutepath;
            arguments << "-translations" << "-norotations" << "-noscaling" << "-noshearing";
            if(transformThePoints) {
                arguments << "-displacements";
                arguments << sourceMeshPath;
                arguments << targetMeshPath;
            }
            else{
                arguments << "-source" << sourceMeshPath;
                arguments << "-target" << targetMeshPath;
            }
            arguments << "-verbose" << "3";

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
                         mitk::IOUtil::GetProgramPath();
        }//_if
    }

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutepath);

    if(!successful) {
        if (_useDockerContainers){
            MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
            setUseDockerContainersOff();
            ExecuteSimpleTranslation(dir, sourceMeshP, targetMeshP, txName, transformThePoints);
        } else{
            MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
        }
    }
}

/***************************************************************************
 ***************************************************************************/
QString CemrgCommandLine::ExecuteMorphologicalOperation(QString operation, QString dir, QString segPath, QString outputPath, int iter){
    MITK_INFO << "[ATTENTION] Attempting Pointset transformation.";
    QString executablePath = "";
    QString executableName;
    QString commandName;
    QStringList arguments;

    if (QString::compare(operation, "dilate", Qt::CaseInsensitive)==0){
        commandName = "dilate-image";
    } else if(QString::compare(operation, "erode", Qt::CaseInsensitive)==0){
        commandName = "erode-image";
    } else if(QString::compare(operation, "open", Qt::CaseInsensitive)==0){
        commandName = "open-image";
    } else if(QString::compare(operation, "close", Qt::CaseInsensitive)==0){
        commandName = "close-image";
    } else{
        MITK_ERROR << ("Morphological operation: " + operation + " misspelled or not supported.").toStdString();
        return "ERROR_IN_PROCESSING";
    }

    QDir mirtkhome(dir);

    QString inputImgFullPath, outAbsolutepath;
    QString prodPath = dir + mitk::IOUtil::GetDirectorySeparator();

    inputImgFullPath = segPath.contains(dir, Qt::CaseSensitive) ? segPath : prodPath + segPath;
    outAbsolutepath = outputPath.contains(dir, Qt::CaseSensitive) ? outputPath : prodPath + outputPath;

    MITK_INFO << ("[...] OPERATION: " + operation).toStdString();
    MITK_INFO << ("[...] INPUT IMAGE: " + inputImgFullPath).toStdString();
    MITK_INFO << ("[...] OUTPUT IMAGE: " + outAbsolutepath).toStdString();

    if(_useDockerContainers){
        MITK_INFO << "Using docker containers.";
        #if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
        #endif

        executableName = executablePath+"docker";
        arguments = getDockerArguments(mirtkhome.absolutePath(), commandName);

        arguments << mirtkhome.relativeFilePath(inputImgFullPath);
        arguments << mirtkhome.relativeFilePath(outAbsolutepath);
        arguments << "-iterations" << QString::number(iter);

    } else{
        MITK_INFO << "Using static MIRTK libraries.";
        #if defined(__APPLE__)
        executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
            mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
            mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
        #endif

        executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
        executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + commandName;
        QDir apathd(executablePath);
        if (apathd.exists()) {
            process->setWorkingDirectory(executablePath);

            arguments << inputImgFullPath;
            arguments << outAbsolutepath;
            arguments << "-iterations" << QString::number(iter);

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
            mitk::IOUtil::GetProgramPath();
        }//_if
    }

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutepath);

    if(!successful) {
        if (_useDockerContainers){
            MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
            setUseDockerContainersOff();
            return ExecuteMorphologicalOperation(operation, dir, segPath, outputPath, iter);
        } else{
            MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
            return "ERROR_IN_PROCESSING";
        }
    } else {
        return outAbsolutepath;
    }
}

QString CemrgCommandLine::ExecuteExtractSurface(QString dir, QString segPath, QString outputPath,float th, int blur){
    MITK_INFO << "[ATTENTION] Attempting Surface extraction.";
    QString executablePath = "";
    QString executableName;
    QString commandName;

    commandName = "extract-surface";

    QDir mirtkhome(dir);

    QString inputImgFullPath, outAbsolutepath;
    QString prodPath = dir + mitk::IOUtil::GetDirectorySeparator();
    QStringList arguments;

    inputImgFullPath = segPath.contains(dir, Qt::CaseSensitive) ? segPath : prodPath + segPath;
    outAbsolutepath = outputPath.contains(dir, Qt::CaseSensitive) ? outputPath : prodPath + outputPath;

    MITK_INFO << ("[...] INPUT IMAGE: " + inputImgFullPath).toStdString();
    MITK_INFO << ("[...] OUTPUT MESH: " + outAbsolutepath).toStdString();

    if(_useDockerContainers){
        MITK_INFO << "Using docker containers.";
        #if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
        #endif

        executableName = executablePath+"docker";
        arguments = getDockerArguments(mirtkhome.absolutePath(), commandName);

        arguments << mirtkhome.relativeFilePath(inputImgFullPath);
        arguments << mirtkhome.relativeFilePath(outAbsolutepath);
        arguments << "-isovalue" << QString::number(th);
        arguments << "-blur" << QString::number(blur);
        arguments << "-ascii";
        arguments << "-verbose" << "3";

    } else{
        MITK_INFO << "Using static MIRTK libraries.";
        #if defined(__APPLE__)
        executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
            mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
            mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
        #endif

        executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
        executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + commandName;
        QDir apathd(executablePath);
        if (apathd.exists()) {
            process->setWorkingDirectory(executablePath);

            arguments << inputImgFullPath;
            arguments << outAbsolutepath;
            arguments << "-isovalue" << QString::number(th);
            arguments << "-blur" << QString::number(blur);
            arguments << "-ascii";
            arguments << "-verbose" << "3";

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
            mitk::IOUtil::GetProgramPath();
        }//_if
    }

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutepath);

    if(!successful) {
        if (_useDockerContainers){
            MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
            setUseDockerContainersOff();
            return ExecuteExtractSurface(dir, segPath, outputPath,th, blur);
        } else{
            MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
            return "ERROR_IN_PROCESSING";
        }
    } else{
        return outAbsolutepath;
    }
}

QString CemrgCommandLine::ExecuteSmoothSurface(QString dir, QString segPath, QString outputPath, int smth){
    MITK_INFO << "[ATTENTION] Attempting Surface extraction.";
    QString executablePath = "";
    QString executableName;
    QString commandName;

    commandName = "smooth-surface";

    QDir mirtkhome(dir);

    QString inputMeshFullPath, outAbsolutepath;
    QString prodPath = dir + mitk::IOUtil::GetDirectorySeparator();
    QStringList arguments;

    inputMeshFullPath = segPath.contains(dir, Qt::CaseSensitive) ? segPath : prodPath + segPath;
    outAbsolutepath = outputPath.contains(dir, Qt::CaseSensitive) ? outputPath : prodPath + outputPath;

    MITK_INFO << ("[...] INPUT IMAGE: " + inputMeshFullPath).toStdString();
    MITK_INFO << ("[...] OUTPUT MESH: " + outAbsolutepath).toStdString();

    if(_useDockerContainers){
        MITK_INFO << "Using docker containers.";
        #if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
        #endif

        executableName = executablePath+"docker";
        arguments = getDockerArguments(mirtkhome.absolutePath(), commandName);

        arguments << mirtkhome.relativeFilePath(inputMeshFullPath);
        arguments << mirtkhome.relativeFilePath(outAbsolutepath);
        arguments << "-iterations" << QString::number(smth);
        arguments << "-verbose" << "3";

    } else{
        MITK_INFO << "Using static MIRTK libraries.";
        #if defined(__APPLE__)
        executablePath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
            mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
            mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
        #endif

        executablePath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
        executableName = executablePath + mitk::IOUtil::GetDirectorySeparator() + commandName;
        QDir apathd(executablePath);
        if (apathd.exists()) {
            process->setWorkingDirectory(executablePath);

            arguments << inputMeshFullPath;
            arguments << outAbsolutepath;
            arguments << "-iterations" << QString::number(smth);
            arguments << "-verbose" << "3";

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
            mitk::IOUtil::GetProgramPath();
        }//_if
    }

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutepath);

    if(!successful) {
        if (_useDockerContainers){
            MITK_WARN << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
            setUseDockerContainersOff();
            return ExecuteSmoothSurface(dir, segPath, outputPath, smth);
        } else{
            MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
            return "ERROR_IN_PROCESSING";
        }
    } else{
        return outAbsolutepath;
    }
}


/***************************************************************************
 ************************** SERVER CONC UTILITIES **************************
 ***************************************************************************/

bool CemrgCommandLine::ConnectToServer(QString userID, QString server) {

    //Setup ssh pass prompt
#if defined(__APPLE__)
    ofstream sshAskPass;
    sshAskPass.open(QDir::homePath().toStdString() + "/.ssh/ssh-askpass");
    sshAskPass << "#!/bin/bash" << endl;
    sshAskPass << "TITLE=\"${SSH_ASKPASS_TITLE:-SSH}\";" << endl;
    sshAskPass << "TEXT=\"$(whoami)'s password:\";" << endl;
    sshAskPass << "IFS=$(printf \"\\n\");" << endl;
    sshAskPass << "CODE=(\"on GetCurrentApp()\");" << endl;
    sshAskPass << "CODE=(${CODE[*]} \"tell application \\\"System Events\\\" to get short name of first process whose frontmost is true\");" << endl;
    sshAskPass << "CODE=(${CODE[*]} \"end GetCurrentApp\");" << endl;
    sshAskPass << "CODE=(${CODE[*]} \"tell application GetCurrentApp()\");" << endl;
    sshAskPass << "CODE=(${CODE[*]} \"activate\");" << endl;
    sshAskPass << "CODE=(${CODE[*]} \"display dialog \\\"${@:-$TEXT}\\\" default answer \\\"\\\" with title \\\"${TITLE}\\\" with icon caution with hidden answer\");" << endl;
    sshAskPass << "CODE=(${CODE[*]} \"text returned of result\");" << endl;
    sshAskPass << "CODE=(${CODE[*]} \"end tell\");" << endl;
    sshAskPass << "SCRIPT=\"/usr/bin/osascript\"" << endl;
    sshAskPass << "for LINE in ${CODE[*]}; do" << endl;
    sshAskPass << "\tSCRIPT=\"${SCRIPT} -e $(printf \"%q\" \"${LINE}\")\";" << endl;
    sshAskPass << "done;" << endl;
    sshAskPass << "eval \"${SCRIPT}\";" << endl;
    sshAskPass.close();
    chmod((QDir::homePath().toStdString() + "/.ssh/ssh-askpass").c_str(), S_IRUSR|S_IRGRP|S_IROTH|S_IWUSR|S_IXUSR|S_IXGRP|S_IXOTH);
    QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
    env.insert("SSH_ASKPASS", (QDir::homePath().toStdString() + "/.ssh/ssh-askpass").c_str());
    process->setProcessEnvironment(env);
#endif

    //Setup ssh config file
    ofstream sshConfigFile;
    sshConfigFile.open(QDir::homePath().toStdString() + "/.ssh/config");
    sshConfigFile << "Host " + server.toStdString() << endl;
    sshConfigFile << "ControlMaster auto" << endl;
    sshConfigFile << "ControlPath ~/.ssh/%r@%h:%p";
    sshConfigFile.close();

    //Setup connection
    QStringList arguments;
    QString connection = "ssh";
    QString usernameID = userID;
    QString serverName = server;
    arguments << usernameID + "@" + serverName;

    completion = false;
    process->start(connection, arguments);
    if (!process->waitForStarted(20000)) {
        process->close();
        return false;
    }
    if (!process->waitForReadyRead(20000)) {
        process->close();
        return false;
    }//_if_connection

    //User logged in
    process->write("mkdir ~/CEMRG-GPUReconstruction\n");
    process->write("echo\n"); process->write("echo\n");
    process->write("echo 'Festive Connection Established!'\n");
    process->write("echo\n"); process->write("echo\n");
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
        if (panel->toPlainText().contains("Festive Connection Established!")) return true;
        if (panel->toPlainText().contains("ssh Completed!")) return false;
    }//_while

    return false;
}

bool CemrgCommandLine::TransferTFServer(QString directory, QString fname, QString userID, QString server, bool download) {

    //Setup transfer command
    QStringList arguments;
    QString transfer = "scp";
    QString usernameID = userID;
    QString serverName = server;

    //Setup download/upload
    if (download == false) {

        //Clear remote host first
        arguments << usernameID + "@" + serverName;
        arguments << "rm -rf" << "~/CEMRG-GPUReconstruction/" + fname;
        process->start("ssh", arguments);
        if (!process->waitForFinished(60000)) {
            process->close();
            return false;
        }//_if_logged
        arguments.clear();
        arguments << "-r" << directory + mitk::IOUtil::GetDirectorySeparator() + fname;
        arguments << usernameID + "@" + serverName + ":~/CEMRG-GPUReconstruction";
        completion = false;
        process->start(transfer, arguments);
        if (!process->waitForStarted(20000)) {
            process->close();
            return false;
        }//_if_logged
        while (!completion) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
            if (panel->toPlainText().contains("lost connection")) return false;
            if (panel->toPlainText().contains("scp Completed!")) return true;
        }//_while

    } else {

        arguments.clear();
        arguments << usernameID + "@" + serverName + ":~/CEMRG-GPUReconstruction/" + fname;
        arguments << directory;
        completion = false;
        process->start(transfer, arguments);
        if (!process->waitForStarted(20000)) {
            process->close();
            return false;
        }//_if_logged
        while (!completion) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
            if (panel->toPlainText().contains("No such file or directory")) return false;
            if (panel->toPlainText().contains("scp Completed!")) return true;
        }//_while

    }//_if_upload
    return false;
}

void CemrgCommandLine::GPUReconstruction(QString userID, QString server, QStringList imgsList, QString targetImg, double resolution, double delta, int package, QString out) {

    //Setup remote commands
    QStringList arguments;
    QString usernameID = userID;
    QString serverName = server;
    QString cwdCommand = "cd ~/CEMRG-GPUReconstruction/Transfer;";

    //Set order of images
    QString packageList = "";
    imgsList.removeAt(imgsList.indexOf("Mask.nii.gz"));
    imgsList.removeAt(imgsList.indexOf(targetImg));
    imgsList.insert(0, targetImg);
    for (int i=0; i<imgsList.size(); i++)
        packageList = packageList + QString::number(package) + " ";

    //Setup reconstruction command
    arguments << usernameID + "@" + serverName;
    arguments << cwdCommand;
    arguments << "reconstruction_GPU2";
    arguments << "-o" << "../" + out;
    arguments << "-i" << imgsList;
    arguments << "-m" << "Mask.nii.gz";
    arguments << "-d" << QString::number(0);
    arguments << "--resolution" << QString::number(resolution);
    arguments << "--delta" << QString::number(delta);
    arguments << "--packages" << packageList.trimmed();

    completion = false;
    process->start("ssh", arguments);
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }//_while
}

//Docker - ML
QString CemrgCommandLine::dockerCemrgNetPrediction(QString mra) {

    MITK_INFO << "[CEMRGNET] Attempting prediction using Docker";
    QString aPath = "";
#if defined(__APPLE__)
    aPath = "/usr/local/bin/";
#endif

    QFileInfo finfo(mra);
    QDir cemrgnethome(finfo.absolutePath());
    QString inputfilepath = cemrgnethome.absolutePath() + mitk::IOUtil::GetDirectorySeparator() + "test.nii";
    QString tempfilepath = cemrgnethome.absolutePath() + mitk::IOUtil::GetDirectorySeparator() + "output.nii";
    QString outputfilepath = cemrgnethome.absolutePath() + mitk::IOUtil::GetDirectorySeparator() + "LA-cemrgnet.nii";

    bool test;
    if (QFile::exists(inputfilepath)) {
        MITK_INFO << "[CEMRGNET] File test.nii exists.";
        test = true;
    } else {
        MITK_INFO << "[CEMRGNET] Copying file to test.nii";
        test = QFile::copy(finfo.absoluteFilePath(), inputfilepath);
    }

    QString res;

    if(test) {

        QString inputRelativePath = cemrgnethome.relativeFilePath(inputfilepath);
        process->setWorkingDirectory(cemrgnethome.absolutePath());

        //Setup docker
        QString docker = aPath+"docker";
        QString dockerimage = "orodrazeghi/cemrgnet";
        //QString dockerexe  = "cemrgnet";
        QStringList arguments;

        arguments << "run" << "--rm";
        arguments << "--volume="+cemrgnethome.absolutePath()+":/data";
        arguments << dockerimage;
        //arguments << dockerexe;

        bool debugvar=true;
        if(debugvar) {
            MITK_INFO << "[DEBUG] Input path:";
            MITK_INFO << inputfilepath.toStdString();

            MITK_INFO << "[DEBUG] Docker command to run:";
            MITK_INFO << printFullCommand(docker, arguments);
        }

        completion = false;
        process->start(docker, arguments);
        checkForStartedProcess();
        while (!completion) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
        }
        mitk::ProgressBar::GetInstance()->Progress();

        bool test2 = QFile::rename(tempfilepath, outputfilepath);
        if(test2) {
            MITK_INFO << "[CEMRGNET] Prediction and output creation - successful.";
            res = outputfilepath;
        } else if(isOutputSuccessful(tempfilepath)) {
            MITK_INFO << "[CEMRGNET] Prediction - successful.";
            res = tempfilepath;
        } else {
            MITK_WARN << "[CEMRGNET] Problem with prediction.";
            res = "";
        }
    } else {
        MITK_WARN << "Copying input file to 'test.nii' was unsuccessful.";
        res = "";
    }//_if

    return res;
}

//Helper functions
bool CemrgCommandLine::ExecuteCommand(QString executableName, QStringList arguments, QString outputPath){
    MITK_INFO << printFullCommand(executableName, arguments);
    completion = false;
    process->start(executableName, arguments);
    bool successful = false;
    bool processStarted = checkForStartedProcess();
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();

    if (processStarted){
        successful = isOutputSuccessful(outputPath);
    }
    return successful;
}

bool CemrgCommandLine::isOutputSuccessful(QString outputfullpath) {

    MITK_INFO << "[ATTENTION] Checking for successful output on path:";
    MITK_INFO << outputfullpath.toStdString();
    QFileInfo finfo(outputfullpath);
    bool res = finfo.exists();
    MITK_INFO << (res ? "Successful output" : "Output file not found.");
    return res;
}

void CemrgCommandLine::ExecuteTouch(QString filepath) {
#ifdef _WIN32
    MITK_INFO << "[ATTENTION] touch command only necessary on macOS systems. Step ignored.";
#else
    QString commandName;
    QStringList arguments;
    commandName = "touch"; // touch filepath
    arguments << filepath;
    ExecuteCommand(commandName, arguments, filepath);
#endif
}

std::string CemrgCommandLine::printFullCommand(QString command, QStringList arguments) {

    QString argumentList = "";
    for (int ix=0; ix < arguments.size(); ix++) {
        argumentList += arguments.at(ix) + " ";
    }

    bool debugging = true;
    if (debugging) {
        QString prodPath = QString::fromStdString(mitk::IOUtil::GetProgramPath());
        MITK_INFO << ("Program path: " + prodPath).toStdString();
        ofstream prodFile1;
        prodFile1.open((prodPath + "dockerDebug.txt").toStdString(), ofstream::out | ofstream::app);
        prodFile1 << (command + " " + argumentList).toStdString() << "\n";
        prodFile1.close();
    }//_if
    return (command + " " + argumentList).toStdString();
}

/* CHECK FOR STARTED PROCESS
this function prevents freezing of the app when something goes wrong with the
Qt process.
*/
bool CemrgCommandLine::checkForStartedProcess() {
    bool startedProcess = false;
    bool debugvar = false;
    if(debugvar) {
        QStringList errinfo = QProcess::systemEnvironment();
        QString errorInfoString = "";
        for (int ix=0; ix < errinfo.size(); ix++) {
            errorInfoString += errinfo.at(ix) + " ";
        }
        MITK_INFO << "SYSTEM ENVIRONMENT:";
        MITK_INFO << errorInfoString.toStdString();
    }

    if(process->waitForStarted()) {
        MITK_INFO << "Starting process";
        startedProcess = true;
    } else {
        completion=true;
        MITK_WARN << "[ATTENTION] Process error!";
        MITK_INFO << "STATE:";
        MITK_INFO << process->state();
        MITK_INFO << "ERROR:";
        MITK_INFO << process->error();
    }
    return startedProcess;
}

void CemrgCommandLine::setUseDockerContainers(bool dockerContainersOnOff) {
    QString onoff = dockerContainersOnOff ? "ON" : "OFF";
    MITK_INFO << ("[...] Setting _useDockerContainers variable to: " + onoff).toStdString();
    _useDockerContainers = dockerContainersOnOff;
}

QStringList CemrgCommandLine::getDockerArguments(QString volume, QString dockerexe){
    bool mirtkTest = QString::compare(_dockerimage, "biomedia/mirtk:v1.1.0", Qt::CaseSensitive);

    QStringList argumentList;
    argumentList << "run" << "--rm"  << "--volume="+volume+":/data";
    argumentList << _dockerimage;
    if(mirtkTest == 0){
        argumentList << dockerexe;
    }

	return argumentList;
}

/***************************************************************************
 ************************** COMMANDLINE UTILITIES **************************
 ***************************************************************************/

void CemrgCommandLine::UpdateStdText() {

    QByteArray data = process->readAllStandardOutput();
    panel->append(QString(data));
}

void CemrgCommandLine::UpdateErrText() {

    QByteArray data = process->readAllStandardError();
    panel->append(QString(data));
}

void CemrgCommandLine::FinishedAlert() {

    completion = true;
    QString data = process->program() + " Completed!";
    panel->append(data);
}

QDialog* CemrgCommandLine::GetDialog() {

    return dial;
}
