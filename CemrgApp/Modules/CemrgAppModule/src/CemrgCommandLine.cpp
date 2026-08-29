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

// Qmitk
#include <mitkIOUtil.h>
#include <mitkProgressBar.h>

// Qt
#include <QFileDialog>
#include <QFileInfo>
#include <QFile>
#include <QCoreApplication>
#include <QDebug>
#include <QDir>
#include <QMessageBox>

// C++ Standard
#include <thread>
#include <chrono>
#include <sys/stat.h>
#include "CemrgCommandLine.h"
#include "CemrgCommonUtils.h"

CemrgCommandLine::CemrgCommandLine() {

    _useDockerContainers = true;
    _debugvar = false;
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
    dial = new QDialog(0, Qt::WindowFlags());
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

QDialog* CemrgCommandLine::GetDialog() {

    return dial;
}

/***************************************************************************
 ****************** Execute Plugin Specific Functions **********************
 ***************************************************************************/

QString CemrgCommandLine::ExecuteSurf(QString dir, QString segPath, float thresh, int blur, int smooth) {
    MITK_INFO << "[ATTENTION] SURFACE CREATION: Surface -> Smooth";

    // Load input image into memory
    QString inputPath = segPath.contains(dir) ? segPath : dir + "/" + segPath;
    mitk::ProgressBar::GetInstance()->Progress();
    mitk::Image::Pointer inputImage = mitk::IOUtil::Load<mitk::Image>(inputPath.toStdString());

    // TODO: check semantics of thresh/blur/smooth are identical
    mitk::Surface::Pointer smoothedSurface = CemrgCommonUtils::ExtractSurfaceFromSegmentation(inputImage, (double)thresh, (double)blur, (double)smooth);
    CemrgCommonUtils::FlipXYPlane(smoothedSurface, "", "");

    // Calling code expects a QString of the file path as return value, not the Surface pointer itself.
    QString outAbsolutePath = dir + "/segmentation.vtk";
    mitk::IOUtil::Save(smoothedSurface, outAbsolutePath.toStdString());
    mitk::ProgressBar::GetInstance()->Progress(2);

    return outAbsolutePath;
}

QString CemrgCommandLine::ExecuteCreateCGALMesh(QString dir, QString outputName, QString paramsFullPath, QString segmentationName) {

    MITK_INFO << "[ATTENTION] Attempting MeshTools3D libraries.";

    QString segmentationDirectory = dir + "/";
    QString outputDirectory = segmentationDirectory + "CGALMeshDir";
    QString outAbsolutePath = outputDirectory + "/" + outputName + ".vtk"; // many outputs are created with meshtools3d. .vtk is the one used in CemrgApp

    MITK_INFO << "Using static MeshTools3D libraries.";
    QString executablePath = QCoreApplication::applicationDirPath() + "/M3DLib";
    QString executableName = executablePath + "/meshtools3d";
    QDir apathd(executablePath);
    QStringList arguments;

    if (apathd.exists()) {

        process->setWorkingDirectory(executablePath);
        arguments << "-f" << paramsFullPath;
        arguments << "-seg_dir" << segmentationDirectory;;
        arguments << "-seg_name" << segmentationName;
        arguments << "-out_dir" << outputDirectory;
        arguments << "-out_name" << outputName;

    } else {
        QMessageBox::warning(NULL, "Please check the LOG", "MeshTools3D libraries not found");
        MITK_WARN << "MeshTools3D libraries not found. Please make sure the M3DLib folder is inside the directory:\n\t" + mitk::IOUtil::GetProgramPath();
    }//_if

    //Setup EnVariable - in windows TBB_NUM_THREADS should be set in the system environment variables
#ifndef _WIN32
    QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
    env.insert("TBB_NUM_THREADS","12");
    process->setProcessEnvironment(env);
#endif

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);
    if (!successful) {
        MITK_WARN << "MeshTools3D did not produce a good outcome.";
        return "ERROR_IN_PROCESSING";
    } else {
        return outAbsolutePath;
    }
}

void CemrgCommandLine::ExecuteTracking(QString dir, QString imgTimes, QString param, QString output) {

    MITK_INFO << "[ATTENTION] Attempting Registration.";

    QString commandName = "register";
    QString imgTimesFilePath, outAbsolutePath;
    QString prodPath = dir + "/";

    imgTimesFilePath = imgTimes.contains(dir, Qt::CaseSensitive) ? imgTimes : prodPath + imgTimes;
    outAbsolutePath = output.contains(dir, Qt::CaseSensitive) ? output : prodPath + output;

    if (!outAbsolutePath.contains(".dof", Qt::CaseSensitive)) outAbsolutePath += ".dof";

    MITK_INFO << ("[...] IMAGE FILES PATH: " + imgTimesFilePath).toStdString();
    MITK_INFO << ("[...] OUTPUT DOF: " + outAbsolutePath).toStdString();

    MITK_INFO << "Using static MIRTK libraries.";
    QString executablePath = QCoreApplication::applicationDirPath() + "/MLib";
    QString executableName = executablePath + "/" + commandName;
    QDir apathd(executablePath);
    QStringList arguments;

    if (apathd.exists()) {

        process->setWorkingDirectory(executablePath);
        arguments << "-images" << imgTimesFilePath;
        if (!param.isEmpty()) arguments << "-parin" << param;
        arguments << "-dofout" << outAbsolutePath;
        arguments << "-threads" << "12";
        arguments << "-verbose" << "3";

    } else {
        QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
        MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
                        mitk::IOUtil::GetProgramPath();
    }//_if

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);
    if (!successful)
        MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
}

void CemrgCommandLine::ExecuteApplying(QString dir, QString inputMesh, double iniTime, QString dofin, int noFrames, int smooth) {

    //Time Smoothness
    int fctTime = 10;
    noFrames *= smooth;
    if (smooth == 2) {
        fctTime = 5;
    } else if (smooth == 5) {
        fctTime = 2;
    }

    QString output = dir + "/transformed-";
    QString thisOutput;
    for (int i=0; i<noFrames; i++) {
        thisOutput = output + QString::number(i)+".vtk";
        ExecuteTransformationOnPoints(dir, inputMesh, thisOutput, dofin, iniTime);
        iniTime += fctTime;
        mitk::ProgressBar::GetInstance()->Progress();
    }
}

void CemrgCommandLine::ExecuteRegistration(QString dir, QString fixed, QString moving, QString transformFileName, QString modelname) {

    MITK_INFO << "[ATTENTION] Attempting Registration.";

    //lge : fixed   ||   mra : moving
    QString commandName = "register";
    QString fixedfullpath, movingfullpath, outAbsolutePath;
    QString prodPath = dir + "/";

    fixedfullpath = fixed.contains(dir, Qt::CaseSensitive) ? fixed : prodPath + fixed;
    movingfullpath = moving.contains(dir, Qt::CaseSensitive) ? moving : prodPath + moving;
    outAbsolutePath = transformFileName.contains(dir, Qt::CaseSensitive) ? transformFileName : prodPath + transformFileName;

    if (!fixedfullpath.contains(".nii", Qt::CaseSensitive)) fixedfullpath += ".nii";
    if (!movingfullpath.contains(".nii", Qt::CaseSensitive)) movingfullpath += ".nii";
    if (!outAbsolutePath.contains(".dof", Qt::CaseSensitive)) movingfullpath += ".dof";

    MITK_INFO << ("[...] MOVING (source): " + movingfullpath).toStdString();
    MITK_INFO << ("[...] FIXED (target): " + fixedfullpath).toStdString();
    MITK_INFO << ("[...] OUTPUT DOF: " + outAbsolutePath).toStdString();

    MITK_INFO << "Using static MIRTK libraries.";
    QString executablePath = QCoreApplication::applicationDirPath() + "/MLib";
    QString executableName = executablePath + "/" + commandName;
    QDir apathd(executablePath);
    QStringList arguments;

    if (apathd.exists()) {

        process->setWorkingDirectory(executablePath);
        arguments << movingfullpath;
        arguments << fixedfullpath;
        arguments << "-dofout" << outAbsolutePath;
        arguments << "-model" << modelname;
        arguments << "-verbose" << "3";

    } else {
        QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
        MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+mitk::IOUtil::GetProgramPath();
    }//_if

    MITK_INFO << ("Performing a " + modelname + " registration").toStdString();
    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);
    if (!successful)
        MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
}

void CemrgCommandLine::ExecuteTransformation(QString dir, QString imgname, QString regname, QString transformFileFullPath) {

    MITK_INFO << "[ATTENTION] Attempting Image Transformation.";

    QString commandName = "transform-image";
    QString dofpath, imgNamefullpath, outAbsolutePath;
    QString prodPath = dir + "/";

    imgNamefullpath = imgname.contains(dir, Qt::CaseSensitive) ? imgname : prodPath + imgname;
    outAbsolutePath = regname.contains(dir, Qt::CaseSensitive) ? regname : prodPath + regname;
    dofpath = transformFileFullPath.contains(dir, Qt::CaseSensitive) ? transformFileFullPath : prodPath + transformFileFullPath;

    MITK_INFO << ("[...] INPUT IMAGE: " + imgNamefullpath).toStdString();
    MITK_INFO << ("[...] INPUT DOF: " + dofpath).toStdString();
    MITK_INFO << ("[...] OUTPUT IMAGE: " + outAbsolutePath).toStdString();

    MITK_INFO << "Using static MIRTK libraries.";
    QString executablePath = QCoreApplication::applicationDirPath() + "/MLib";
    QString executableName = executablePath + "/" + commandName;
    QDir apathd(executablePath);
    QStringList arguments;

    if (apathd.exists()) {

        process->setWorkingDirectory(executablePath);
        arguments << imgNamefullpath; //input
        arguments << outAbsolutePath; //output
        arguments << "-dofin" << dofpath;
        arguments << "-verbose" << "3";

    } else {
        QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
        MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+mitk::IOUtil::GetProgramPath();
    }//_if

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);
    if (!successful)
        MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
}

void CemrgCommandLine::ExecuteSimpleTranslation(QString dir, QString sourceMeshP, QString targetMeshP, QString transformFileName, bool transformThePoints) {

    MITK_INFO << "[ATTENTION] Attempting INIT-DOF.";

    QString executablePath, executableName, commandName, sourceMeshPath, targetMeshPath, outAbsolutePath;
    QString prodPath = dir + "/";
    QStringList arguments;

    commandName = "init-dof"; //simple translation
    sourceMeshPath = sourceMeshP.contains(dir, Qt::CaseSensitive) ? sourceMeshP : prodPath + sourceMeshP;
    targetMeshPath = targetMeshP.contains(dir, Qt::CaseSensitive) ? targetMeshP : prodPath + targetMeshP;
    outAbsolutePath = transformFileName.contains(dir, Qt::CaseSensitive) ? transformFileName : prodPath + transformFileName;

    MITK_INFO << "Using static MIRTK libraries.";
    executablePath = QCoreApplication::applicationDirPath() + "/MLib";
    executableName = executablePath + "/" + commandName;
    QDir apathd(executablePath);

    if (apathd.exists()) {

        process->setWorkingDirectory(executablePath);
        arguments << outAbsolutePath;
        arguments << "-translations" << "-norotations" << "-noscaling" << "-noshearing";
        if (transformThePoints) {
            arguments << "-displacements";
            arguments << sourceMeshPath;
            arguments << targetMeshPath;
        } else {
            arguments << "-source" << sourceMeshPath;
            arguments << "-target" << targetMeshPath;
        }
        arguments << "-verbose" << "3";

    } else {
        QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
        MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
                        mitk::IOUtil::GetProgramPath();
    }//_if

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);
    if (!successful)
        MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
}

/***************************************************************************
 ****************** Execute MIRTK Specific Functions **********************
 ***************************************************************************/

void CemrgCommandLine::ExecuteTransformationOnPoints(QString dir, QString meshFullPath, QString outputMeshFullPath, QString transformFileFullPath, double applyingIniTime) {

    MITK_INFO << "[ATTENTION] Attempting Pointset transformation.";

    QString commandName = "transform-points";
    QString dofpath, inputMeshFullPath, outAbsolutePath;
    QString prodPath = dir + "/";

    inputMeshFullPath = meshFullPath.contains(dir, Qt::CaseSensitive) ? meshFullPath : prodPath + meshFullPath;
    outAbsolutePath = outputMeshFullPath.contains(dir, Qt::CaseSensitive) ? outputMeshFullPath : prodPath + outputMeshFullPath;
    dofpath = transformFileFullPath.contains(dir, Qt::CaseSensitive) ? transformFileFullPath : prodPath + transformFileFullPath;

    MITK_INFO << ("[...] INPUT IMAGE: " + inputMeshFullPath).toStdString();
    MITK_INFO << ("[...] INPUT DOF: " + dofpath).toStdString();
    MITK_INFO << ("[...] OUTPUT IMAGE: " + outAbsolutePath).toStdString();

    MITK_INFO << "Using static MIRTK libraries.";
    QString executablePath = QCoreApplication::applicationDirPath() + "/MLib";
    QString executableName = executablePath + "/" + commandName;
    QDir apathd(executablePath);
    QStringList arguments;

    if (apathd.exists()) {

        process->setWorkingDirectory(executablePath);
        arguments << inputMeshFullPath; // input
        arguments << outAbsolutePath; // output
        arguments << "-dofin" << dofpath;
        arguments << "-ascii";
        if (applyingIniTime != -100) {
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

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);
    if (!successful)
        MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
}

void CemrgCommandLine::ExecuteResamplingOnNifti(QString niiFullPath, QString outputNiiFullPath, int isovalue) {

    MITK_INFO << "[ATTENTION] Attempting Image Transformation.";

    QString commandName = "resample-image";
    QString imgNamefullpath = niiFullPath;
    QString outAbsolutePath = outputNiiFullPath;

    MITK_INFO << ("[...] INPUT IMAGE: " + imgNamefullpath).toStdString();
    MITK_INFO << ("[...] OUTPUT IMAGE: " + outAbsolutePath).toStdString();

    MITK_INFO << "Using static MIRTK libraries.";
    QString executablePath = QCoreApplication::applicationDirPath() + "/MLib";
    QString executableName = executablePath + "/" + commandName;
    QDir apathd(executablePath);
    QStringList arguments;

    if (apathd.exists()) {

        process->setWorkingDirectory(executablePath);
        arguments << imgNamefullpath; //input
        arguments << outAbsolutePath; //output
        arguments << "-isotropic" << QString::number(isovalue);
        arguments << "-interp" << "CSpline";
        arguments << "-verbose" << "3";

    } else {
        QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
        MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
                        mitk::IOUtil::GetProgramPath();
    }//_if

    bool successful = ExecuteCommand(executableName, arguments, outAbsolutePath);
    if (!successful)
        MITK_WARN << "Local MIRTK libraries did not produce a good outcome.";
}

/***************************************************************************
 ****************** Execute Docker Specific Functions **********************
 ***************************************************************************/

QString CemrgCommandLine::DockerCemrgNetPrediction(QString mra) {

    MITK_INFO << "[CEMRGNET] Attempting prediction using Docker";

    QString executablePath;
#if defined(__APPLE__)
    executablePath = "/usr/local/bin/";
#endif
    QFileInfo finfo(mra);
    QDir cemrgnethome(finfo.absolutePath());
    QString inputfilepath = cemrgnethome.absolutePath() + "/test.nii";
    QString tempfilepath = cemrgnethome.absolutePath() + "/output.nii";
    QString outputfilepath = cemrgnethome.absolutePath() + "/LA-cemrgnet.nii";
    bool test;
    QString res;

    if (QFile::exists(inputfilepath)) {
        MITK_INFO << "[CEMRGNET] File test.nii exists.";
        test = true;
    } else {
        MITK_INFO << "[CEMRGNET] Copying file to test.nii";
        test = QFile::copy(finfo.absoluteFilePath(), inputfilepath);
    }//_if

    if (test) {

        process->setWorkingDirectory(cemrgnethome.absolutePath());

        //Setup docker
        QString docker = executablePath+"docker";
        QString dockerimage = "orodrazeghi/cemrgnet";
        QStringList arguments;
        arguments << "run" << "--rm";
        arguments << "--volume="+cemrgnethome.absolutePath()+":/data";
        arguments << dockerimage;

        if (_debugvar) {
            MITK_INFO << "[DEBUG] Input path:";
            MITK_INFO << inputfilepath.toStdString();
            MITK_INFO << "[DEBUG] Docker command to run:";
        }
        MITK_INFO << PrintFullCommand(docker, arguments);

        completion = false;
        process->start(docker, arguments);
        CheckForStartedProcess();
        while (!completion) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
        }

        bool test2 = QFile::rename(tempfilepath, outputfilepath);
        if (test2) {
            MITK_INFO << "[CEMRGNET] Prediction and output creation - successful.";
            res = outputfilepath;
        } else if (IsOutputSuccessful(tempfilepath)) {
            MITK_INFO << "[CEMRGNET] Prediction - successful.";
            res = tempfilepath;
        } else {
            MITK_WARN << "[CEMRGNET] Problem with prediction.";
            res = "";
        }//_if

    } else {
        MITK_WARN << "Copying input file to 'test.nii' was unsuccessful.";
        res = "";
    }//_if

    return res;
}

QString CemrgCommandLine::DockerDicom2Nifti(QString path2dicomfolder) {

    MITK_INFO << "[ATTENTION] Attempting alternative DICOM to NIFTI conversion.";

    QDir dicomhome(path2dicomfolder);
    QString outAbsolutePath = "ERROR_IN_PROCESSING";
    QString executablePath;
    QString outPath;
    bool successful = false;

    MITK_INFO << ("[...] OUTPUT IMAGE: " + outAbsolutePath).toStdString();

    if (_useDockerContainers) {

        MITK_INFO << "Using docker containers.";
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
        outPath = dicomhome.absolutePath() + "/NIIs";
        QString executableName = executablePath+"docker";

        QStringList arguments;
        arguments << "run" << "--rm"  << "--volume="+dicomhome.absolutePath()+":/Data";
        arguments << "orodrazeghi/dicom-converter" << ".";
        arguments << "--gantry" << "--inconsistent";

        successful = ExecuteCommand(executableName, arguments, outPath, false);

    } else {
        MITK_WARN << "Docker must be running for this feature to be used.";
    }//_if

    if (successful) {
        MITK_INFO << "Conversion successful.";
        outAbsolutePath = outPath;
    } else {
        MITK_WARN << "Error with DICOM2NIFTI Docker container.";
    }

    return outAbsolutePath;
}

QString CemrgCommandLine::DockerUniversalAtrialCoordinates(QString dir, QString uaccmd, QStringList fibreAtlas, QString meshname, QStringList cmdargs, QStringList landmarks, QString outnameext){
    SetDockerImageUac();
    QString executablePath;
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
    QString executableName = executablePath+"docker";
    QString outAbsolutePath = "ERROR_IN_PROCESSING";

    QDir home(dir);
    QStringList arguments = GetDockerArguments(home.absolutePath());

    arguments << uaccmd;
    arguments << fibreAtlas; // append list
    arguments << meshname;
    arguments << cmdargs; // append list

    for (int ix = 0; ix < landmarks.size(); ix++) {
        arguments << home.relativeFilePath(landmarks.at(ix));
    }

    // output filename checked when running ExecuteCommand
    QString outPath = home.absolutePath() + "/" + outnameext;
    bool successful = ExecuteCommand(executableName, arguments, outPath, outnameext.isEmpty());

    if(successful){
        MITK_INFO << ("UAC command: " + uaccmd + " successful").toStdString();
        outAbsolutePath = outPath;
    } else{
        MITK_WARN << ("Error running UAC command: " + uaccmd).toStdString();
    }

    return outAbsolutePath;
}

QString CemrgCommandLine::DockerUacMainMode(QString dir, QString stage, QString atrium, QString layer, QString fibre, QString meshname, QStringList tags, QStringList landmarks, bool fourch, bool noraa, int scale){
   SetDockerImageUac("3.0-beta");
   QString executablePath;
#if defined(__APPLE__)
   executablePath = "/usr/local/bin/";
#endif
   QString executableName = executablePath + "docker";
   QString outAbsolutePath = "ERROR_IN_PROCESSING";

   QDir home(dir);
   QStringList arguments = GetDockerArguments(home.absolutePath());

   QString uaccmd = "uac";
   arguments << uaccmd;
   arguments << "--uac-stage" << stage;
   arguments << "--atrium" << atrium;
   arguments << "--layer" << layer;
   arguments << "--fibre" << fibre;
   arguments << "--msh" << meshname;

   arguments << "--tags";
   for (int ix = 0; ix < tags.size(); ix++) {
       arguments << tags.at(ix);
   }

   if (landmarks.size() > 0) {
       arguments << "--landmarks" << home.relativeFilePath(landmarks.at(0));
       if (landmarks.size() > 1) {
           arguments << "--regions" << home.relativeFilePath(landmarks.at(1));
       }
   }

   if (fourch) {
       arguments << "--fourch";
   }

   if (noraa) {
       arguments << "--noraa";
   }

   arguments << QString::number(scale);

   QStringList outputs;
   if (stage.compare("1")) {
       outputs << "LSbc1.vtx" << "LSbc2.vtx" << "PAbc1.vtx" << "PAbc2.vtx";

   } else if (stage.compare("2a")){
       outputs << "AnteriorMesh.elem" << "AnteriorMesh.pts" << "PosteriorMesh.elem" << "PosteriorMesh.pts";

   } else if (stage.compare("2b")){
       outputs << "Labelled_Coords_2D_Rescaling_v3_C.elem" << "Labelled_Coords_2D_Rescaling_v3_C.pts";
   }

   QString outPath = home.absolutePath() + "/" + outputs.at(0);
   bool successful = ExecuteCommand(executableName, arguments, outPath);

    if(successful){
        MITK_INFO << ("UAC command: " + uaccmd + " successful").toStdString();
        outAbsolutePath = outPath;
    } else{
        MITK_WARN << ("Error running UAC command: " + uaccmd).toStdString();
    }

    return outAbsolutePath;

}

QString CemrgCommandLine::DockerUacFibreMappingMode(QString dir, QString atrium, QString layer, QString fibre, QString meshname, bool msh_endo_epi, QString output, bool fourch, QString tags, QString biproj){
   SetDockerImageUac("3.0-beta");
   QString executablePath;
#if defined(__APPLE__)
   executablePath = "/usr/local/bin/";
#endif
   QString executableName = executablePath + "docker";
   QString outAbsolutePath = "ERROR_IN_PROCESSING";

   QDir home(dir);
   QStringList arguments = GetDockerArguments(home.absolutePath());


   QString uaccmd = "fibremap";

   arguments << uaccmd ;
   arguments << "--atrium" << atrium;
   arguments << "--layer" << layer;
   arguments << "--fibre" << fibre;
   arguments << "--msh" << meshname;
   arguments << "--output" << output;

   if (fourch) {
       arguments << "--fourch";
   }

   if (msh_endo_epi){
       arguments << "--msh-endo" << "Labelled";
       arguments << "--msh-epi" << "Labelled";
   }

   arguments << "--tags" << tags;
   arguments << "--fibre-biproj" << biproj;

   QString omsh = output;
   omsh += layer.contains("bilayer") ? "Bilayer" : "";

   QStringList outputs;
   outputs << omsh+".pts" << omsh+".elem";

   QString outPath = home.absolutePath() + "/" + outputs.at(0);
   bool successful = ExecuteCommand(executableName, arguments, outPath);

    if(successful){
        MITK_INFO << ("UAC command: " + uaccmd + " successful").toStdString();
        outAbsolutePath = outPath;
    } else{
        MITK_WARN << ("Error running UAC command: " + uaccmd).toStdString();
    }

    return outAbsolutePath;
}

QString CemrgCommandLine::DockerSurfaceFromMesh(QString dir, QString meshname, QString outname, QString op, QString outputSuffix){
    // Method equivalent to:  meshtool extract surface
    SetDockerImage("alonsojasl/cemrg-meshtool:v1.0");
    QString executablePath;
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
    QString executableName = executablePath+"docker";
    QString outAbsolutePath = "ERROR_IN_PROCESSING";

    QDir home(dir);

    outname += (outputSuffix.at(0)=="_") ? outputSuffix : ("_"+outputSuffix);

    QStringList arguments = GetDockerArguments(home.absolutePath());
    arguments << "extract" << "surface";
    arguments << ("-msh="+meshname);
    arguments << ("-op="+op);
    arguments << ("-surf="+outname);

    QString outPath = home.absolutePath() + "/" + outname + ".surf.vtx";

    bool successful = ExecuteCommand(executableName, arguments, outPath);

    if (successful) {
        MITK_INFO << "Surface extraction successful.";
        outAbsolutePath = outPath;
    } else{
        MITK_WARN << "Error with MESHTOOL Docker container.";
    }
    return outAbsolutePath;
}

QString CemrgCommandLine::DockerExtractGradient(QString dir, QString meshname, QString idatName, QString odatName, bool elemGrad){
    // Method equivalent to:  meshtool extract surface
    SetDockerImage("alonsojasl/cemrg-meshtool:v1.0");
    QString executablePath;
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
    QString executableName = executablePath+"docker";
    QString outAbsolutePath = "ERROR_IN_PROCESSING";

    QDir home(dir);

    QStringList arguments = GetDockerArguments(home.absolutePath());
    arguments << "extract" << "gradient";
    arguments << ("-msh="+meshname);
    arguments << ("-idat="+home.relativeFilePath(idatName));
    arguments << ("-odat="+home.relativeFilePath(odatName));
    if(elemGrad){
        arguments << "-mode=1"; // compute element gradient
    }

    QString outPath = home.absolutePath() + "/" + odatName + ".grad.vec";

    bool successful = ExecuteCommand(executableName, arguments, outPath);

    if (successful) {
        MITK_INFO << "Gradient extraction successful.";
        outAbsolutePath = outPath;
    } else{
        MITK_WARN << "Error with MESHTOOL Docker container.";
    }
    return outAbsolutePath;
}

QString CemrgCommandLine::DockerRemeshSurface(QString dir, QString meshname, QString outname, double hmax, double hmin, double havg, double surfCorr){
    // Method equivalent to: meshtool resample surfmesh
    SetDockerImage("alonsojasl/cemrg-meshtool:v1.0");
    QString executablePath;
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
    QString executableName = executablePath+"docker";
    QString outAbsolutePath = "ERROR_IN_PROCESSING";

    QDir home(dir);

    QStringList arguments = GetDockerArguments(home.absolutePath());
    arguments << "resample" << "surfmesh";
    arguments << ("-msh="+meshname);
    arguments << ("-ifmt=vtk");
    arguments << ("-outmsh="+outname);
    arguments << ("-ofmt=vtk_polydata");
    arguments << ("-max="+QString::number(hmax));
    arguments << ("-avrg="+QString::number(havg));
    arguments << ("-min="+QString::number(hmin));

    if(surfCorr>0){
        arguments << ("-surf_corr="+QString::number(surfCorr));
    }

    QString outPath = home.absolutePath() + "/" + outname + ".vtk";

    bool successful = ExecuteCommand(executableName, arguments, outPath);

    if (successful) {
        MITK_INFO << "Surface remeshing successful.";
        outAbsolutePath = outPath;
    } else{
        MITK_WARN << "Error with MESHTOOL Docker container.";
    }
    return outAbsolutePath;
}

QString CemrgCommandLine::DockerInterpolateData(QString dir, QString meshname, QString outmesh, QString idatExt, QString odatExt, QString dataType){
    // Method equivalent to: meshtool interpolate dataType
    SetDockerImage("alonsojasl/cemrg-meshtool:v1.0");
    QString executablePath = "";
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
    QString executableName = executablePath+"docker";
    QString outAbsolutePath = "ERROR_IN_PROCESSING";

    QDir home(dir);
    if(!dataType.contains("elemdata") && !dataType.contains("nodedata") && !dataType.contains("clouddata")){
        MITK_ERROR << "Incorrect parameter seleted";
    } else{
        QStringList arguments = GetDockerArguments(home.absolutePath());
        arguments << "interpolate" << dataType;
        if(dataType.contains("clouddata")){
            arguments << ("-pts="+meshname);
        } else{
            arguments << ("-imsh="+meshname);
        }
        arguments << ("-omsh="+outmesh);
        arguments << ("-idat="+idatExt);
        arguments << ("-odat="+odatExt);

        QString outPath = home.absolutePath() + "/" + odatExt;

        bool successful = ExecuteCommand(executableName, arguments, outPath);

        if (successful) {
            MITK_INFO << "Interpolating data successful.";
            outAbsolutePath = outPath;
        } else{
            MITK_WARN << "Error with MESHTOOL Docker container.";
        }
    }
    return outAbsolutePath;
}

QString CemrgCommandLine::DockerConvertMeshFormat(QString dir, QString imsh, QString ifmt, QString omsh, QString ofmt, double scale){
    // Method equivalent to: meshtool convert
    SetDockerImage("alonsojasl/cemrg-meshtool:v1.0");
    QString executablePath = "";
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
    QString executableName = executablePath+"docker";
    QString outAbsolutePath = "ERROR_IN_PROCESSING";

    QDir home(dir);

    QStringList arguments = GetDockerArguments(home.absolutePath());
    arguments << "convert";
    arguments << ("-imsh="+imsh);
    arguments << ("-ifmt="+ifmt);
    arguments << ("-omsh="+omsh);
    arguments << ("-ofmt="+ofmt);

    if(scale>0){
        arguments << ("-scale="+QString::number(scale));
    }


    QString outPath = home.absolutePath() + "/" + omsh;
    bool isConvertToCarp = ofmt.contains("carp", Qt::CaseInsensitive);
    outPath += (isConvertToCarp) ? ".pts" : ".vtk";

    bool successful = ExecuteCommand(executableName, arguments, outPath, !isConvertToCarp);

    if (successful) {
        MITK_INFO << "Surface remeshing successful.";
        outAbsolutePath = outPath;
    } else{
        MITK_WARN << "Error with MESHTOOL Docker container.";
    }
    return outAbsolutePath;
}

void CemrgCommandLine::DockerCleanMeshQuality(QString dir, QString meshname, QString outMesh, double qualityThres, QString ifmt, QString ofmt){
    // Method equivalent to: meshtool clean quality
    SetDockerImage("alonsojasl/cemrg-meshtool:v1.0");
    QString executablePath = "";
#if defined(__APPLE__)
        executablePath = "/usr/local/bin/";
#endif
    QString executableName = executablePath+"docker";


    QDir home(dir);

    double smth=0.75;
    int iter=200;

    QStringList arguments = GetDockerArguments(home.absolutePath());
    arguments << "clean" << "quality";
    arguments << ("-msh="+meshname);
    arguments << ("-ifmt="+ifmt);
    arguments << ("-outmsh="+outMesh);
    arguments << ("-ofmt="+ofmt);
    arguments << ("-thr="+QString::number(qualityThres));
    arguments << ("-smth="+QString::number(smth));
    arguments << ("-iter="+QString::number(iter));


    QString outPath = home.absolutePath() + "/" + outMesh;
    outPath += (ofmt.contains("carp", Qt::CaseInsensitive)) ? ".pts" : ".vtk";

    bool successful = ExecuteCommand(executableName, arguments, outPath);

    if (successful) {
        MITK_INFO << "Surface remeshing successful.";

    } else{
        MITK_WARN << "Error with MESHTOOL Docker container.";
    }
    // return outAbsolutePath;
}


/***************************************************************************
 *********************** Docker Helper Functions ***************************
 ***************************************************************************/

void CemrgCommandLine::SetUseDockerContainers(bool dockerContainersOnOff) {

    QString onoff = dockerContainersOnOff ? "ON" : "OFF";
    MITK_INFO << ("[...] Setting _useDockerContainers variable to: " + onoff).toStdString();
    _useDockerContainers = dockerContainersOnOff;
}

QStringList CemrgCommandLine::GetDockerArguments(QString volume, QString dockerexe) {

    bool mirtkTest = QString::compare(_dockerimage, "biomedia/mirtk:v1.1.0", Qt::CaseSensitive);
    QStringList argumentList;
    argumentList << "run" << "--rm"  << "--volume="+volume+":/data";
    argumentList << _dockerimage;
    if (mirtkTest == 0)
        argumentList << dockerexe;
    return argumentList;
}

QStringList CemrgCommandLine::GetOpenCarpDockerDefaultArguments(QString volume){

    QString petscPath = QCoreApplication::applicationDirPath() + "/petsc_opts";
    QString parab="ilu_cg_opts", ellip="amg_cg_opts";

    QString parabFile = petscPath+"/"+parab;
    QString ellipFile = petscPath+"/"+ellip;

    QString parabDestination = volume+"/"+parab;
    QString ellipDestination = volume+"/"+ellip;

    if(!QFile::exists(parabDestination)){
        MITK_INFO << ("Copying: ["+parabFile+ "]").toStdString();
        MITK_INFO << ("Into: ["+parabDestination+ "]").toStdString();
        MITK_INFO(QFile::copy(parabFile, parabDestination)) << "Success!";
    }

    if(!QFile::exists(ellipDestination)){
        MITK_INFO << "Copying: [" << ellipFile.toStdString();
        MITK_INFO << "Into: [" << ellipDestination.toStdString();
        MITK_INFO(QFile::copy(ellipFile, ellipDestination)) << "Success!";
    }

    QDir home(volume);
    QStringList defaultArguments;
    defaultArguments << "run" << "--rm" << ("--volume="+volume+":/shared:z") << "--workdir=/shared";
    defaultArguments << "docker.opencarp.org/opencarp/opencarp:latest";
    defaultArguments << "openCARP";
    defaultArguments << "-ellip_use_pt" << "0" << "-parab_use_pt" << "0";
    defaultArguments << "-parab_options_file";
    defaultArguments << home.relativeFilePath(parab);
    defaultArguments << "-ellip_options_file";
    defaultArguments << home.relativeFilePath(ellip);

    return defaultArguments;
}

QString CemrgCommandLine::OpenCarpDockerLaplaceSolves(QString dir, QString meshName, QString outName, QStringList zeroName, QStringList oneName, QStringList regionLabels){
    SetDockerImage("docker.opencarp.org/opencarp/opencarp:latest");
    QString executablePath;
    #if defined(__APPLE__)
            executablePath = "/usr/local/bin/";
    #endif
        QString executableName = executablePath+"docker";
        QString outAbsolutePath = "ERROR_IN_PROCESSING";

        QDir home(dir);
        QString outPath = home.absolutePath() + "/fibres_pt_1/" + outName;
        QDir outDir(outPath);
        QString outIgbFile = outPath + "/phie.igb";

        MITK_INFO(outDir.mkpath(outPath)) << "Output directory created.";
        if(!outDir.exists()){
            MITK_INFO << ("Error creating directory: " + outPath).toStdString();
        } else{
            QStringList arguments = GetOpenCarpDockerDefaultArguments(home.absolutePath());
            arguments << "-simID" << home.relativeFilePath(outPath);
            arguments << "-meshname" << meshName;
            arguments << "-experiment" << "2";
            arguments << "-bidomain" << "1";

            MITK_INFO << "[openCarp] setup arguments computed";

            int numStim = zeroName.size() + oneName.size();
            arguments << "-num_stim" << QString::number(numStim);
            for (int ix = 0; ix < zeroName.size(); ix++) {
                arguments << ("-stimulus[" + QString::number(ix) + "].vtx_file") << home.relativeFilePath(zeroName.at(ix));
                arguments << ("-stimulus[" + QString::number(ix) + "].stimtype") << "3"; // =0
            }
            for (int jx = zeroName.size(); jx < numStim; jx++) {
                arguments << ("-stimulus[" + QString::number(jx) + "].vtx_file") << home.relativeFilePath(oneName.at(jx-zeroName.size()));
                arguments << ("-stimulus[" + QString::number(jx) + "].stimtype") << "2"; // =1
            }

            arguments << ("-stimulus[" + QString::number(numStim-1) + "].start") << "0";
            arguments << ("-stimulus[" + QString::number(numStim-1) + "].duration") << "1";
            arguments << ("-stimulus[" + QString::number(numStim-1) + "].strength") << "1";

            MITK_INFO << "[openCarp] stimulus arguments computed";

            arguments << "-num_gregions" << "1";
            arguments << "-gregion[0].name" << "myo";
            arguments << "-gregion[0].g_il" << "1";
            arguments << "-gregion[0].g_it" << "1";
            arguments << "-gregion[0].g_in" << "1";
            arguments << "-gregion[0].g_el" << "1";
            arguments << "-gregion[0].g_et" << "1";
            arguments << "-gregion[0].g_en" << "1";
            arguments << "-gregion[0].num_IDs"<< QString::number(regionLabels.size());
            for (int ix = 0; ix < regionLabels.size(); ix++) {
                arguments << "-gregion[0].ID[" +QString::number(ix)+ "]"<< regionLabels.at(ix);
            }

            bool successful = ExecuteCommand(executableName, arguments, outIgbFile);

            if (successful) {
                MITK_INFO << "Laplace solves generation successful. Creating .dat file";
                QString outPathFile = "/" + meshName + "_" + outName + "_potential.dat";

                arguments.clear();
                arguments << "run" << "--rm" << ("--volume="+home.absolutePath()+":/shared:z") << "--workdir=/shared";
                arguments << "docker.opencarp.org/opencarp/opencarp:latest";
                arguments << "igbextract" << home.relativeFilePath(outIgbFile) << "-O";
                arguments << home.relativeFilePath(outPathFile) << "-o" << "ascii_1pL";
                successful = ExecuteCommand(executableName, arguments, outPathFile);
                if(successful){
                    outAbsolutePath =  outPathFile;
                }
            } else{
                MITK_WARN << "Error with openCARP LAPLACE SOLVES Docker container.";
            }
        }

        return outAbsolutePath;
}

QString CemrgCommandLine::OpenCarpDocker(QString dir, QString paramfile, QString simID){
    SetDockerImage("docker.opencarp.org/opencarp/opencarp:latest");
    QString executablePath;
    #if defined(__APPLE__)
            executablePath = "/usr/local/bin/";
    #endif
        QString executableName = executablePath+"docker";
        QString outAbsolutePath = "ERROR_IN_PROCESSING";

        QDir home(dir);
        QString outPath = home.absolutePath() + "/" + simID;
        QString outPhieFilePath = outPath; // + "/phie.igb";
        QDir outDir(outPath);

        MITK_INFO(outDir.mkpath(outPath)) << "Output directory created.";
        if(!outDir.exists()){
            MITK_INFO << ("Error creating directory: " + outPath).toStdString();
        } else{
            QStringList arguments = GetOpenCarpDockerDefaultArguments(home.absolutePath());

            arguments << "+F" << home.relativeFilePath(paramfile);
            arguments << "-simID" << home.relativeFilePath(outPath);

            bool successful = ExecuteCommand(executableName, arguments, outPhieFilePath, false);
            if (successful) {
                outAbsolutePath =  outPhieFilePath;
            } else{
                MITK_WARN << "Error with openCARP LAPLACE SOLVES Docker container.";
            }
        }

        return outAbsolutePath;
}

/***************************************************************************
 **************************** Helper Functions *****************************
 ***************************************************************************/

bool CemrgCommandLine::CheckForStartedProcess() {

    //CHECK FOR STARTED PROCESS
    //This function prevents freezing of the app when something goes wrong with the Qt process.
    bool startedProcess = false;

    if (_debugvar) {
        QStringList errinfo = QProcess::systemEnvironment();
        QString errorInfoString = "";
        for (int ix=0; ix < errinfo.size(); ix++)
            errorInfoString += errinfo.at(ix) + " ";
        MITK_INFO << "SYSTEM ENVIRONMENT:";
        MITK_INFO << errorInfoString.toStdString();
    }

    if (process->waitForStarted()) {

        MITK_INFO << "Starting process";
        startedProcess = true;

    } else {

        completion = true;
        MITK_WARN << "[ATTENTION] Process error!";
        MITK_INFO << "STATE:";
        MITK_INFO << process->state();
        MITK_INFO << "ERROR:";
        MITK_INFO << process->error();

    }//_if

    return startedProcess;
}

void CemrgCommandLine::ExecuteTouch(QString filepath) {

#ifdef _WIN32
    MITK_INFO << "[ATTENTION] touch command is not necessary on Windows systems. Step ignored.";
#else

    QString commandName;
    QStringList arguments;
    commandName = "touch"; // touch filepath
    arguments << filepath;

    completion = false;
    process->start(commandName, arguments);
    bool processStarted = CheckForStartedProcess();
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    MITK_INFO(!processStarted) << "[ATTENTION] TOUCH Process never started.";

#endif
}

bool CemrgCommandLine::IsOutputSuccessful(QString outputFullPath) {

    MITK_INFO << "[ATTENTION] Checking for successful output on path:";
    MITK_INFO << outputFullPath.toStdString();

    QFileInfo finfo(outputFullPath);
    bool fileExists = finfo.exists();
    bool result = false;

    MITK_INFO << (fileExists ? "File exists." : "Output file not found.");

    if (fileExists) {
        bool fileSizeTest = finfo.size() > 0;
        if (fileSizeTest) {
            MITK_INFO << ("File size: " + QString::number(finfo.size())).toStdString();
            result = true;
        } else {
            MITK_INFO << "File empty. Output unsuccessful";
        }
    }//_if

    return result;
}

std::string CemrgCommandLine::PrintFullCommand(QString command, QStringList arguments) {
    SetDebugOn();
    QString argumentList = "";
    for (int ix=0; ix < arguments.size(); ix++)
        argumentList += arguments.at(ix) + " ";

    if (_debugvar) {
        QString prodPath = QString::fromStdString(mitk::IOUtil::GetProgramPath());
        MITK_INFO << ("Program path: " + prodPath).toStdString();
        std::ofstream prodFile1;
        prodFile1.open((prodPath + "dockerDebug.txt").toStdString(), std::ofstream::out | std::ofstream::app);
        prodFile1 << (command + " " + argumentList).toStdString() << "\n";
        prodFile1.close();
    }//_if
    SetDebugOff();

    return (command + " " + argumentList).toStdString();
}

bool CemrgCommandLine::ExecuteCommand(QString executableName, QStringList arguments, QString outputPath, bool isOutputFile) {

    MITK_INFO << PrintFullCommand(executableName, arguments);

    if(isOutputFile){ // if false, the output is a folder and does not need touch
        MITK_INFO << ("[ExecuteCommand] Creating empty file at output:" + outputPath).toStdString();
        ExecuteTouch(outputPath);
    }

    completion = false;
    process->start(executableName, arguments);
    bool successful = false;
    bool processStarted = CheckForStartedProcess();

    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }

    if (processStarted)
        successful = IsOutputSuccessful(outputPath);

    return successful;
}

/***************************************************************************
 ************************** Protected Slots ********************************
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
