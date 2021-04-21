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
 * Commandline Tools for MITK
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
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
#include <QCoreApplication>
#include <QDebug>

#include <thread>
#include <chrono>
#include <sys/stat.h>
#include "CemrgCommandLine.h"


CemrgCommandLine::CemrgCommandLine() {

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

    //Absolute path
    QString aPath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
    aPath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
            mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
            mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
    process->setWorkingDirectory(aPath);

    //Dilation
    QStringList arguments;
    QString input  = segPath;
    QString output = dir + mitk::IOUtil::GetDirectorySeparator() + "segmentation.d.nii";
    QString mirtk  = aPath + mitk::IOUtil::GetDirectorySeparator() + "dilate-image";
    arguments << input;
    arguments << output;
    arguments << "-iterations" << QString::number(iter);
    arguments << "-verbose" << "3";
    completion = false;
    process->start(mirtk, arguments);
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();

    //Erosion
    arguments.clear();
    input  = output;
    output = dir + mitk::IOUtil::GetDirectorySeparator() + "segmentation.s.nii";
    mirtk  = aPath + mitk::IOUtil::GetDirectorySeparator() + "erode-image";
    arguments << input;
    arguments << output;
    arguments << "-iterations" << QString::number(iter);
    arguments << "-verbose" << "3";
    completion = false;
    process->start(mirtk, arguments);
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();

    //Marching Cubes
    arguments.clear();
    input  = output;
    output = dir + mitk::IOUtil::GetDirectorySeparator() + "segmentation.vtk";
    mirtk  = aPath + mitk::IOUtil::GetDirectorySeparator() + "extract-surface";
    arguments << input;
    arguments << output;
    arguments << "-isovalue" << QString::number(th);
    arguments << "-blur" << QString::number(blur);
    arguments << "-ascii";
    arguments << "-verbose" << "3";
    completion = false;
    process->start(mirtk, arguments);
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();

    //Smoothing
    arguments.clear();
    input  = output;
    output = output;
    mirtk  = aPath + mitk::IOUtil::GetDirectorySeparator() + "smooth-surface";
    arguments << input;
    arguments << output;
    arguments << "-iterations" << QString::number(smth);
    arguments << "-verbose" << "3";
    completion = false;
    process->start(mirtk, arguments);
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();

    //Return path to output mesh
    remove((dir + mitk::IOUtil::GetDirectorySeparator() + "segmentation.d.nii").toStdString().c_str());
    remove((dir + mitk::IOUtil::GetDirectorySeparator() + "segmentation.s.nii").toStdString().c_str());
    return output;
}

QString CemrgCommandLine::ExecuteCreateCGALMesh(QString dir, QString fileName, QString templatePath) {

    //Absolute path
    QString aPath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "M3DLib";
#if defined(__APPLE__)
    aPath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
            mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
            mitk::IOUtil::GetDirectorySeparator() + QString("M3DLib");
#endif
    process->setWorkingDirectory(aPath);

    //Setup EnVariable
    QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
    env.insert("TBB_NUM_THREADS","12");
    process->setProcessEnvironment(env);

    //Setup Mesh3DTool
    QStringList arguments;
    QString input  = templatePath;
    QString output = dir + mitk::IOUtil::GetDirectorySeparator() + "CGALMeshDir";
    QString mesh3D = aPath + mitk::IOUtil::GetDirectorySeparator() + "meshtools3d";

    arguments << "-f" << input;
    arguments << "-seg_dir" << dir;
    arguments << "-seg_name" << "converted.inr";
    arguments << "-out_dir" << output;
    arguments << "-out_name" << fileName;

    completion = false;
    process->start(mesh3D, arguments);
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();

    //Return path to output CGAL mesh
    return output + mitk::IOUtil::GetDirectorySeparator() + fileName + ".vtk";
}

/***************************************************************************
 **************************** TRACKING UTILITIES ***************************
 ***************************************************************************/

void CemrgCommandLine::ExecuteTracking(QString dir, QString imgTimes, QString param) {

    //Absolute path
    QString aPath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
    aPath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
            mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
            mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
    process->setWorkingDirectory(aPath);

    //Setup
    QStringList arguments;
    QString output = dir + mitk::IOUtil::GetDirectorySeparator() + "tsffd.dof";
    QString mirtk  = aPath + mitk::IOUtil::GetDirectorySeparator() + "register";

    arguments << "-images" << imgTimes;
    if (!param.isEmpty()) arguments << "-parin" << param;
    arguments << "-dofout" << output;
    arguments << "-threads" << "12";
    arguments << "-verbose" << "3";

    completion = false;
    process->start(mirtk, arguments);
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();
}

void CemrgCommandLine::ExecuteApplying(QString dir, QString inputMesh, double iniTime, QString dofin, int noFrames, int smooth) {

    //Absolute path
    QString aPath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
    aPath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
            mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
            mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
    process->setWorkingDirectory(aPath);

    //Setup
    QStringList arguments;
    QString input  = inputMesh;
    QString output = dir + mitk::IOUtil::GetDirectorySeparator() + "transformed-";
    QString mirtk  = aPath + mitk::IOUtil::GetDirectorySeparator() + "transform-points";
    int fctTime = 10;
    noFrames *= smooth;
    if (smooth == 2)
        fctTime = 5;
    else if (smooth == 5)
        fctTime = 2;

    for (int i=0; i<noFrames; i++) {

        arguments.clear();
        arguments << input;
        arguments << output + QString::number(i) + ".vtk";
        arguments << "-dofin" << dofin;
        arguments << "-ascii";
        arguments << "-St";
        arguments << QString::number(iniTime);
        arguments << "-verbose" << "3";

        completion = false;
        process->start(mirtk, arguments);
        while (!completion) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
        }
        iniTime += fctTime;
        mitk::ProgressBar::GetInstance()->Progress();
    }
}

void CemrgCommandLine::ExecuteRegistration(QString dir, QString lge, QString mra) {

    //Absolute path
    QString aPath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
    aPath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
            mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
            mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
    process->setWorkingDirectory(aPath);

    //Setup registration
    QStringList arguments;
    QString input1 = dir + mitk::IOUtil::GetDirectorySeparator() + mra + ".nii";
    QString input2 = dir + mitk::IOUtil::GetDirectorySeparator() + lge + ".nii";
    QString output = dir + mitk::IOUtil::GetDirectorySeparator() + "rigid.dof";
    QString mirtk  = aPath + mitk::IOUtil::GetDirectorySeparator() + "register";

    arguments << input1;
    arguments << input2;
    arguments << "-dofout" << output;
    arguments << "-model" << "Rigid";
    arguments << "-verbose" << "3";

    completion = false;
    process->start(mirtk, arguments);
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();
}

void CemrgCommandLine::ExecuteTransformation(QString dir, QString imgName, QString regImgName) {

    //Absolute path
    QString aPath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
    aPath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
            mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
            mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
    process->setWorkingDirectory(aPath);

    //Setup transformation
    QStringList arguments;
    QString input  = dir + mitk::IOUtil::GetDirectorySeparator() + imgName;
    QString output = dir + mitk::IOUtil::GetDirectorySeparator() + regImgName;
    QString mirtk  = aPath + mitk::IOUtil::GetDirectorySeparator() + "transform-image";

    arguments << input;
    arguments << output;
    arguments << "-dof" << dir + mitk::IOUtil::GetDirectorySeparator() + "rigid.dof";
    arguments << "-verbose" << "3";

    completion = false;
    process->start(mirtk, arguments);
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();
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
