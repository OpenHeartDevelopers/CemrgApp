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
    useDockerContainers = true;

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
    useDockerContainers = true;

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

    QString retOutput;
    QString dockerOutput = dockerSurf(dir, segPath, iter, th, blur, smth);
    bool successful = useDockerContainers ? isOutputSuccessful(dockerOutput) : false;

    if(!useDockerContainers || !successful){
        MITK_INFO(!useDockerContainers) << "Using static MIRTK libraries.";
        MITK_WARN(useDockerContainers && (!successful)) << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";

        //Absolute path
        QString aPath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
        aPath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
        process->setWorkingDirectory(aPath);

        QDir apathd(aPath);
        if (apathd.exists()) {
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
            output = output+"";
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
            retOutput = output;
        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+  mitk::IOUtil::GetProgramPath();
            retOutput = "";
        }
    } else {
        retOutput = dockerOutput;
    }
    return retOutput;
}

QString CemrgCommandLine::ExecuteCreateCGALMesh(QString dir, QString fileName, QString templatePath) {
    bool debugvar = false;
    MITK_INFO(debugvar) << ("DIRECTORY: " + dir).toStdString();
    MITK_INFO(debugvar) << ("FILENAME: " + fileName).toStdString();
    MITK_INFO(debugvar) << ("TEMPLATE: " + templatePath).toStdString();

    QString dockerOutput = dockerCreateCGALMesh(dir, fileName, templatePath);
    QString retOutput;
    bool successful = useDockerContainers ? isOutputSuccessful(dockerOutput) : false;

    if(!useDockerContainers || !successful){
        MITK_INFO(!useDockerContainers) << "Using static MESHTOOLS3D libraries.";
        MITK_WARN(useDockerContainers && (!successful)) << "Docker did not produce a good outcome. Trying with local MESHTOOLS3D libraries.";
        //Absolute path
        QString aPath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "M3DLib";
#if defined(__APPLE__)
        aPath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("M3DLib");
#endif

        //Setup EnVariable - in windows TBB_NUM_THREADS should be set in the system environment variables
#ifndef _WIN32
        //	QString setenv = "set TBB_NUM_THREADS=4";
        //	process->start(setenv);
        //#else
        QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
        env.insert("TBB_NUM_THREADS","12");
        process->setProcessEnvironment(env);
#endif

        QDir apathd(aPath);
        if(apathd.exists()) {
            process->setWorkingDirectory(aPath);
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
            retOutput = output + mitk::IOUtil::GetDirectorySeparator() + fileName + ".vtk";

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "Meshtools3D libraries not found");
            MITK_WARN << "Meshtools3D libraries not found. Please make sure the MLib folder is inside the directory;\n\t" +
                         mitk::IOUtil::GetProgramPath();
            retOutput = "";
        }
    } else{
        retOutput = dockerOutput;
    }
    return retOutput;
}

/***************************************************************************
 **************************** TRACKING UTILITIES ***************************
 ***************************************************************************/

void CemrgCommandLine::ExecuteTracking(QString dir, QString imgTimes, QString param) {

    bool successful = useDockerContainers ? dockerTracking(dir, imgTimes, param) : false;

    if(!useDockerContainers || !successful) {
        MITK_INFO(!useDockerContainers) << "Using static MIRTK libraries.";
        MITK_WARN(useDockerContainers && (!successful)) << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
        //Absolute path
        QString aPath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
        aPath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif

        QDir apathd(aPath);
        if(apathd.exists()) {
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

            MITK_INFO << printFullCommand(mirtk, arguments);
            process->start(mirtk, arguments);
            while (!completion) {
                std::this_thread::sleep_for(std::chrono::seconds(1));
                QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
            }
            mitk::ProgressBar::GetInstance()->Progress();

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
                         mitk::IOUtil::GetProgramPath();
        }//_if
    }
}

void CemrgCommandLine::ExecuteApplying(QString dir, QString inputMesh, double iniTime, QString dofin, int noFrames, int smooth) {

    bool successful = useDockerContainers ? dockerApplying(dir, inputMesh, iniTime, dofin, noFrames, smooth) : false;

    if(!useDockerContainers || !successful) {
        MITK_INFO(!useDockerContainers) << "Using static MIRTK libraries.";
        MITK_WARN(useDockerContainers && (!successful)) << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
        //Absolute path
        QString aPath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
        aPath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif

        QDir apathd(aPath);
        if(apathd.exists()) {
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
        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
                         mitk::IOUtil::GetProgramPath();
        }//_if
    }
}

void CemrgCommandLine::ExecuteRegistration(QString dir, QString fixed, QString moving, QString txname, QString modelname) {

    //lge : fixed
    //mra : moving
    QString fixedfullpath, movingfullpath, output;
    QString prodPath = dir + mitk::IOUtil::GetDirectorySeparator();

    fixedfullpath = fixed.contains(dir, Qt::CaseSensitive) ? fixed : prodPath + fixed;
    movingfullpath = moving.contains(dir, Qt::CaseSensitive) ? moving : prodPath + moving;
    output = txname.contains(dir, Qt::CaseSensitive) ? txname : prodPath + txname;

    if(!fixedfullpath.contains(".nii", Qt::CaseSensitive)) fixedfullpath += ".nii";
    if(!movingfullpath.contains(".nii", Qt::CaseSensitive)) movingfullpath += ".nii";
    if(!output.contains(".dof", Qt::CaseSensitive)) movingfullpath += ".dof";

    bool successful = useDockerContainers ? dockerRegistration(dir, fixedfullpath, movingfullpath, output, modelname) : false;

    if(!useDockerContainers || !successful) {
        MITK_INFO(!useDockerContainers) << "Using static MIRTK libraries.";
        MITK_WARN(useDockerContainers && (!successful)) << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
        //Absolute path
        QString aPath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
        aPath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
        QDir apathd(aPath);
        if (apathd.exists()) {
            process->setWorkingDirectory(aPath);

            //Setup registration
            QStringList arguments;
            QString input1 = movingfullpath;
            QString input2 = fixedfullpath;
            QString mirtk  = aPath + mitk::IOUtil::GetDirectorySeparator() + "register";

            arguments << input1;
            arguments << input2;
            arguments << "-dofout" << output;
            arguments << "-model" << modelname;
            arguments << "-verbose" << "3";

            completion = false;
            printFullCommand(mirtk, arguments);
            process->start(mirtk, arguments);
            MITK_INFO << "Executing a " + modelname.toStdString() + " registration.";
            while (!completion) {
                std::this_thread::sleep_for(std::chrono::seconds(1));
                QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
            }
            mitk::ProgressBar::GetInstance()->Progress();

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
                         mitk::IOUtil::GetProgramPath();
        }//_if
    }
}

void CemrgCommandLine::ExecuteTransformation(QString dir, QString imgname, QString regname, QString txfullpath) {

    QString dofpath, imgNamefullpath, regImgNamefullpath;
    QString prodPath = dir + mitk::IOUtil::GetDirectorySeparator();

    imgNamefullpath = imgname.contains(dir, Qt::CaseSensitive) ? imgname : prodPath + imgname;
    regImgNamefullpath = regname.contains(dir, Qt::CaseSensitive) ? regname : prodPath + regname;
    dofpath = txfullpath.contains(dir, Qt::CaseSensitive) ? txfullpath : prodPath + txfullpath;

    bool successful = useDockerContainers ? dockerTranformation(dir, imgNamefullpath, regImgNamefullpath, dofpath) : false;

    if(!useDockerContainers || !successful) {
        MITK_INFO(!useDockerContainers) << "Using static MIRTK libraries.";
        MITK_WARN(useDockerContainers && (!successful)) << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
        //Absolute path
        QString aPath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
        aPath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
        QDir apathd(aPath);
        if (apathd.exists()) {
            process->setWorkingDirectory(aPath);

            //Setup transformation
            QStringList arguments;
            QString input  = imgNamefullpath;
            QString output = regImgNamefullpath;
            QString mirtk  = aPath + mitk::IOUtil::GetDirectorySeparator() + "transform-image";

            arguments << input;
            arguments << output;
            arguments << "-dof" << dofpath;
            arguments << "-verbose" << "3";

            completion = false;
            MITK_INFO << printFullCommand(mirtk, arguments);
            process->start(mirtk, arguments);
            while (!completion) {
                std::this_thread::sleep_for(std::chrono::seconds(1));
                QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
            }
            mitk::ProgressBar::GetInstance()->Progress();

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
                         mitk::IOUtil::GetProgramPath() + "\nLooked for folder in: \n\t" + aPath.toStdString();
        }//_if
    }
}

void CemrgCommandLine::ExecuteResamplingOnNifti(
        QString niifullpath, QString outputtniifullpath, int isovalue) {

    // /resample-image niifullpath outputtniifullpath -isotropic 0.5 -interp CSpline -verbose 3
    bool successful = useDockerContainers ? dockerResamplingOmNifti(niifullpath, outputtniifullpath, isovalue) : false;

    if(!useDockerContainers || !successful) {
        MITK_INFO(!useDockerContainers) << "Using static MIRTK libraries.";
        MITK_WARN(useDockerContainers && (!successful)) << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
        //Absolute path
        QString aPath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
        aPath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
        QDir apathd(aPath);
        if (apathd.exists()) {
            process->setWorkingDirectory(aPath);

            //Setup transformation
            QStringList arguments;
            QString input  = niifullpath;
            QString output = outputtniifullpath;
            QString mirtk  = aPath + mitk::IOUtil::GetDirectorySeparator() + "resample-image";

            arguments << input;
            arguments << output;
            arguments << "-isotropic" << QString::number(isovalue);
            arguments << "-interp" << "CSpline";
            arguments << "-verbose" << "3";

            completion = false;
            process->start(mirtk, arguments);
            while (!completion) {
                std::this_thread::sleep_for(std::chrono::seconds(1));
                QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
            }
            mitk::ProgressBar::GetInstance()->Progress();

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
                         mitk::IOUtil::GetProgramPath();
        }//_if
    }
}

void CemrgCommandLine::ExecuteTransformationOnPoints(
        QString dir, QString meshfullpath, QString outputtmeshfullpath, QString txfullpath) {

    bool successful = useDockerContainers ? dockerTransformationOnPoints(dir, meshfullpath,  outputtmeshfullpath, txfullpath) : false;

    if(!useDockerContainers || !successful) {
        MITK_INFO(!useDockerContainers) << "Using static MIRTK libraries.";
        MITK_WARN(useDockerContainers && (!successful)) << "Docker did not produce a good outcome. Trying with local MIRTK libraries.";
        //Absolute path
        QString aPath = QString::fromStdString(mitk::IOUtil::GetProgramPath()) + mitk::IOUtil::GetDirectorySeparator() + "MLib";
#if defined(__APPLE__)
        aPath = mitk::IOUtil::GetDirectorySeparator() + QString("Applications") +
                mitk::IOUtil::GetDirectorySeparator() + QString("CemrgApp") +
                mitk::IOUtil::GetDirectorySeparator() + QString("MLib");
#endif
        QDir apathd(aPath);
        if (apathd.exists()) {
            process->setWorkingDirectory(aPath);

            //Setup transformation
            QStringList arguments;
            QString input  = meshfullpath;
            QString output = outputtmeshfullpath;
            QString mirtk  = aPath + mitk::IOUtil::GetDirectorySeparator() + "transform-points";

            arguments << input;
            arguments << output;
            arguments << "-dofin" << txfullpath;
            arguments << "-ascii";
            arguments << "-verbose" << "3";

            completion = false;
            process->start(mirtk, arguments);
            while (!completion) {
                std::this_thread::sleep_for(std::chrono::seconds(1));
                QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
            }
            mitk::ProgressBar::GetInstance()->Progress();

        } else {
            QMessageBox::warning(NULL, "Please check the LOG", "MIRTK libraries not found");
            MITK_WARN << "MIRTK libraries not found. Please make sure the MLib folder is inside the directory;\n\t"+
                         mitk::IOUtil::GetProgramPath();
        }//_if
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


/* DOCKER UTILITIES.
The principal workflow of the following routines is:

- Run the specific function needed from the particular image (MIRTK, MESHTOOLS3D, ...)
- Search for the output in the corresponding directory and
- Return a boolean value as to this search was successful.

The boolean value will determine the principal functions to assess whether the
static libraries need to be used.

Normal docker command looks similar to:

> docker run --rm  --volume=/absolute/path/to/directory:/data IMAGE [FUNCTION] [OPTIONS] [FLAGS]

The --volume flag is the absolute path to the working directory and any other
path specified must be relative to the volume specified.
*/

bool CemrgCommandLine::dockerRegistration(
        QString directory, QString fixed, QString moving, QString txname, QString modelname) {

    MITK_INFO << "[ATTENTION] Attempting MIRTK REGISTRATION using Docker.";
    QString aPath = "";
#if defined(__APPLE__)
    aPath = "/usr/local/bin/";
#endif

    QDir mirtkhome(directory);
    QString fixedRelativePath = mirtkhome.relativeFilePath(fixed);
    QString movingRelativePath = mirtkhome.relativeFilePath(moving);
    QString dofRelativePath = mirtkhome.relativeFilePath(txname);

    process->setWorkingDirectory(mirtkhome.absolutePath());

    //Setup docker
    QStringList arguments;
    QString docker = aPath+"docker";
    QString dockerimage = "biomedia/mirtk:v1.1.0";

    arguments << "run" << "--rm";
    arguments << "--volume="+mirtkhome.absolutePath()+":/data";

    //Setup registration
    QString input1 = movingRelativePath;
    QString input2 = fixedRelativePath;
    QString output = dofRelativePath;
    QString dockerexe  = "register";

    arguments << dockerimage;
    arguments << dockerexe;
    arguments << input1;
    arguments << input2;
    arguments << "-dofout" << output;
    arguments << "-model" << modelname;
    arguments << "-verbose" << "3";
    arguments << "-color";

    MITK_INFO << "Executing a " + modelname.toStdString() + " registration.";
    MITK_INFO << ("Performing a " + modelname + " registration").toStdString();
    MITK_INFO << printFullCommand(docker, arguments);
    MITK_INFO << "\n";

    completion = false;
    process->start(docker, arguments);
    checkForStartedProcess();
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();

    QString  outAbsolutepath = mirtkhome.absolutePath() +
            mitk::IOUtil::GetDirectorySeparator() + dofRelativePath;
    bool successful = isOutputSuccessful(outAbsolutepath);
    if (!successful)
        MITK_WARN << "Docker unsuccessful. Check your configuration.";

    return successful;
}

bool CemrgCommandLine::dockerTranformation(
        QString directory, QString imgNamefullpath, QString regImgNamefullpath, QString txfullpath) {

    MITK_INFO << "[ATTENTION] Attempting MIRTK IMAGE TRANSFORMATION using Docker.";
    QString aPath = "";
#if defined(__APPLE__)
    aPath = "/usr/local/bin/";
#endif

    QDir mirtkhome(directory);
    QString inputRelativePath = mirtkhome.relativeFilePath(imgNamefullpath);
    QString outputRelativePath = mirtkhome.relativeFilePath(regImgNamefullpath);
    QString dofRelativePath = mirtkhome.relativeFilePath(txfullpath);

    process->setWorkingDirectory(mirtkhome.absolutePath());

    //Setup docker
    QStringList arguments;
    QString docker = aPath+"docker";
    QString dockerimage = "biomedia/mirtk:v1.1.0";

    arguments << "run" << "--rm";
    arguments << "--volume="+mirtkhome.absolutePath()+":/data";

    //Setup transformation
    QString dockerexe  = "transform-image";

    arguments << dockerimage;
    arguments << dockerexe;
    arguments << inputRelativePath;
    arguments << outputRelativePath;
    arguments << "-dofin" << dofRelativePath;
    arguments << "-verbose" << "3";
    MITK_INFO << printFullCommand(docker, arguments);

    completion = false;
    process->start(docker, arguments);
    checkForStartedProcess();
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();

    QString  outAbsolutepath = mirtkhome.absolutePath() +
            mitk::IOUtil::GetDirectorySeparator() + outputRelativePath;

    bool successful = isOutputSuccessful(outAbsolutepath);
    if (!successful)
        MITK_WARN << "Docker unsuccessful. Check your configuration.";

    return successful;
}

bool CemrgCommandLine::dockerTransformationOnPoints(
        QString directory, QString meshfullpath, QString outputtmeshfullpath, QString txfullpath) {

    MITK_INFO << "[ATTENTION] Attempting MIRTK MESH TRANSFORMATION using Docker.";
    QString aPath = "";
#if defined(__APPLE__)
    aPath = "/usr/local/bin/";
#endif

    QDir mirtkhome(directory);
    QString inputRelativePath = mirtkhome.relativeFilePath(meshfullpath);
    QString outputRelativePath = mirtkhome.relativeFilePath(outputtmeshfullpath);
    QString dofRelativePath = mirtkhome.relativeFilePath(txfullpath);

    process->setWorkingDirectory(mirtkhome.absolutePath());

    //Setup docker
    QStringList arguments;
    QString docker = aPath+"docker";
    QString dockerimage = "biomedia/mirtk:v1.1.0";

    arguments << "run" << "--rm";
    arguments << "--volume="+mirtkhome.absolutePath()+":/data";

    //Setup transformation
    QString dockerexe  = "transform-points";

    arguments << dockerimage;
    arguments << dockerexe;
    arguments << inputRelativePath;
    arguments << outputRelativePath;
    arguments << "-dofin" << dofRelativePath;
    arguments << "-nocompress";
    arguments << "-ascii";
    arguments << "-verbose" << "3";
    MITK_INFO << printFullCommand(docker, arguments);

    completion = false;
    process->start(docker, arguments);
    checkForStartedProcess();
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();

    QString  outAbsolutepath = mirtkhome.absolutePath() +
            mitk::IOUtil::GetDirectorySeparator() + outputRelativePath;

    bool successful = isOutputSuccessful(outAbsolutepath);
    if (!successful)
        MITK_WARN << "Docker unsuccessful. Check your configuration.";

    return successful;
}

QString CemrgCommandLine::dockerExpandSurf(
        QString dir, QString segPath, int iter, float th, int blur, int smth) {

    MITK_INFO << "[ATTENTION] Attempting SURFACE CREATION using Docker.";
    QString aPath = "";
#if defined(__APPLE__)
    aPath = "/usr/local/bin/";
#endif

    QDir mirtkhome(dir);
    QString inputRelativePath = mirtkhome.relativeFilePath(segPath);
    QString outputRelativePath = mirtkhome.relativeFilePath("segmentation.s.nii");

    process->setWorkingDirectory(mirtkhome.absolutePath());

    //Setup docker
    QString docker = aPath+"docker";
    QString dockerimage = "biomedia/mirtk:v1.1.0";
    QStringList arguments;

    arguments << "run" << "--rm";
    arguments << "--volume="+mirtkhome.absolutePath()+":/data";

    QString dockerexe  = "dilate-image"; //Dilation followed by Erosion

    arguments << dockerimage;
    arguments << dockerexe;
    arguments << inputRelativePath;
    arguments << outputRelativePath;
    arguments << "-iterations" << QString::number(iter);
    MITK_INFO << printFullCommand(docker, arguments);

    completion = false;
    process->start(docker, arguments);
    checkForStartedProcess();
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();

    //Reset arguments
    arguments.clear();
    arguments << "run" << "--rm";
    arguments << "--volume="+mirtkhome.absolutePath()+":/data";
    arguments << dockerimage;

    //Marching Cubes
    inputRelativePath  = outputRelativePath;
    outputRelativePath = mirtkhome.relativeFilePath("segmentation.vtk");
    dockerexe = "extract-surface";

    arguments << dockerexe;
    arguments << inputRelativePath;
    arguments << outputRelativePath;
    arguments << "-isovalue" << QString::number(th);
    arguments << "-blur" << QString::number(blur);
    arguments << "-ascii";
    arguments << "-verbose" << "3";
    completion = false;
    MITK_INFO << printFullCommand(docker, arguments);

    process->start(docker, arguments);
    checkForStartedProcess();
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();

    //Reset arguments
    arguments.clear();
    arguments << "run" << "--rm";
    arguments << "--volume="+mirtkhome.absolutePath()+":/data";
    arguments << dockerimage;

    //Smoothing
    inputRelativePath  = outputRelativePath;
    outputRelativePath = outputRelativePath+"";
    dockerexe = "smooth-surface";

    arguments << dockerexe;
    arguments << inputRelativePath;
    arguments << outputRelativePath;
    arguments << "-iterations" << QString::number(smth);
    arguments << "-verbose" << "3";
    completion = false;
    MITK_INFO << printFullCommand(docker, arguments);

    process->start(docker, arguments);
    checkForStartedProcess();
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();

    //Return path to output mesh
    remove((dir + mitk::IOUtil::GetDirectorySeparator() + "segmentation.s.nii").toStdString().c_str());
    QString  outAbsolutepath = mirtkhome.absolutePath() +
            mitk::IOUtil::GetDirectorySeparator() + outputRelativePath;
    return outAbsolutepath;
}

QString CemrgCommandLine::dockerSurf(
        QString dir, QString segPath, int iter, float th, int blur, int smth) {

    MITK_INFO << "[ATTENTION] Attempting SURFACE CREATION using Docker.";
    QString aPath = "";
#if defined(__APPLE__)
    aPath = "/usr/local/bin/";
#endif
    QDir mirtkhome(dir);

    QString inputRelativePath = mirtkhome.relativeFilePath(segPath);
    QString outputRelativePath = mirtkhome.relativeFilePath("segmentation.s.nii");

    process->setWorkingDirectory(mirtkhome.absolutePath());

    //Setup docker
    QString docker = aPath+"docker";
    QString dockerimage = "biomedia/mirtk:v1.1.0";
    QStringList arguments;

    arguments << "run" << "--rm";
    arguments << "--volume="+mirtkhome.absolutePath()+":/data";

    QString dockerexe  = "close-image"; //Dilation followed by Erosion

    arguments << dockerimage;
    arguments << dockerexe;
    arguments << inputRelativePath;
    arguments << outputRelativePath;
    arguments << "-iterations" << QString::number(iter);
    arguments << "-verbose" << "3";

    MITK_INFO << printFullCommand(docker, arguments);

    completion = false;
    process->start(docker, arguments);
    checkForStartedProcess();
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();

    //Reset arguments
    arguments.clear();
    arguments << "run" << "--rm";
    arguments << "--volume="+mirtkhome.absolutePath()+":/data";
    arguments << dockerimage;

    //Marching Cubes
    inputRelativePath  = outputRelativePath;
    outputRelativePath = mirtkhome.relativeFilePath("segmentation.vtk");
    dockerexe = "extract-surface";

    arguments << dockerexe;
    arguments << inputRelativePath;
    arguments << outputRelativePath;
    arguments << "-isovalue" << QString::number(th);
    arguments << "-blur" << QString::number(blur);
    arguments << "-ascii";
    arguments << "-verbose" << "3";
    MITK_INFO << printFullCommand(docker, arguments);

    completion = false;
    process->start(docker, arguments);
    checkForStartedProcess();
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();

    //Reset arguments
    arguments.clear();
    arguments << "run" << "--rm";
    arguments << "--volume="+mirtkhome.absolutePath()+":/data";
    arguments << dockerimage;

    //Smoothing
    inputRelativePath  = outputRelativePath;
    outputRelativePath = outputRelativePath+"";
    dockerexe = "smooth-surface";

    arguments << dockerexe;
    arguments << inputRelativePath;
    arguments << outputRelativePath;
    arguments << "-iterations" << QString::number(smth);
    arguments << "-verbose" << "3";
    MITK_INFO << printFullCommand(docker, arguments);

    completion = false;
    process->start(docker, arguments);
    checkForStartedProcess();
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();

    //Return path to output mesh
    remove((dir + mitk::IOUtil::GetDirectorySeparator() + "segmentation.s.nii").toStdString().c_str());
    QString  outAbsolutepath = mirtkhome.absolutePath() +
            mitk::IOUtil::GetDirectorySeparator() + outputRelativePath;
    return outAbsolutepath;
}

bool CemrgCommandLine::dockerResamplingOmNifti(
        QString niifullpath, QString outputtniifullpath, int isovalue) {

    // /resample-image niifullpath outputtniifullpath -isotropic 0.5 -interp CSpline -verbose 3
    MITK_INFO << "[ATTENTION] Attempting RESAMPLING IMAGE using Docker.";
    QString aPath = "";
#if defined(__APPLE__)
    aPath = "/usr/local/bin/";
#endif

    QFileInfo inputnii(niifullpath);
    QDir mirtkhome(inputnii.absolutePath());
    QString inputRelativePath = mirtkhome.relativeFilePath(niifullpath);
    QString outputRelativePath = mirtkhome.relativeFilePath(outputtniifullpath);

    process->setWorkingDirectory(mirtkhome.absolutePath());

    //Setup docker
    QString docker = aPath+"docker";
    QString dockerimage = "biomedia/mirtk:v1.1.0";
    QStringList arguments;

    arguments << "run" << "--rm";
    arguments << "--volume="+mirtkhome.absolutePath()+":/data";
    arguments << dockerimage;

    QString dockerexe  = "resample-image"; //Dilation followed by Erosion

    arguments << dockerexe;
    arguments << inputRelativePath;
    arguments << outputRelativePath;
    arguments << "-isotropic" << QString::number(isovalue);
    arguments << "-interp" << "CSpline";
    arguments << "-verbose" << "3";

    MITK_INFO << printFullCommand(docker, arguments);

    completion = false;
    process->start(docker, arguments);
    checkForStartedProcess();
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();

    bool successful = isOutputSuccessful(outputtniifullpath);
    if (!successful)
        MITK_WARN << "Docker unsuccessful. Check your configuration.";
    return successful;
}

//Tracking Utilities - Docker
bool CemrgCommandLine::dockerTracking(QString dir, QString imgTimes, QString param) {

    MITK_INFO << "[ATTENTION] Attempting TRACKING (registration) using Docker";
    QString aPath = "";
#if defined(__APPLE__)
    aPath = "/usr/local/bin/";
#endif

    QDir mirtkhome(dir);
    QString inputRelativePath = mirtkhome.relativeFilePath(imgTimes);
    QString outputRelativePath = mirtkhome.relativeFilePath("tsffd.dof");

    process->setWorkingDirectory(mirtkhome.absolutePath());

    //Setup docker
    QString docker = aPath+"docker";
    QString dockerimage = "biomedia/mirtk:v1.1.0";
    QString dockerexe  = "register";
    QStringList arguments;

    arguments << "run" << "--rm";
    arguments << "--volume="+mirtkhome.absolutePath()+":/data";
    arguments << dockerimage;
    arguments << dockerexe;
    arguments << "-images" << inputRelativePath;

    if (!param.isEmpty()) {
        QString paramRelativePath = mirtkhome.relativeFilePath(param);
        arguments << "-parin" << paramRelativePath;
    }
    arguments << "-dofout" << outputRelativePath;
    arguments << "-verbose" << "3";

    completion = false;

    MITK_INFO << printFullCommand(docker, arguments);

    process->start(docker, arguments);
    checkForStartedProcess();
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();
    QString  outAbsolutepath = mirtkhome.absolutePath() +
            mitk::IOUtil::GetDirectorySeparator() + outputRelativePath;

    bool successful = isOutputSuccessful(outAbsolutepath);
    MITK_WARN(!successful) << "Docker unsuccessful. Check your configuration.";
    return successful;
}

bool CemrgCommandLine::dockerApplying(
        QString dir, QString inputMesh, double iniTime, QString dofin, int noFrames, int smooth) {

    MITK_INFO << "[ATTENTION] Attempting APPLYING (registration) using Docker";
    QString aPath = "";
#if defined(__APPLE__)
    aPath = "/usr/local/bin/";
#endif

    QDir mirtkhome(dir);
    QString inputRelativePath = mirtkhome.relativeFilePath(inputMesh);
    QString outputRelativePath_ = mirtkhome.relativeFilePath("transformed-");
    QString dofRelativePath = mirtkhome.relativeFilePath(dofin);

    process->setWorkingDirectory(mirtkhome.absolutePath());

    //Setup docker
    QString docker = aPath+"docker";
    QString dockerimage = "biomedia/mirtk:v1.1.0";
    QString dockerexe  = "transform-points";
    QStringList arguments;

    arguments << "run" << "--rm";
    arguments << "--volume="+mirtkhome.absolutePath()+":/data";
    arguments << dockerimage;
    arguments << dockerexe;

    int fctTime = 10;
    noFrames *= smooth;
    int suxs = 0;
    if (smooth == 2)
        fctTime = 5;
    else if (smooth == 5)
        fctTime = 2;

    QString  outAbsolutepath = mirtkhome.absolutePath() +
            mitk::IOUtil::GetDirectorySeparator() + outputRelativePath_;
    for (int i=0; i<noFrames; i++) {

        arguments.clear();
        arguments << "run" << "--rm";
        arguments << "--volume="+mirtkhome.absolutePath()+":/data";
        arguments << dockerimage;
        arguments << dockerexe;
        arguments << inputRelativePath;
        arguments << outputRelativePath_ + QString::number(i) + ".vtk";
        arguments << "-dofin" << dofRelativePath;
        arguments << "-ascii";
        arguments << "-St";
        arguments << QString::number(iniTime);
        arguments << "-verbose" << "3";

        MITK_INFO << printFullCommand(docker, arguments);

        completion = false;
        process->start(docker, arguments);
        checkForStartedProcess();
        while (!completion) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
        }

        suxs += (isOutputSuccessful(outAbsolutepath + QString::number(i) + ".vtk")) ? 1 : 0;
        iniTime += fctTime;
        mitk::ProgressBar::GetInstance()->Progress();
    }

    bool successful = (suxs == noFrames);
    return successful;
}

bool CemrgCommandLine::dockerSimpleTranslation(
        QString dir, QString sourceMeshP, QString targetMeshP, QString outputPath, bool transformThePoints) {

    MITK_INFO << "[ATTENTION] Attempting INIT-DOF + TRANSFORM-POINTS using Docker.";
    QString aPath = "";
#if defined(__APPLE__)
    aPath = "/usr/local/bin/";
#endif

    QDir mirtkhome(dir);
    QString sourceRelativePath = mirtkhome.relativeFilePath(sourceMeshP);
    QString targetRelativePath = mirtkhome.relativeFilePath(targetMeshP);
    QString txRelativePath = mirtkhome.relativeFilePath(".init-tx.dof");
    QString outputRelativePath = mirtkhome.relativeFilePath(outputPath);

    process->setWorkingDirectory(mirtkhome.absolutePath());

    //Setup docker
    QString docker = aPath+"docker";
    QString dockerimage = "biomedia/mirtk:v1.1.0";
    QStringList arguments;

    arguments << "run" << "--rm";
    arguments << "--volume="+mirtkhome.absolutePath()+":/data";

    QString dockerexe  = "init-dof"; //simple translation

    arguments << dockerimage;
    arguments << dockerexe;
    arguments << txRelativePath;
    arguments << "-translations" << "-norotations" << "-noscaling" << "-noshearing";
    if(transformThePoints) {
        arguments << "-displacements";
        arguments << sourceRelativePath;
        arguments << targetRelativePath;
    }
    else{
        arguments << "-source" << sourceRelativePath;
        arguments << "-target" << targetRelativePath;
    }

    arguments << "-verbose" << "3";

    MITK_INFO << printFullCommand(docker, arguments);

    completion = false;
    process->start(docker, arguments);
    checkForStartedProcess();
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();

    QString outAbsolutepath;

    if (transformThePoints) {
        MITK_INFO << "Performing the transformation of the points";
        //Reset arguments
        arguments.clear();
        arguments << "run" << "--rm";
        arguments << "--volume="+mirtkhome.absolutePath()+":/data";
        arguments << dockerimage;

        //Transformation
        dockerexe = "transform-points";
        arguments << dockerexe;
        arguments << sourceRelativePath;
        arguments << outputRelativePath;
        arguments << "-dofin" << txRelativePath;
        arguments << "-ascii";
        arguments << "-verbose" << "3";

        MITK_INFO << printFullCommand(docker, arguments);

        completion = false;
        process->start(docker, arguments);
        checkForStartedProcess();
        while (!completion) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
        }
        mitk::ProgressBar::GetInstance()->Progress();

        outAbsolutepath = mirtkhome.absolutePath() +
                mitk::IOUtil::GetDirectorySeparator() + outputRelativePath;

    } else {
        MITK_INFO << "Image detected, not performing transformation";
        outAbsolutepath = mirtkhome.absolutePath() +
                mitk::IOUtil::GetDirectorySeparator() + txRelativePath;
    }

    bool successful = isOutputSuccessful(outAbsolutepath);
    if (!successful)
        MITK_WARN << "Docker unsuccessful. Check your configuration.";
    return successful;
}

//Docker - meshtools3d
QString CemrgCommandLine::dockerCreateCGALMesh(QString dir, QString fileName, QString templatePath) {

    MITK_INFO << "[ATTENTION] Attempting CreateCGALMesh (Meshtools3D) using Docker.";
    QString aPath = "";
#if defined(__APPLE__)
    aPath = "/usr/local/bin/";
#endif

    QDir meshtools3dhome(dir);
    QString templateRelativePath = meshtools3dhome.relativeFilePath(templatePath);
    QString segRelativePath = meshtools3dhome.relativeFilePath(dir);
    //QString outnameRelativePath = fileName;
    QString outdirRelativePath = meshtools3dhome.relativeFilePath("CGALMeshDir");

    process->setWorkingDirectory(meshtools3dhome.absolutePath());

    //Setup docker
    QString docker = aPath+"docker";
    QString dockerimage = "alonsojasl/meshtools3d:v1.0";
    QStringList arguments;

    arguments << "run" << "--rm";
    arguments << "--volume="+meshtools3dhome.absolutePath()+":/data";
    arguments << dockerimage;
    arguments << "-f" << templateRelativePath;
    arguments << "-seg_dir" << segRelativePath;
    arguments << "-seg_name" << "converted.inr";
    arguments << "-out_dir" << outdirRelativePath;
    arguments << "-out_name" << fileName;

    MITK_INFO << printFullCommand(docker, arguments);

    completion = false;
    process->start(docker, arguments);
    checkForStartedProcess();
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();

    QString  outAbsolutepath = meshtools3dhome.absolutePath() +
            mitk::IOUtil::GetDirectorySeparator() + "CGALMeshDir" +
            mitk::IOUtil::GetDirectorySeparator() + fileName + ".vtk";

    bool successful = isOutputSuccessful(outAbsolutepath);
    MITK_WARN((!successful)) << "Docker unsuccessful. Check your configuration.";
    return outAbsolutepath;
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
bool CemrgCommandLine::isOutputSuccessful(QString outputfullpath) {

    MITK_INFO << "[ATTENTION] Checking for successful output on path:";
    MITK_INFO << outputfullpath.toStdString();
    QFileInfo finfo(outputfullpath);
    bool res = finfo.exists();
    MITK_INFO << (res ? "Successful output" : "Output file not found.");
    return res;
}

void CemrgCommandLine::ExecuteTouch(QString filepath) {
    QString commandName;
    QStringList arguments;
    #ifdef _WIN32
    commandName = "echo"; // echo . > filepath
    arguments << "." << ">";
    #else
    commandName = "touch"; // touch filepath
    #endif
    arguments << filepath;
    completion = false;

    MITK_INFO << printFullCommand(commandName, arguments);
    process->start(commandName, arguments);
    checkForStartedProcess();
    while (!completion) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
    }
    mitk::ProgressBar::GetInstance()->Progress();
}

std::string CemrgCommandLine::printFullCommand(QString command, QStringList arguments) {

    QString teststr = "";
    for (int ix=0; ix < arguments.size(); ix++) {
        teststr += arguments.at(ix) + " ";
    }

    bool debugging = true;
    if (debugging) {
        QString prodPath = QString::fromStdString(mitk::IOUtil::GetProgramPath());
        MITK_INFO << ("Program path: " + prodPath).toStdString();
        ofstream prodFile1;
        prodFile1.open((prodPath + "dockerDebug.txt").toStdString(), ofstream::out | ofstream::app);
        prodFile1 << (command + " " + teststr).toStdString() << "\n";
        prodFile1.close();
    }//_if
    return (command + " " + teststr).toStdString();
}

/* CHECK FOR STARTED PROCESS
this function prevents freezing of the app when something goes wrong with the
Qt process.
*/
void CemrgCommandLine::checkForStartedProcess() {
    bool debugvar = false;
    if(debugvar) {
        QStringList errinfo = QProcess::systemEnvironment();
        QString teststr = "";
        for (int ix=0; ix < errinfo.size(); ix++) {
            teststr += errinfo.at(ix) + " ";
        }
        MITK_INFO << "SYSTEM ENVIRONMENT:";
        MITK_INFO << teststr.toStdString();
    }

    if(process->waitForStarted()) {
        MITK_INFO << "Starting process";
    } else {
        completion=true;
        MITK_WARN << "[ATTENTION] Process error!";
        MITK_INFO << "STATE:";
        MITK_INFO << process->state();
        MITK_INFO << "ERROR:";
        MITK_INFO << process->error();
    }
}

void CemrgCommandLine::setUseDockerContainers(bool dockerContainersOnOff) {

    QString onoff = dockerContainersOnOff ? "ON" : "OFF";
    MITK_INFO << ("[...] Setting useDockerContainers variable to: " + onoff).toStdString();
    useDockerContainers = dockerContainersOnOff;
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
