/*=========================================================================

Program:   Medical Imaging & Interaction Toolkit
Module:    $RCSfile$
Language:  C++
Date:      $Date$
Version:   $Revision: 13820 $

Copyright (c) German Cancer Research Center, Division of Medical and
Biological Informatics. All rights reserved.
See MITKCopyright.txt or http://www.mitk.org/ for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include <mitkBaseApplication.h>

#include <QStringList>
#include <QString>
#include <QVariant>
#include <QPixmap>
#include <QTimer>
#include <QSplashScreen>
#include <QFile>
#include <QTextStream>

int main(int argc, char** argv) {

    QString version, sha;
    QFile file(":/splash/version.txt");
    if (file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QTextStream in(&file);
        version = in.readLine();
        sha = in.readLine();
    }

    mitk::BaseApplication myApp(argc, argv);
    myApp.setSingleMode(true);
    myApp.setApplicationName("CemrgApp " + version);
    myApp.setOrganizationName("KCL");
    myApp.initializeQt();

    QPixmap pixmap(":/splash/splashscreen.png");
    QSplashScreen splash(pixmap);
    splash.setWindowFlags(Qt::SplashScreen | Qt::WindowStaysOnTopHint | Qt::FramelessWindowHint);
    splash.show();
    if (!version.isEmpty())
        splash.showMessage("CemrgApp " + version + " (" + sha + ")\nPowered by: MITK v2018.04.2", Qt::AlignLeft, Qt::white);
    QTimer::singleShot(4000, &splash, SLOT(close()));

    // -------------------------------------------------------------------
    // Here you can switch to your customizable application:
    // -------------------------------------------------------------------

    QStringList preloadLibs;
    preloadLibs << "liborg_mitk_gui_qt_common";
    myApp.setPreloadLibraries(preloadLibs);
    //myApp.setProperty(mitk::BaseApplication::PROP_APPLICATION, "org.mitk.qt.extapplication");
    myApp.setProperty(mitk::BaseApplication::PROP_APPLICATION, "kcl.cemrgapp.mainapp");
    return myApp.run();
}