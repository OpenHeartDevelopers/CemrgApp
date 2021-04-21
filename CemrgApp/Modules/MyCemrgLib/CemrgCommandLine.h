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

#ifndef CemrgCommandLine_h
#define CemrgCommandLine_h

// Qt
#include <memory>
#include <QProcess>
#include <QTextEdit>
#include <QVBoxLayout>
#include <MyCemrgLibExports.h>


class MyCemrgLib_EXPORT CemrgCommandLine : public QObject {

    Q_OBJECT
    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)

public:

    CemrgCommandLine();
    ~CemrgCommandLine();
    QDialog* GetDialog();

    QString ExecuteSurf(QString dir, QString segPath, int iter, float th, int blur, int smth);
    QString ExecuteCreateCGALMesh(QString dir, QString fileName, QString templatePath);

    void ExecuteTracking(QString dir, QString imgTimes, QString param);
    void ExecuteApplying(QString dir, QString inputMesh, double iniTime, QString dofin, int noFrames, int smooth);
    void ExecuteRegistration(QString dir, QString lge, QString mra);
    void ExecuteTransformation(QString dir, QString imgName, QString regImgName);

    bool ConnectToServer(QString userID, QString server);
    bool TransferTFServer(QString directory, QString fname, QString userID, QString server, bool download);
    void GPUReconstruction(QString userID, QString server, QStringList imgsList, QString targetImg, double resolution, double delta, int package, QString out);

protected slots:

    void UpdateStdText();
    void UpdateErrText();
    void FinishedAlert();

protected:


private:

    //QProcess, dial and panel
    QDialog* dial;
    QTextEdit* panel;
    QVBoxLayout* layout;
    std::unique_ptr<QProcess> process;
    bool completion;
};

#endif // CemrgCommandLine_h
