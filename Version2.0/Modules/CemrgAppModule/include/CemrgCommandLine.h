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

#ifndef CemrgCommandLine_h
#define CemrgCommandLine_h

// Qt
#include <memory>
#include <QProcess>
#include <QTextEdit>
#include <QVBoxLayout>
#include <MitkCemrgAppModuleExports.h>

class MITKCEMRGAPPMODULE_EXPORT CemrgCommandLine : public QObject {

    Q_OBJECT
    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)

public:

    CemrgCommandLine();
    CemrgCommandLine(bool cmd);
    ~CemrgCommandLine();
    QDialog* GetDialog();

    QString ExecuteSurf(QString dir, QString segPath, int iter, float th, int blur, int smth);
    QString ExecuteCreateCGALMesh(QString dir, QString fileName, QString templatePath);
    void ExecuteTracking(QString dir, QString imgTimes, QString param);
    void ExecuteApplying(QString dir, QString inputMesh, double iniTime, QString dofin, int noFrames, int smooth);
    void ExecuteRegistration(QString dir, QString fixed, QString moving, QString txname="rigid.dof", QString modelname="Rigid");
    void ExecuteTransformation(QString dir, QString imgNamefullpath, QString regImgNamefullpath, QString txfullpath="rigid.dof");
    void ExecuteResamplingOnNifti(QString niifullpath, QString outputtniifullpath, int isovalue);
    void ExecuteTransformationOnPoints(QString dir, QString meshfullpath, QString outputtmeshfullpath, QString txfullpath);

    bool ConnectToServer(QString userID, QString server);
    bool TransferTFServer(QString directory, QString fname, QString userID, QString server, bool download);
    void GPUReconstruction(QString userID, QString server, QStringList imgsList, QString targetImg, double resolution, double delta, int package, QString out);

    //Docker
    bool dockerRegistration(QString directory, QString fixed, QString moving, QString txname, QString modelname);
    bool dockerTranformation(QString directory, QString imgNamefullpath, QString regImgNamefullpath, QString txfullpath);
    bool dockerTransformationOnPoints(QString directory, QString meshfullpath, QString outputtmeshfullpath, QString txfullpath);
    QString dockerSurf(QString dir, QString segPath, int iter, float th, int blur, int smth);
    QString dockerExpandSurf(QString dir, QString segPath, int iter, float th, int blur, int smth);
    bool dockerResamplingOmNifti(QString niifullpath, QString outputtniifullpath, int isovalue);
    bool dockerTracking(QString dir, QString imgTimes, QString param);
    bool dockerApplying(QString dir, QString inputMesh, double iniTime, QString dofin, int noFrames, int smooth);
    bool dockerSimpleTranslation(QString dir, QString sourceMeshP, QString targetMeshP, QString outputPath, bool transformThePoints = true);
    QString dockerCreateCGALMesh(QString dir, QString fileName, QString templatePath);
    QString dockerCemrgNetPrediction(QString mra);

    // Helper functions
    bool isOutputSuccessful(QString outputfullpath);
    void ExecuteTouch(QString filepath);
    std::string printFullCommand(QString command, QStringList arguments);
    void checkForStartedProcess();
    void setUseMIRKTDocker(bool dockerOnMirtk);

protected slots:

    void UpdateStdText();
    void UpdateErrText();
    void FinishedAlert();

private:

    //QProcess, dial and panel
    QDialog* dial;
    QTextEdit* panel;
    QVBoxLayout* layout;
    std::unique_ptr<QProcess> process;
    bool completion;
    bool isUI;
    bool useMirtkDocker;
};

#endif // CemrgCommandLine_h
