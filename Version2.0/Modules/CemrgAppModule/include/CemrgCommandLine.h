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
    CemrgCommandLine(std::string dockerimage);
    ~CemrgCommandLine();
    QDialog* GetDialog();

    // Execute functions
    QString ExecuteSurf(QString dir, QString segPath, int iter, float th, int blur, int smth);
    QString ExecuteCreateCGALMesh(QString dir, QString outputName, QString paramsFullPath, QString segmentationName="converted.inr");
    void ExecuteTracking(QString dir, QString imgTimes, QString param, QString output="tsffd.dof");
    void ExecuteApplying(QString dir, QString inputMesh, double iniTime, QString dofin, int noFrames, int smooth);
    void ExecuteRegistration(QString dir, QString fixed, QString moving, QString txname="rigid.dof", QString modelname="Rigid");
    void ExecuteTransformation(QString dir, QString imgNamefullpath, QString regImgNamefullpath, QString txfullpath="rigid.dof");
    void ExecuteResamplingOnNifti(QString niifullpath, QString outputtniifullpath, int isovalue);
    void ExecuteTransformationOnPoints(QString dir, QString meshfullpath, QString outputtmeshfullpath, QString txfullpath, double applyingIniTime=-100);
    QString ExecuteExpandSurf(QString dir, QString segPath, int iter=1, float th=0.5, int blur=0, int smth=10);
    void ExecuteSimpleTranslation(QString dir, QString sourceMeshP, QString targetMeshP, QString txName="init.dof", bool transformThePoints=true);

    QString ExecuteMorphologicalOperation(QString operation, QString dir, QString segPath, QString outputPath = "segmentation.s.nii", int iter=1);
    QString ExecuteExtractSurface(QString dir, QString segPath, QString outputPath = "segmentation.vtk", float th=0.5, int blur=0);
    QString ExecuteSmoothSurface(QString dir, QString segPath, QString outputPath, int smth=10);

    // Server functions
    bool ConnectToServer(QString userID, QString server);
    bool TransferTFServer(QString directory, QString fname, QString userID, QString server, bool download);
    void GPUReconstruction(QString userID, QString server, QStringList imgsList, QString targetImg, double resolution, double delta, int package, QString out);

    // Docker specific functions
    QString dockerCemrgNetPrediction(QString mra);

    // Helper functions
    bool ExecuteCommand(QString executableName, QStringList arguments, QString outputPath);
    void ExecuteTouch(QString filepath);
    bool isOutputSuccessful(QString outputfullpath);
    std::string printFullCommand(QString command, QStringList arguments);
    bool checkForStartedProcess();

    // Docker helper functions
    QStringList getDockerArguments(QString volume, QString dockerexe = "");
    void setUseDockerContainers(bool dockerContainersOnOff);

    inline void setUseDockerContainersOn(){setUseDockerContainers(true);};
    inline void setUseDockerContainersOff(){setUseDockerContainers(false);};
    inline void setDockerImage(QString dockerimage){_dockerimage = dockerimage;};
    inline void setDockerImage(std::string dockerimage){_dockerimage = QString::fromStdString(dockerimage);};
    inline QString getDockerImage(){return _dockerimage;};

protected slots:

    void UpdateStdText();
    void UpdateErrText();
    void FinishedAlert();

private:

    //QProcess, dial and panel
    QDialog* dial;
    QTextEdit* panel;
    QVBoxLayout* layout;
    QString _dockerimage;
    std::unique_ptr<QProcess> process;
    bool completion;
    bool _useDockerContainers;
};

#endif // CemrgCommandLine_h
