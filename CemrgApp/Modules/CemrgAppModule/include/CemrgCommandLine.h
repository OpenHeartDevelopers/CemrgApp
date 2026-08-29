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

class MITKCEMRGAPPMODULE_EXPORT CemrgCommandLine: public QObject {

    Q_OBJECT
        // this is needed for all Qt objects that should have a Qt meta-object
        // (everything that derives from QObject and wants to have signal/slots)

public:

    CemrgCommandLine();
    ~CemrgCommandLine();
    QDialog* GetDialog();

    //Execute Plugin Specific Functions
    QString ExecuteSurf(QString dir, QString segPath, float th = 0.5, int blur = 0, int smth = 10);
    QString ExecuteCreateCGALMesh(QString dir, QString outputName, QString paramsFullPath, QString segmentationName = "converted.inr");
    void ExecuteTracking(QString dir, QString imgTimes, QString param, QString output = "tsffd.dof");
    void ExecuteApplying(QString dir, QString inputMesh, double iniTime, QString dofin, int noFrames, int smooth);
    void ExecuteRegistration(QString dir, QString fixed, QString moving, QString transformFileName = "rigid.dof", QString modelname = "Rigid");
    void ExecuteTransformation(QString dir, QString imgNamefullpath, QString regImgNamefullpath, QString transformFileFullPath = "rigid.dof");
    void ExecuteSimpleTranslation(QString dir, QString sourceMeshP, QString targetMeshP, QString transformFileName = "init.dof", bool transformThePoints = true);

    //Execute MIRTK Specific Functions
    void ExecuteTransformationOnPoints(QString dir, QString meshFullPath, QString outputMeshFullPath, QString transformFileFullPath, double applyingIniTime = -100);
    void ExecuteResamplingOnNifti(QString niiFullPath, QString outputNiiFullPath, int isovalue);

    //Execute Docker Specific Functions
    QString DockerCemrgNetPrediction(QString mra);
    QString DockerDicom2Nifti(QString path2dicomfolder);
    QString DockerUniversalAtrialCoordinates(QString dir, QString uaccmd, QStringList fibreAtlas, QString meshname, QStringList cmdargs, QStringList landmarks, QString outnameext="");
    QString DockerUacFibreMapping(QString dir, QString uaccmd, QStringList fibreAtlas, QString meshname, QStringList cmdargs, QString landmarks, QString outnameext);
    QString DockerUacMainMode(QString dir, QString stage, QString atrium, QString layer, QString fibre, QString meshname, QStringList tags, QStringList landmarks, bool fourch=false, bool noraa=false, int scale=1000);
    QString DockerUacFibreMappingMode(QString dir, QString atrium, QString layer, QString fibre, QString meshname, bool msh_endo_epi, QString output, bool fourch=false, QString tags="11", QString biproj="100");

    // openCARP docker
    QString OpenCarpDockerLaplaceSolves(QString dir, QString meshName, QString outName, QStringList zeroNames, QStringList oneNames, QStringList regionLabels);
    QString OpenCarpDocker(QString dir, QString paramfile, QString simID);

    // meshtool routines
    QString DockerSurfaceFromMesh(QString dir, QString meshname, QString outname, QString op, QString outputSuffix); // extract surface
    QString DockerExtractGradient(QString dir, QString meshname, QString idatName, QString odatName, bool elemGrad=true); // extract gradient
    QString DockerRemeshSurface(QString dir, QString meshname, QString outname, double hmax=1, double hmin=0.98, double havg=0.3, double surfCorr=0.95); // resample surfmesh
    QString DockerConvertMeshFormat(QString dir, QString imsh, QString ifmt, QString omsh, QString ofmt, double scale=-1);
    QString DockerInterpolateData(QString dir, QString meshname, QString outmesh, QString idatExt, QString odatExt, QString dataType);
    void DockerCleanMeshQuality(QString dir, QString meshname, QString outMesh, double qualityThres, QString ifmt="vtk", QString ofmt="vtk_polydata");

    inline QString DockerInterpolatePoint(QString dir, QString meshname, QString outmesh, QString idatExt, QString odatExt){return DockerInterpolateData(dir, meshname, outmesh, idatExt, odatExt,"nodedata");}; // meshtool interpolate nodedata
    inline QString DockerInterpolateCell(QString dir, QString meshname, QString outmesh, QString idatExt, QString odatExt){return DockerInterpolateData(dir, meshname, outmesh, idatExt, odatExt,"elemdata");}; // meshtool interpolate elemdata
    inline QString DockerInterpolateCloudPoint(QString dir, QString ptsname, QString outmesh, QString idatExt, QString odatExt){return DockerInterpolateData(dir, ptsname, outmesh, idatExt, odatExt,"clouddata");}; // meshtool interpolate clouddata

    inline void SetDebug(bool b) { _debugvar = b; };
    inline void SetDebugOn() { SetDebug(true); };
    inline void SetDebugOff() { SetDebug(false); };

    //Docker Helper Functions
    void SetUseDockerContainers(bool dockerContainersOnOff);
    inline void SetUseDockerContainersOn() {SetUseDockerContainers(true);};
    inline void SetUseDockerContainersOff() {SetUseDockerContainers(false);};
    inline void SetDockerImage(QString dockerimage) {_dockerimage = dockerimage;};
    inline QString GetDockerImage() {return _dockerimage;};
    inline void SetDockerImageOpenCarp(){_dockerimage = "docker.opencarp.org/opencarp/opencarp:latest";};
    inline void SetDockerImageUac(QString uac_tag="latest"){_dockerimage = "cemrg/uac:"+uac_tag;};// modify when docker image has been pushed to hub
    QStringList GetDockerArguments(QString volume, QString dockerexe = "");
    QStringList GetOpenCarpDockerDefaultArguments(QString volume);

    //Helper Functions
    bool CheckForStartedProcess();
    void ExecuteTouch(QString filepath);
    bool IsOutputSuccessful(QString outputFullPath);
    std::string PrintFullCommand(QString command, QStringList arguments);
    bool ExecuteCommand(QString executableName, QStringList arguments, QString outputPath, bool isOutputFile = true);

protected slots:

    void UpdateStdText();
    void UpdateErrText();
    void FinishedAlert();

private:

    //Dial and panel
    QDialog* dial;
    QTextEdit* panel;
    QVBoxLayout* layout;

    //QProcess
    bool completion;
    QString _dockerimage;
    bool _useDockerContainers, _debugvar;
    std::unique_ptr<QProcess> process;
};

#endif // CemrgCommandLine_h
