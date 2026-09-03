/*=========================================================================BSD 3-Clause License

Copyright (c) 2020, OpenHeartDevelopers
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=========================================================================*/

/*=========================================================================
*
* CemrgApp Plugin - Cemrg Atrial Strain Motion
*
* Cardiac Electromechanics Research Group
* http://www.cemrgapp.com
* info@cemrgapp.com 
*
* This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
* 
=========================================================================*/



#ifndef AtrialStrainMotionView_h
#define AtrialStrainMotionView_h

#include <berryISelectionListener.h>

#include <QmitkAbstractView.h>

#include "ui_AtrialStrainMotionViewControls.h"
#include "ui_AtrialFibresViewUIAnalysisSelector.h"
#include "ui_AtrialFibresViewUIMeshing.h"
#include "ui_AtrialFibresViewUIUacSelector.h"
#include "ui_AtrialFibresViewUIEditLabels.h"


#include "CemrgAtrialTools.h"

// RunAndCheckDockerAtrialStrainMotionSubcommand takes a CemrgCommandLine pointer.
// For a pointer, the compiler only needs the name of the class, not the full definition.
// An include of CemrgCommandLine.h also works, but is unnecessary in this header.
class CemrgCommandLine;

/**
  \brief AtrialStrainMotionView

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
class AtrialStrainMotionView : public QmitkAbstractView {
  // this is needed for all Qt objects that should have a Qt meta-object
  // (everything that derives from QObject and wants to have signal/slots)
  Q_OBJECT

public:
  static const std::string VIEW_ID;
  bool RequestProjectDirectoryFromUser();
  bool GetUserAnalysisSelectorInputs();
  void AutomaticAnalysis();
  bool GetUserMeshingInputs();

  inline QString Path(QString fnameExt=""){return (directory+"/"+fnameExt);};
  inline std::string StdStringPath(QString fnameExt=""){return (Path(fnameExt).toStdString());};
  inline void SetAutomaticPipeline(bool isAuto){automaticPipeline=isAuto;};

  int Ask(std::string title, std::string msg);
  bool LoadSurfaceChecks();
  void SetLgeAnalysis(bool b);
  bool UserLoadSurface();
  QString GetFilePath(QString type, QString extension);
  bool CheckLoadedMeshQuality();
  void SetTagNameFromPath(QString path);
  QString UserIncludeLgeAnalysis(QString segPath, ImageType::Pointer seg);
  bool GetUserUacOptionsInputs(bool enableFullUiOptions=true);
  bool GetUserEditLabelsInputs();
  bool IsOutputFileCorrect(QString dir, QStringList filenames);
  bool RunAndCheckDockerAtrialStrainMotionSubcommand(CemrgCommandLine* cmd, QString function, QString expectedOutput);

protected:
  virtual void CreateQtPartControl(QWidget *parent) override;

  virtual void SetFocus() override;

  /// \brief Overridden from QmitkAbstractView
  virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer source,
                                  const QList<mitk::DataNode::Pointer> &nodes) override;

  Ui::AtrialStrainMotionViewControls m_Controls;
  Ui::AtrialFibresViewUIAnalysisSelector m_UISelector;
  Ui::AtrialFibresViewUIMeshing m_UIMeshing;
  Ui::AtrialFibresViewUIUacSelector m_UIUac;
  Ui::AtrialFibresViewUIEditLabels m_UIEditLabels;

  /// \brief Called when the user clicks the GUI button


protected slots:
  // Model creation pipeline
  void SegmentExtract();
  void SegmentationPostprocessing();
  void IdentifyPV();
  void CreateLabelledMesh();
  void MeshPreprocessing();
  void ClipperPV();
  void ClipperMV();
  void UacCalculationVerifyLabels();
  void MeshImprovement();
  void AutoLandMark();
  void UAC_Stage1();
  void UAC_Stage2();
  void CreateModel();
  // motion tracking pipeline
  void CropImage();
  void Registration();
  void DeformMesh();
  void GenerateCellAreaStrains();
  void JacobianThreshold();
  void PlotAreaStrain();
  void CalcFiberStrains();
  void PlotFibersTrains();
  void Reset();

private:
   void SetDefaultPanelState();
   bool LoadSegmentationIntoStorage();

  double uiMesh_th, uiMesh_bl, uiMesh_smth, uiMesh_iter;
  QString directory, tagName, cnnPath;
  std::unique_ptr<CemrgAtrialTools> atrium;
  int uiSelector_pipeline;
  bool uiUac_meshtype_labelled;
  int uiUac_whichAtriumIndex, uiUac_fibreFileIndex, uiUac_surftypeIndex;
  bool automaticPipeline, analysisOnLge, resurfaceMesh, userHasSetMeshingParams;
  bool uiSelector_imgauto_skipCemrgNet, uiSelector_imgauto_skipLabel, uiSelector_img_scar, uiSelector_man_useCemrgNet;
  QStringList uiUac_fibreFile, uiUac_whichAtrium, uiUac_surftype;
  QStringList uiLabels;
  QString uac_whichAtrium;

};

#endif // AtrialStrainMotionView_h
