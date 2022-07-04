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
 * Eikonal Activation Simulation (EASI) Plugin for MITK
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/


#ifndef powertransViewPlot_h
#define powertransViewPlot_h

#include <berryISelectionListener.h>
#include <QmitkAbstractView.h>
#include <QmitkPlotWidget.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkSmartPointer.h>
#include <vtkColorTransferFunction.h>
#include <vtkRenderWindowInteractor.h>
#include "CemrgStrains.h"
#include "ui_powertransViewPlotControls.h"

#include "QmitkRenderWindow.h"
#include "mitkCommon.h"
#include "mitkDataStorage.h"
#include "mitkDataNode.h"
#include "mitkSurface.h"
#include "vtkRenderer.h"
#include "vtkTextActor.h"

/**
  \brief powertransViewPlot

  \warning  This class is not yet documented. Use "git blame" and ask the author to provide basic documentation.

  \sa QmitkAbstractView
  \ingroup ${plugin_target}_internal
*/
class powertransViewPlot: public QmitkAbstractView {
    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)
    Q_OBJECT

public:
    static const std::string VIEW_ID;
    static void SetDirectory(const QString directory);
    static void SetRibSpacing(int ribSpacing);

    powertransViewPlot();

protected:
    virtual void CreateQtPartControl(QWidget *parent) override;

    virtual void SetFocus() override;

    /// \brief called by QmitkFunctionality when DataManager's selection has changed
    virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer source, const QList<mitk::DataNode::Pointer>& nodes) override;

    Ui::powertransViewPlotControls m_Controls;

    /// \brief Called when the user clicks the GUI button
    void PlotData();
    void BullPlot();
    void FilePlot();
    void ColourAHASegments(int);

private:
    void HandleBullPlot(bool global);
    void HandleCurvPlot();
    void DrawAHALines();
    void DrawAHASegments(int frame, double* range);
    void DrawAHATextInfo();
    void WritePlotToCSV(QString dir);
    void WritePlotToVTK(QString dir);
    vtkSmartPointer<vtkColorTransferFunction> GetLookupTable(double *range);

    int cardiCycle;
    static int noFrames;
    static int smoothness;
    static QString directory;
    std::unique_ptr<QwtLegend> legend;
    std::unique_ptr<CemrgStrains> strain;
    vtkSmartPointer<vtkCamera> AHA_camera;
    vtkSmartPointer<vtkRenderer> AHA_renderer;
    vtkSmartPointer<vtkRenderWindowInteractor> AHA_interactor;
    std::vector<std::vector<double>> plotValueVectors;
    std::vector<vtkSmartPointer<vtkFloatArray>> flatPlotScalars;
    std::map<int, int> AHA;
};

#endif // powertransViewPlot_h
