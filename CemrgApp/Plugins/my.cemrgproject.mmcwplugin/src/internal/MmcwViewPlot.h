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
 * Motion Measurement of Cardiac Wall (MMCW) Plugin for MITK
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef MmcwViewPlot_h
#define MmcwViewPlot_h

#include <berryISelectionListener.h>
#include <QmitkAbstractView.h>
#include <QmitkPlotWidget.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkSmartPointer.h>
#include <vtkColorTransferFunction.h>
#include "CemrgStrains.h"
#include "ui_MmcwViewControlsPlot.h"


/*!
  \brief MmcwViewPlot

  \warning  This application module is not yet documented. Use "svn blame/praise/annotate" and ask the author to provide basic documentation.

  \sa QmitkFunctionality
  \ingroup Functionalities
*/
class MmcwViewPlot : public QmitkAbstractView {

    // this is needed for all Qt objects that should have a Qt meta-object
    // (everything that derives from QObject and wants to have signal/slots)
    Q_OBJECT

public:

    static const std::string VIEW_ID;
    virtual void CreateQtPartControl(QWidget *parent);
    static void SetDirectory(const QString directory);
    static void SetNoFrames(int frames, int smoothness);

protected slots:

    /// \brief Called when the user clicks the GUI button
    void PlotData();
    void BullPlot();
    void FilePlot();
    void ColourAHASegments(int);

protected:

    virtual void SetFocus();
    /// \brief called by QmitkAbstractView when DataManager's selection has changed
    virtual void OnSelectionChanged(berry::IWorkbenchPart::Pointer source, const QList<mitk::DataNode::Pointer>& nodes);
    Ui::MmcwViewControlsPlot m_Controls;

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
    std::map<int,int> AHA;
};

#endif // MmcwViewPlot_h
