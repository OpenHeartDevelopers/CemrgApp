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
 * Morphological Quantification
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>
#include <berryIWorkbenchPage.h>

// Qmitk
#include <mitkImage.h>
#include <mitkIOUtil.h>
#include <mitkProgressBar.h>
#include <mitkNodePredicateProperty.h>
#include <mitkManualSegmentationToSurfaceFilter.h>
#include "AtrialFibresVisualiseView.h"
#include "AtrialFibresView.h"

// VTK
#include <vtkGlyph3D.h>
#include <vtkSphere.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkProperty.h>
#include <vtkCellPicker.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCell.h>
#include <vtkMath.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkFloatArray.h>
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>
#include <vtkCellDataToPointData.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkClipPolyData.h>
#include <vtkHedgeHog.h>
#include <vtkNamedColors.h>

//ITK
#include <itkAddImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include "itkLabelObject.h"
#include "itkLabelMap.h"
#include "itkLabelImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelSelectionLabelMapFilter.h"

// Qt
#include <QMessageBox>
#include <QDesktopWidget>
#include <QFileDialog>

// CemrgAppModule
#include <CemrgAtriaClipper.h>
#include <CemrgCommandLine.h>
#include <CemrgCommonUtils.h>

QString AtrialFibresVisualiseView::fileName;
QString AtrialFibresVisualiseView::directory;

const std::string AtrialFibresVisualiseView::VIEW_ID = "org.mitk.views.AtrialFibresVisualiseView";

void AtrialFibresVisualiseView::CreateQtPartControl(QWidget *parent) {

    // create GUI widgets from the Qt Designer's .ui file
    m_Controls.setupUi(parent);
    connect(m_Controls.button_guide1, SIGNAL(clicked()), this, SLOT(Help()));
    connect(m_Controls.button_loadfibres, SIGNAL(clicked()), this, SLOT(LoadFibreFile()));

    //Setup renderer
    surfActor = vtkSmartPointer<vtkActor>::New();
    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(0.5,0.5,0.5);
    renderer->AutomaticLightCreationOn();
    renderer->LightFollowCameraOn();
    // renderer->TwoSidedLightingOn();
    // renderer->UpdateLightsGeometryToFollowCamera();
    vtkSmartPointer<vtkTextActor> txtActor = vtkSmartPointer<vtkTextActor>::New();
    std::string shortcuts = "CLICK ON Load Fibre File BUTTON";
    txtActor->SetInput(shortcuts.c_str());
    txtActor->GetTextProperty()->SetFontSize(14);
    txtActor->GetTextProperty()->SetColor(1.0, 1.0, 1.0);
    renderer->AddActor2D(txtActor);

    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow =
            vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    m_Controls.widget_1->SetRenderWindow(renderWindow);
    m_Controls.widget_1->GetRenderWindow()->AddRenderer(renderer);

    //Setup keyboard interactor
    // callBack = vtkSmartPointer<vtkCallbackCommand>::New();
    // callBack->SetCallback(KeyCallBackFunc);
    // callBack->SetClientData(this);
    interactor = m_Controls.widget_1->GetRenderWindow()->GetInteractor();
    interactor->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New());
    interactor->GetInteractorStyle()->KeyPressActivationOff();
    // interactor->GetInteractorStyle()->AddObserver(vtkCommand::KeyPressEvent, callBack);

    //Initialisation
    iniPreSurf();
    if (surface.IsNotNull()) {
        Visualiser();
    }
}

void AtrialFibresVisualiseView::SetFocus() {
    m_Controls.button_guide1->setFocus();
}

void AtrialFibresVisualiseView::OnSelectionChanged(
        berry::IWorkbenchPart::Pointer /*src*/, const QList<mitk::DataNode::Pointer>& /*nodes*/) {
}

AtrialFibresVisualiseView::~AtrialFibresVisualiseView() {
}

void AtrialFibresVisualiseView::LoadFibreFile(){
    QString fibresFilePath = QFileDialog::getOpenFileName(NULL, "Open Fibres file",
        directory.toStdString().c_str(), tr("Vector Space (*.vpts, *.vec)"));

    vtkSmartPointer<vtkPolyData> pd = surface->GetVtkPolyData();
    int numElems = pd->GetNumberOfCells(); // should match number of fibres in file
    MITK_INFO << ("Number of elements: " + QString::number(numElems)).toStdString();

    // read fibre file: std::vector<double> vectorfield=CemrgCommonUtils::ReadVectorField
    std::vector<double> vectorField=CemrgCommonUtils::ReadVectorField(fibresFilePath, true);
    MITK_INFO << ("Fibres read. Number of fibres: " + QString::number(vectorField.size())).toStdString();

    // rearrange vectorfield into vtkFlostArray fibres
    vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();
    vtkSmartPointer<vtkFloatArray> fibres = vtkSmartPointer<vtkFloatArray>::New();
    fibres->SetNumberOfComponents(3);
    fibres->SetNumberOfTuples(vectorField.size());
    double v[3];
    for (unsigned int ix = 0; ix < vectorField.size(); ix++) {
        v[0] = vectorField.at(3*ix + 0);
        v[1] = vectorField.at(3*ix + 1);
        v[2] = vectorField.at(3*ix + 2);

        fibres->InsertTuple(ix, v);
    }

    pd->GetCellData()->SetVectors(fibres);

    vtkSmartPointer<vtkHedgeHog> fibresHh = vtkSmartPointer<vtkHedgeHog>::New();
    fibresHh->SetInputData(pd);
    fibresHh->SetScaleFactor(0.1); // dunno

    vtkNew<vtkPolyDataMapper> fibresMapper;
    fibresMapper->SetInputConnection(fibresHh->GetOutputPort());

    vtkNew<vtkActor> fibresActor;
    fibresActor->SetMapper(fibresMapper);
    fibresActor->GetProperty()->SetColor(colors->GetColor3d("Gold").GetData());

    renderer->AddActor(fibresActor);
}

// slots
void AtrialFibresVisualiseView::Help(){
    std::string msg = "";
    msg = "HELP\n Visualise fibres on your mesh.\n ";
    msg += "Click Load Fibre File and select the corresponding vector file (vpts)";

    QMessageBox::information(NULL, "Help", msg.c_str());
}

void AtrialFibresVisualiseView::SetDirectoryFile(const QString directory, const QString fileName) {
    AtrialFibresVisualiseView::fileName = fileName;
    AtrialFibresVisualiseView::directory = directory;
}

void AtrialFibresVisualiseView::iniPreSurf() {
    //Find the selected node
    QString path = AtrialFibresVisualiseView::directory + "/" + AtrialFibresVisualiseView::fileName;
    // mitk::Surface::Pointer shell = mitk::IOUtil::Load<mitk::Surface>(path.toStdString());
    mitk::Surface::Pointer shell = CemrgCommonUtils::LoadVTKMesh(path.toStdString());
    CemrgCommonUtils::FlipXYPlane(shell, "", "");
    mitk::IOUtil::Save(shell, path.toStdString());

    surface = shell;
}

void AtrialFibresVisualiseView::Visualiser(double opacity){
    MITK_INFO << "[Visualiser]";
    double max_scalar=-1, min_scalar=1e9;
    // vtkFloatArray *scalars = vtkFloatArray::New();
    vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
    svtkFloatArray *scalars = vtkFloatArray::SafeDownCast(surface->GetVtkPolyData()->GetCellData()->GetScalars());
    for (vtkIdType i=0;i<surface->GetVtkPolyData()->GetNumberOfCells();i++) {
        double s = scalars->GetTuple1(i);
        if (s > max_scalar)
            max_scalar = s;
        if (s < min_scalar)
            min_scalar = s;
    }
    this->maxScalar = max_scalar;
    this->minScalar = min_scalar;

    //Create a mapper and actor for surface
    vtkSmartPointer<vtkPolyDataMapper> surfMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    surfMapper->SetInputData(surface->GetVtkPolyData());
    surfMapper->SetScalarRange(min_scalar, max_scalar);
    surfMapper->SetScalarModeToUseCellData();
    surfMapper->ScalarVisibilityOn();
    lut->SetTableRange(min_scalar, max_scalar);
    lut->SetHueRange(0.5, 0.0);  // this is the way_neighbourhood_size you tell which colors you want to be displayed.
    lut->Build();     // this is important

    vtkSmartPointer<vtkScalarBarActor> scalarBar =
            vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBar->SetLookupTable(surfMapper->GetLookupTable());
    scalarBar->SetWidth(0.075);
    scalarBar->SetHeight(0.3);
    scalarBar->SetTextPositionToPrecedeScalarBar();
    scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    scalarBar->GetPositionCoordinate()->SetValue( 0.9, 0.01 );
    scalarBar->SetNumberOfLabels(5);

    std::string scalarBarTitle = "Atrium Labels";
    scalarBar->SetTitle(scalarBarTitle.c_str());
    scalarBar->SetVerticalTitleSeparation(15);

    surfMapper->SetLookupTable(lut);
    scalarBar->SetLookupTable(lut);

    vtkSmartPointer<vtkActor> surfActor = vtkSmartPointer<vtkActor>::New();
    surfActor->SetMapper(surfMapper);
    surfActor->GetProperty()->SetOpacity(opacity);
    renderer->AddActor(surfActor);
    renderer->AddActor2D(scalarBar);
}

void AtrialFibresVisualiseView::SphereSourceVisualiser(vtkSmartPointer<vtkPolyData> pointSources, QString colour, double scaleFactor){
    // e.g colour = "0.4,0.1,0.0" - values for r,g, and b separated by commas.
    double r, g, b;
    bool okr, okg, okb;
    r = colour.section(',',0,0).toDouble(&okr);
    g = colour.section(',',1,1).toDouble(&okg);
    b = colour.section(',',2,2).toDouble(&okb);

    if(!okr) r=1.0;
    if(!okg) g=0.0;
    if(!okb) b=0.0;

    vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
    vtkSmartPointer<vtkSphereSource> glyphSource = vtkSmartPointer<vtkSphereSource>::New();
    glyph3D->SetInputData(pointSources);
    glyph3D->SetSourceConnection(glyphSource->GetOutputPort());
    glyph3D->SetScaleModeToDataScalingOff();
    glyph3D->SetScaleFactor(surface->GetVtkPolyData()->GetLength()*scaleFactor);
    glyph3D->Update();

    //Create a mapper and actor for glyph
    vtkSmartPointer<vtkPolyDataMapper> glyphMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    glyphMapper->SetInputConnection(glyph3D->GetOutputPort());
    vtkSmartPointer<vtkActor> glyphActor = vtkSmartPointer<vtkActor>::New();
    glyphActor->SetMapper(glyphMapper);
    glyphActor->GetProperty()->SetColor(r,g,b);
    glyphActor->PickableOff();
    renderer->AddActor(glyphActor);
}
