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
 * Atrial Clipper
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// Qmitk
#include <mitkProgressBar.h>
#include <mitkIOUtil.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkMorphologicalOperations.h>
#include "CemrgAtriaClipper.h"

// VTK
#include <vtkPointLocator.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkClipPolyData.h>
#include <vtkImplicitBoolean.h>
#include <vtkCenterOfMass.h>
#include <vtkPlane.h>
#include <vtkSphere.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkLinearExtrusionFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolygon.h>
#include <vtkCellArray.h>

// ITK
#include <itkSubtractImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkLabelShapeKeepNObjectsImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkImageDuplicator.h>
#include <itkResampleImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

// Qt
#include <QDebug>
#include <QString>
#include <QFileInfo>

// CemrgApp
#include "CemrgCommonUtils.h"
#include "CemrgMeasure.h"


CemrgAtriaClipper::CemrgAtriaClipper(QString directory, mitk::Surface::Pointer surface) {

    this->directory = directory;
    this->surface = surface;
    this->clippedSurface = surface;
    this->clippedSegImage = mitk::Image::New();
    this->ctrlnOrientation = false;
}

bool CemrgAtriaClipper::ComputeCtrLines(std::vector<int> pickedSeedLabels, vtkSmartPointer<vtkIdList> pickedSeedIds, bool autoLines) {

    try {

        MITK_INFO << "Producibility test. ";
        QString prodPath = directory + "/";
        mitk::IOUtil::Save(surface, (prodPath + "prodLineSurface.vtk").toStdString());
        ofstream prodFile1;
        prodFile1.open((prodPath + "prodSeedLabels.txt").toStdString());
        for (unsigned int i = 0; i < pickedSeedLabels.size(); i++)
            prodFile1 << pickedSeedLabels.at(i) << "\n";
        prodFile1.close();
        ofstream prodFile2;
        prodFile2.open((prodPath + "prodSeedIds.txt").toStdString());
        for (unsigned int i = 0; i < pickedSeedIds->GetNumberOfIds(); i++)
            prodFile2 << pickedSeedIds->GetId(i) << "\n";
        prodFile2.close();
        ofstream prodFile3;
        prodFile3.open((prodPath + "prodLineFlip.txt").toStdString());
        prodFile3 << autoLines << "\n";
        prodFile3.close();

        if (centreLines.size() == 0) {

            //Prepare source and target seeds
            vtkSmartPointer<vtkIdList> inletSeedIds = vtkSmartPointer<vtkIdList>::New();
            vtkSmartPointer<vtkIdList> outletSeedIds = vtkSmartPointer<vtkIdList>::New();

            MITK_INFO << "Determining centre lines' orientation.";
            vtkIdType centreOfMassId = CentreOfMass(surface);
            inletSeedIds->InsertNextId(centreOfMassId);
            MITK_INFO << "Number of pickedSeedLabels: ";
            MITK_INFO << pickedSeedLabels.size();
            MITK_INFO(!autoLines) << "Centre lines orientation set manually.";
            MITK_INFO(autoLines) << "Centre lines orientation set automatically.";

            for (unsigned int i = 0; i < pickedSeedLabels.size(); i++) {

                //Compute Centre Lines
                outletSeedIds->InsertNextId(pickedSeedIds->GetId(i));
                vtkSmartPointer<vtkvmtkPolyDataCenterlines> centreLineFilter = vtkSmartPointer<vtkvmtkPolyDataCenterlines>::New();
                centreLineFilter->SetInputData(surface->GetVtkPolyData());
                centreLineFilter->SetSourceSeedIds(inletSeedIds);
                centreLineFilter->SetTargetSeedIds(outletSeedIds);
                centreLineFilter->SetRadiusArrayName("MaximumInscribedSphereRadius");
                centreLineFilter->SetCostFunction("1/R");
                centreLineFilter->SetFlipNormals(autoLines ? ctrlnOrientation : !ctrlnOrientation);
                centreLineFilter->SetAppendEndPointsToCenterlines(0);
                centreLineFilter->SetSimplifyVoronoi(0);
                centreLineFilter->SetCenterlineResampling(1);
                centreLineFilter->SetResamplingStepLength(clSpacing);
                centreLineFilter->Update();
                outletSeedIds->DeleteId(pickedSeedIds->GetId(i));

                //Centrelines labels
                vtkSmartPointer<vtkIntArray> label = vtkSmartPointer<vtkIntArray>::New();
                label->SetNumberOfComponents(1);
                label->SetName("PickedSeedLabels");
                label->InsertNextValue(pickedSeedLabels.at(i));
                centreLineFilter->GetOutput()->GetFieldData()->AddArray(label);
                centreLines.push_back(centreLineFilter);

            }//_for
        }//_if

    } catch (...) {
        return false;
    }//_try
    return true;
}

bool CemrgAtriaClipper::ComputeCtrLinesClippers(std::vector<int> pickedSeedLabels) {

    //Compute centreline cut points
    manuals.clear();
    normalPlAngles.clear();
    centreLineVeinPlanes.clear();
    centreLinePolyPlanes.clear();
    centreLinePointPlanes.clear();
    vtkSmartPointer<vtkDoubleArray> areaMeter = vtkSmartPointer<vtkDoubleArray>::New();

    try {

        for (unsigned int i = 0; i < pickedSeedLabels.size(); i++) {

            int pointID;
            double slope;
            int clipPointID;
            int highCount = 0;
            vtkSmartPointer<vtkPolyData> line = centreLines.at(i)->GetOutput();
            int noBumpCriterion = round(criterion * line->GetNumberOfPoints());

            //Initialisation
            manuals.push_back(0);
            std::vector<double> tempVec = {0.0, 0.0};
            normalPlAngles.push_back(tempVec);

            //Diameters and area changes
            vtkSmartPointer<vtkvmtkPolyDataCenterlineSections> ctrLineSects = vtkSmartPointer<vtkvmtkPolyDataCenterlineSections>::New();
            ctrLineSects->SetInputData(surface->GetVtkPolyData());
            ctrLineSects->SetCenterlines(line);
            ctrLineSects->SetCenterlineSectionAreaArrayName("CentrelineSectionAreaArrayName");
            ctrLineSects->SetCenterlineSectionMinSizeArrayName("CenterlineSectionMinSizeArrayName");
            ctrLineSects->SetCenterlineSectionMaxSizeArrayName("CenterlineSectionMaxSizeArrayName");
            ctrLineSects->SetCenterlineSectionShapeArrayName("CenterlineSectionShapeArrayName");
            ctrLineSects->SetCenterlineSectionClosedArrayName("CenterlineSectionClosedArrayName");
            ctrLineSects->Update();
            centreLineVeinPlanes.push_back(ctrLineSects->GetOutput());

            //Slope calculations
            areaMeter = vtkDoubleArray::SafeDownCast(line->GetPointData()->GetArray("CentrelineSectionAreaArrayName"));
            for (pointID = 1; pointID < areaMeter->GetNumberOfTuples(); pointID++) {
                if (areaMeter->GetValue(pointID - 1) == 0) continue;
                slope = areaMeter->GetValue(pointID) - areaMeter->GetValue(pointID - 1);
                slope > highSlope ? highCount += 1 : highCount = 0;
                if (slope > maxiSlope) break;
                else if (slope > highSlope && highCount == noBumpCriterion) break;
            }//_for
            clipPointID = (highCount == 0) ? pointID - 1 : pointID - highCount;

            //Create a circle
            vtkSmartPointer<vtkRegularPolygonSource> polygonSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
            polygonSource->SetNumberOfSides(50);
            CalcParamsOfPlane(polygonSource, i, clipPointID);
            polygonSource->Update();
            centreLinePolyPlanes.push_back(polygonSource);

            //Fill in manual cutter
            centreLinePointPlanes.push_back(vtkSmartPointer<vtkPoints>::New());
        }//_for

    } catch (...) {
        return false;
    }//_try
    return true;
}

void CemrgAtriaClipper::ClipVeinsMesh(std::vector<int> pickedSeedLabels) {

    for (unsigned int i = 0; i < pickedSeedLabels.size(); i++) {

        //Label is not appendage-uncut
        if (pickedSeedLabels.at(i) != APPENDAGEUNCUT) {

            //Clipper
            vtkSmartPointer<vtkImplicitBoolean> implicitCircle = vtkSmartPointer<vtkImplicitBoolean>::New();
            vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
            plane->SetOrigin(centreLinePolyPlanes.at(i)->GetCenter());
            plane->SetNormal(centreLinePolyPlanes.at(i)->GetNormal());
            vtkSmartPointer<vtkSphere> sphere = vtkSmartPointer<vtkSphere>::New();
            sphere->SetCenter(centreLinePolyPlanes.at(i)->GetCenter());
            sphere->SetRadius(centreLinePolyPlanes.at(i)->GetRadius());
            implicitCircle->SetOperationTypeToIntersection();
            implicitCircle->AddFunction(plane);
            implicitCircle->AddFunction(sphere);
            vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
            clipper->SetClipFunction(implicitCircle);
            clipper->SetInputData(surface->GetVtkPolyData());
            clipper->InsideOutOff();
            clipper->Update();

            //Extract and clean surface mesh
            vtkSmartPointer<vtkDataSetSurfaceFilter> surfer = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
            surfer->SetInputData(clipper->GetOutput());
            surfer->Update();
            vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
            cleaner->SetInputConnection(surfer->GetOutputPort());
            cleaner->Update();
            vtkSmartPointer<vtkPolyDataConnectivityFilter> lrgRegion = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
            lrgRegion->SetInputConnection(cleaner->GetOutputPort());
            lrgRegion->SetExtractionModeToLargestRegion();
            lrgRegion->Update();
            cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
            cleaner->SetInputConnection(lrgRegion->GetOutputPort());
            cleaner->Update();
            clippedSurface->SetVtkPolyData(cleaner->GetOutput());

        }//_if
    }//_for

    //Save clipped mesh
    QString path = directory + "/segmentation.vtk";
    mitk::IOUtil::Save(clippedSurface, path.toStdString());
}

void CemrgAtriaClipper::ClipVeinsImage(std::vector<int> pickedSeedLabels, mitk::Image::Pointer segImage, bool morphAnalysis) {

    //Type definitions for new cut seg images
    typedef itk::Image<short, 3> ImageType;
    typedef itk::ImageRegionIteratorWithIndex<ImageType> ItType;
    typedef itk::SubtractImageFilter<ImageType, ImageType, ImageType> SubtractFilterType;
    typedef itk::BinaryBallStructuringElement<ImageType::PixelType, 3> BallType;
    typedef itk::GrayscaleDilateImageFilter<ImageType, ImageType, BallType> DilationFilterType;
    typedef itk::ImageDuplicator<ImageType> DuplicatorType;

    //Cast Seg to ITK formats
    ImageType::Pointer segItkImage = ImageType::New();
    CastToItkImage(segImage, segItkImage);
    ImageType::Pointer orgSegItkImage = ImageType::New();
    CastToItkImage(segImage, orgSegItkImage);
    ImageType::Pointer pvLblsItkImage = ImageType::New();
    CastToItkImage(segImage, pvLblsItkImage);
    std::vector<ImageType::Pointer> cutRegions;

    for (unsigned int i = 0; i < pickedSeedLabels.size(); i++) {

        //Find the right vein section by removing unwanted ones
        vtkSmartPointer<vtkPolyData> line = centreLines.at(i)->GetOutput();
        int position = line->FindPoint(centreLinePolyPlanes.at(i)->GetCenter());
        vtkSmartPointer<vtkPolyData> centreVeinPlane = centreLineVeinPlanes.at(i);
        centreVeinPlane->BuildLinks();
        for (vtkIdType cellID = 0; cellID < centreVeinPlane->GetNumberOfCells(); cellID++)
            if (cellID != line->GetNumberOfPoints() - 1 - position)
                centreVeinPlane->DeleteCell(cellID);
        centreVeinPlane->RemoveDeletedCells();
        vtkSmartPointer<vtkCleanPolyData> cleanPolyData = vtkSmartPointer<vtkCleanPolyData>::New();
        cleanPolyData->SetInputData(centreVeinPlane);
        cleanPolyData->Update();
        centreVeinPlane = cleanPolyData->GetOutput();

        //Morphological Analysis
        if (morphAnalysis) {

            int noPV = pickedSeedLabels.size();
            vtkSmartPointer<vtkDoubleArray> areas = vtkDoubleArray::SafeDownCast(line->GetPointData()->GetArray("CentrelineSectionAreaArrayName"));
            double area = areas->GetValue(position);
            CemrgMeasure::Points points;
            for (vtkIdType j = 0; j < centreVeinPlane->GetNumberOfPoints(); j++) {
                double p[3];
                centreVeinPlane->GetPoints()->GetPoint(j, p);
                points.push_back(CemrgMeasure::Point(p[0], p[1], p[2]));
            }//_for
            std::unique_ptr<CemrgMeasure> morphAnal = std::unique_ptr<CemrgMeasure>(new CemrgMeasure());
            double circ = morphAnal->CalcPerimeter(points);
            vtkSmartPointer<vtkDoubleArray> diams = vtkDoubleArray::SafeDownCast(line->GetPointData()->GetArray("CenterlineSectionMinSizeArrayName"));
            double diam = diams->GetValue(position);

            ofstream morphResult;
            QString morphPath = directory + "/morphResults.txt";
            morphResult.open(morphPath.toStdString(), std::ios_base::app);
            if (i == 0)
                morphResult << "NO " << --noPV << "\n";
            morphResult << pickedSeedLabels.at(i) << " " << area << "\n";
            morphResult << pickedSeedLabels.at(i) << " " << circ << "\n";
            morphResult << pickedSeedLabels.at(i) << " " << diam << "\n";
            morphResult.close();
        }//_if

        //Flip the cutter plane
        vtkSmartPointer<vtkPolyData> circle;
        if (manuals[i] == 2) {
            vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
            polygon->GetPointIds()->SetNumberOfIds(centreLinePointPlanes.at(i)->GetNumberOfPoints());
            for (int j = 0; j < centreLinePointPlanes.at(i)->GetNumberOfPoints(); j++)
                polygon->GetPointIds()->SetId(j, j);
            vtkSmartPointer<vtkCellArray> polygons = vtkSmartPointer<vtkCellArray>::New();
            polygons->InsertNextCell(polygon);
            vtkSmartPointer<vtkPolyData> polygonPolyData = vtkSmartPointer<vtkPolyData>::New();
            polygonPolyData->SetPoints(centreLinePointPlanes.at(i));
            polygonPolyData->SetPolys(polygons);
            circle = polygonPolyData;
            QString path = directory + "/manualType2Clipper.vtk";
            vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
            writer->SetInputData(circle);
            writer->SetFileName(path.toStdString().c_str());
            writer->Write();
        } else if (manuals[i] == 1)
            circle = centreLinePolyPlanes.at(i)->GetOutput();
        else
            circle = centreLineVeinPlanes.at(i);
        for (int j = 0; j < circle->GetNumberOfPoints(); j++) {
            double* point = circle->GetPoint(j);
            point[0] = -point[0];
            point[1] = -point[1];
            circle->GetPoints()->SetPoint(j, point);
        }//_for

        /*
         * Producibility Test
         **/
        try {
            QString prodPath = directory + "/";
            mitk::Surface::Pointer prodSurf = mitk::Surface::New();
            prodSurf->SetVtkPolyData(circle);
            mitk::IOUtil::Save(prodSurf, (prodPath + "prodCutter" + QString::number(i) + ".vtk").toStdString());
            ofstream prodFile1;
            prodFile1.open((prodPath + "prodCutter" + QString::number(i) + "TNormals.txt").toStdString());
            prodFile1 << centreLinePolyPlanes.at(i)->GetNormal()[0] << "\n";
            prodFile1 << centreLinePolyPlanes.at(i)->GetNormal()[1] << "\n";
            prodFile1 << centreLinePolyPlanes.at(i)->GetNormal()[2] << "\n";
            prodFile1 << manuals[i] << "\n";
            prodFile1.close();
        } catch (...) {};
        /*
         * End Test
         **/

         //Prepare empty image
        vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();
        double spacing[3];
        segImage->GetVtkImageData()->GetSpacing(spacing);
        whiteImage->SetSpacing(spacing);
        int dimensions[3];
        segImage->GetVtkImageData()->GetDimensions(dimensions);
        whiteImage->SetDimensions(dimensions);
        whiteImage->SetExtent(0, dimensions[0] - 1, 0, dimensions[1] - 1, 0, dimensions[2] - 1);
        double origin[3];
        segImage->GetGeometry()->GetOrigin().ToArray(origin);
        whiteImage->SetOrigin(origin);
        whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
        unsigned char otval = 0;
        unsigned char inval = 255;
        vtkIdType count = whiteImage->GetNumberOfPoints();
        for (vtkIdType j = 0; j < count; j++)
            whiteImage->GetPointData()->GetScalars()->SetTuple1(j, inval);

        //Sweep polygonal data to create an image
        vtkSmartPointer<vtkLinearExtrusionFilter> extruder = vtkSmartPointer<vtkLinearExtrusionFilter>::New();
        extruder->SetInputData(circle);
        extruder->SetScaleFactor(1.0);
        extruder->SetExtrusionTypeToNormalExtrusion();
        extruder->SetVector(centreLinePolyPlanes.at(i)->GetNormal());
        extruder->Update();
        vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
        pol2stenc->SetTolerance(0.5);
        pol2stenc->SetInputConnection(extruder->GetOutputPort());
        pol2stenc->SetOutputOrigin(origin);
        pol2stenc->SetOutputSpacing(spacing);
        pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
        pol2stenc->Update();
        vtkSmartPointer<vtkImageStencil> imgstenc = vtkSmartPointer<vtkImageStencil>::New();
        imgstenc->SetInputData(whiteImage);
        imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
        imgstenc->ReverseStencilOff();
        imgstenc->SetBackgroundValue(otval);
        imgstenc->Update();

        // VTK to ITK conversion
        mitk::Image::Pointer cutImg = mitk::Image::New();
        cutImg->Initialize(imgstenc->GetOutput());
        cutImg->SetVolume(imgstenc->GetOutput()->GetScalarPointer());
        ImageType::Pointer cutItkImage = ImageType::New();
        CastToItkImage(cutImg, cutItkImage);
        itk::ResampleImageFilter<ImageType, ImageType>::Pointer resampleFilter;
        resampleFilter = itk::ResampleImageFilter<ImageType, ImageType >::New();
        resampleFilter->SetInput(cutItkImage);
        resampleFilter->SetReferenceImage(segItkImage);
        resampleFilter->SetUseReferenceImage(true);
        resampleFilter->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<ImageType>::New());
        resampleFilter->SetDefaultPixelValue(0);
        resampleFilter->UpdateLargestPossibleRegion();
        cutItkImage = resampleFilter->GetOutput();

        //Image Dilation
        BallType binaryBall;
        binaryBall.SetRadius(1); // before, radius=(manuals[i] == 1 ? 1 : 1.5) change to a larger value if problems arise
        binaryBall.CreateStructuringElement();
        DilationFilterType::Pointer dilationFilter = DilationFilterType::New();
        dilationFilter->SetInput(cutItkImage);
        dilationFilter->SetKernel(binaryBall);
        dilationFilter->UpdateLargestPossibleRegion();
        cutImg = mitk::ImportItkImage(dilationFilter->GetOutput())->Clone();
        CastToItkImage(cutImg, cutItkImage);

        //Subtract images
        SubtractFilterType::Pointer subFilter = SubtractFilterType::New();
        subFilter->SetInput1(segItkImage);
        subFilter->SetInput2(cutItkImage);
        subFilter->UpdateLargestPossibleRegion();

        //Duplicate subtract images
        DuplicatorType::Pointer duplicator = DuplicatorType::New();
        duplicator->SetInputImage(subFilter->GetOutput());
        duplicator->Update();
        cutRegions.push_back(duplicator->GetOutput());

        //Record voxel locations before cut
        ItType itSub(subFilter->GetOutput(), subFilter->GetOutput()->GetRequestedRegion());
        ItType itLbl(pvLblsItkImage, pvLblsItkImage->GetRequestedRegion());
        itLbl.GoToBegin();
        for (itSub.GoToBegin(); !itSub.IsAtEnd(); ++itSub) {
            if ((int)itSub.Get() == -254 || (int)itSub.Get() == -255) {
                if (pickedSeedLabels.at(i) == APPENDAGEUNCUT && (int)itSub.Get() == -255)
                    itSub.Set(0);
                if (pickedSeedLabels.at(i) != APPENDAGEUNCUT)
                    itSub.Set(0);
                itLbl.Set(0);
            }//_if
            ++itLbl;
        }//_for

        //Relabel the components
        typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> ConnectedComponentImageFilterType;
        ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
        connected->SetInput(subFilter->GetOutput());
        connected->Update();
        typedef itk::RelabelComponentImageFilter<ImageType, ImageType> RelabelFilterType;
        RelabelFilterType::Pointer relabeler = RelabelFilterType::New();
        relabeler->SetInput(connected->GetOutput());
        relabeler->Update();

        //Keep the single largest component
        typedef itk::LabelShapeKeepNObjectsImageFilter<ImageType> LabelShapeKeepNObjImgFilterType;
        LabelShapeKeepNObjImgFilterType::Pointer lblShpKpNObjImgFltr = LabelShapeKeepNObjImgFilterType::New();
        lblShpKpNObjImgFltr->SetInput(relabeler->GetOutput());
        lblShpKpNObjImgFltr->SetBackgroundValue(0);
        lblShpKpNObjImgFltr->SetNumberOfObjects(1);
        lblShpKpNObjImgFltr->SetAttribute(LabelShapeKeepNObjImgFilterType::LabelObjectType::NUMBER_OF_PIXELS);
        lblShpKpNObjImgFltr->Update();
        segItkImage = lblShpKpNObjImgFltr->GetOutput();

    }//_for

    //Label individual veins
    typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> ConnectedComponentImageFilterType;
    ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
    connected->SetInput(pvLblsItkImage);
    connected->Update();
    typedef itk::RelabelComponentImageFilter<ImageType, ImageType> RelabelFilterType;
    RelabelFilterType::Pointer relabeler = RelabelFilterType::New();
    relabeler->SetInput(connected->GetOutput());
    relabeler->Update();
    pvLblsItkImage = relabeler->GetOutput();

    //Adjust voxel labels after cut MV
    ItType itCut(segItkImage, segItkImage->GetRequestedRegion());
    ItType itLbl(pvLblsItkImage, pvLblsItkImage->GetRequestedRegion());
    ItType itOrg(orgSegItkImage, orgSegItkImage->GetRequestedRegion());
    itCut.GoToBegin();
    itLbl.GoToBegin();
    for (itOrg.GoToBegin(); !itOrg.IsAtEnd(); ++itOrg) {
        //MV labels
        if ((int)itOrg.Get() == 2) {
            itCut.Set(2);
            itLbl.Set(10);
        }//_if
        ++itCut;
        ++itLbl;
    }//_for

    //Adjust voxel labels after cut PV
    for (unsigned int i = 0; i < pickedSeedLabels.size(); i++) {
        ItType itCutSub(segItkImage, segItkImage->GetRequestedRegion());
        ItType itLblSub(pvLblsItkImage, pvLblsItkImage->GetRequestedRegion());
        ItType itOrgSub(cutRegions.at(i), cutRegions.at(i)->GetRequestedRegion());
        itCutSub.GoToBegin();
        itLblSub.GoToBegin();
        for (itOrgSub.GoToBegin(); !itOrgSub.IsAtEnd(); ++itOrgSub) {
            //PV labels
            if ((int)itOrgSub.Get() == -254) {
                if (pickedSeedLabels.at(i) != APPENDAGEUNCUT)
                    itCutSub.Set(3);
                itLblSub.Set(pickedSeedLabels.at(i));
            }//_if
            ++itCutSub;
            ++itLblSub;
        }//_for
    }//_for

    //Save image with individual veins labelled
    QString path = directory + "/PVeinsLabelled.nii";
    mitk::Image::Pointer pvLabelled = mitk::ImportItkImage(pvLblsItkImage)->Clone();
    mitk::IOUtil::Save(pvLabelled, path.toStdString());
    if (morphAnalysis) {
        path = directory + "/AnalyticBloodpool.nii";
        mitk::IOUtil::Save(pvLabelled, path.toStdString());
    }//_if

    //Save clipped image
    path = directory + "/PVeinsCroppedImage.nii";
    mitk::Image::Pointer pvCropped = mitk::ImportItkImage(segItkImage)->Clone();
    mitk::IOUtil::Save(pvCropped, path.toStdString());
    clippedSegImage = pvCropped;
}

mitk::Image::Pointer CemrgAtriaClipper::CreateImageCutter(vtkSmartPointer<vtkPolyData> circle, double cline_normal[3], mitk::Image::Pointer segImage){
    // Type definitions for new cut seg images
    typedef itk::Image<short, 3> ImageType;
    typedef itk::BinaryBallStructuringElement<ImageType::PixelType, 3> BallType;
    typedef itk::GrayscaleDilateImageFilter<ImageType, ImageType, BallType> DilationFilterType;

    //Cast Seg to ITK formats
    ImageType::Pointer segItkImage = ImageType::New();
    ImageType::Pointer orgSegItkImage = ImageType::New();
    ImageType::Pointer pvLblsItkImage = ImageType::New();
    CastToItkImage(segImage, segItkImage);
    CastToItkImage(segImage, orgSegItkImage);
    CastToItkImage(segImage, pvLblsItkImage);
    std::vector<ImageType::Pointer> cutRegions;
    //Prepare empty image
    vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();
    double spacing[3];
    segImage->GetVtkImageData()->GetSpacing(spacing);
    whiteImage->SetSpacing(spacing);
    int dimensions[3];
    segImage->GetVtkImageData()->GetDimensions(dimensions);
    whiteImage->SetDimensions(dimensions);
    whiteImage->SetExtent(0, dimensions[0] - 1, 0, dimensions[1] - 1, 0, dimensions[2] - 1);
    double origin[3];
    segImage->GetGeometry()->GetOrigin().ToArray(origin);
    whiteImage->SetOrigin(origin);
    whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);
    unsigned char otval = 0;
    unsigned char inval = 255;
    vtkIdType count = whiteImage->GetNumberOfPoints();
    for (vtkIdType i = 0; i < count; ++i)
        whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);

    //Sweep polygonal data to create an image
    vtkSmartPointer<vtkLinearExtrusionFilter> extruder = vtkSmartPointer<vtkLinearExtrusionFilter>::New();
    extruder->SetInputData(circle);
    extruder->SetScaleFactor(1.0);
    extruder->SetExtrusionTypeToNormalExtrusion();
    extruder->SetVector(cline_normal);
    extruder->Update();
    vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
    pol2stenc->SetTolerance(0.5);
    pol2stenc->SetInputConnection(extruder->GetOutputPort());
    pol2stenc->SetOutputOrigin(origin);
    pol2stenc->SetOutputSpacing(spacing);
    pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
    pol2stenc->Update();
    vtkSmartPointer<vtkImageStencil> imgstenc = vtkSmartPointer<vtkImageStencil>::New();
    imgstenc->SetInputData(whiteImage);
    imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
    imgstenc->ReverseStencilOff();
    imgstenc->SetBackgroundValue(otval);
    imgstenc->Update();

    //VTK to ITK conversion
    mitk::Image::Pointer cutImg = mitk::Image::New();
    cutImg->Initialize(imgstenc->GetOutput());
    cutImg->SetVolume(imgstenc->GetOutput()->GetScalarPointer());
    ImageType::Pointer cutItkImage = ImageType::New();
    CastToItkImage(cutImg, cutItkImage);
    itk::ResampleImageFilter<ImageType, ImageType>::Pointer resampleFilter;
    resampleFilter = itk::ResampleImageFilter<ImageType, ImageType >::New();
    resampleFilter->SetInput(cutItkImage);
    resampleFilter->SetReferenceImage(segItkImage);
    resampleFilter->SetUseReferenceImage(true);
    resampleFilter->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<ImageType>::New());
    resampleFilter->SetDefaultPixelValue(0);
    resampleFilter->UpdateLargestPossibleRegion();
    cutItkImage = resampleFilter->GetOutput();

    //Image Dilation
    BallType binaryBall;
    binaryBall.SetRadius(static_cast<unsigned long>(2.0));
    binaryBall.CreateStructuringElement();
    DilationFilterType::Pointer dilationFilter = DilationFilterType::New();
    dilationFilter->SetInput(cutItkImage);
    dilationFilter->SetKernel(binaryBall);
    dilationFilter->UpdateLargestPossibleRegion();
    cutImg = mitk::ImportItkImage(dilationFilter->GetOutput())->Clone();

    return cutImg;
}

mitk::Image::Pointer CemrgAtriaClipper::SubtractAndLabel(mitk::Image::Pointer im, mitk::Image::Pointer cutIm){
    //Type definitions for new cut seg images
    typedef itk::Image<short, 3> ImageType;
    typedef itk::ImageRegionIteratorWithIndex<ImageType> ItType;
    typedef itk::SubtractImageFilter<ImageType, ImageType, ImageType> SubtractFilterType;

    //Cast Seg to ITK formats
    ImageType::Pointer segItkImage = ImageType::New();
    ImageType::Pointer cutItkImage = ImageType::New();
    ImageType::Pointer pvLblsItkImage = ImageType::New();
    CastToItkImage(im, segItkImage);
    CastToItkImage(im, pvLblsItkImage);
    CastToItkImage(cutIm, cutItkImage);

    //Subtract images
    SubtractFilterType::Pointer subFilter = SubtractFilterType::New();
    subFilter->SetInput1(segItkImage);
    subFilter->SetInput2(cutItkImage);
    subFilter->UpdateLargestPossibleRegion();

    //Record voxel locations before cut
    ItType itSub(subFilter->GetOutput(), subFilter->GetOutput()->GetRequestedRegion());
    ItType itLbl(pvLblsItkImage, pvLblsItkImage->GetRequestedRegion());
    itLbl.GoToBegin();
    for (itSub.GoToBegin(); !itSub.IsAtEnd(); ++itSub) {
        if ((int)itSub.Get() == -254 || (int)itSub.Get() == -255) {
                itSub.Set(0);
            itLbl.Set(0);
        }//_if
        ++itLbl;
    }//_for

    //Relabel the components
    typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> ConnectedComponentImageFilterType;
    ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
    connected->SetInput(subFilter->GetOutput());
    connected->Update();
    typedef itk::RelabelComponentImageFilter<ImageType, ImageType> RelabelFilterType;
    RelabelFilterType::Pointer relabeler = RelabelFilterType::New();
    relabeler->SetInput(connected->GetOutput());
    relabeler->Update();

    mitk::Image::Pointer returnIm = mitk::ImportItkImage(relabeler->GetOutput())->Clone();

    return returnIm;
}

void CemrgAtriaClipper::CreateCroppedAndLabelled(QString dir, mitk::Image::Pointer im, mitk::Image::Pointer pieces, mitk::Image::Pointer labcutters, QString croppedName, QString labelledName){
    //Type definitions for new cut seg images
    typedef itk::Image<short, 3> ImageType;
    typedef itk::ImageRegionIteratorWithIndex<ImageType> ItType;

    //Cast Seg to ITK formats
    ImageType::Pointer segItk = ImageType::New();
    ImageType::Pointer pveinsLabelledItk = ImageType::New();
    ImageType::Pointer pveinsCroppedItk = ImageType::New();
    ImageType::Pointer labcuttersItk = ImageType::New();
    CastToItkImage(im->Clone(), segItk);
    CastToItkImage(pieces->Clone(), pveinsLabelledItk);
    CastToItkImage(pieces->Clone(), pveinsCroppedItk);
    CastToItkImage(labcutters, labcuttersItk);

    ItType iterSeg(segItk, segItk->GetRequestedRegion());
    ItType iterCutters(labcuttersItk, labcuttersItk->GetRequestedRegion());
    ItType iterCropped(pveinsCroppedItk, pveinsCroppedItk->GetRequestedRegion());
    ItType iterLabelled(pveinsLabelledItk, pveinsLabelledItk->GetRequestedRegion());

    iterSeg.GoToBegin();
    iterCutters.GoToBegin();
    iterCropped.GoToBegin();
    iterLabelled.GoToBegin();

    while(!iterCutters.IsAtEnd()){
        int cuttersValue = (int) iterCutters.Get();
        int croppedValue = (int) iterCropped.Get();
        int segValue = (int) iterSeg.Get();

        if(cuttersValue > 0 && segValue > 0){
            iterCropped.Set(3);
            iterLabelled.Set(cuttersValue);
        } else if(croppedValue > 1){
            iterCropped.Set(0);
        }

        ++iterSeg;
        ++iterCutters;
        ++iterCropped;
        ++iterLabelled;
    }

    mitk::Image::Pointer croppedIm = mitk::ImportItkImage(pveinsCroppedItk);
    mitk::Image::Pointer labelledIm = mitk::ImportItkImage(pveinsLabelledItk);

    croppedName += (!croppedName.contains(".nii")) ? ".nii" : "";
    labelledName += (!labelledName.contains(".nii")) ? ".nii" : "";

    mitk::IOUtil::Save(croppedIm, (dir+"/"+croppedName).toStdString());
    mitk::IOUtil::Save(labelledIm, (dir+"/"+labelledName).toStdString());

}

void CemrgAtriaClipper::ClipVeinsImgFromFileList(QStringList cutters, QStringList normals, mitk::Image::Pointer im, QString outname){
    CemrgCommonUtils::Binarise(im);
    QFileInfo qfi(cutters.at(0));
    QString path2files = qfi.absolutePath();

    mitk::Image::Pointer cuttersAdd = CemrgCommonUtils::Zeros(im);
    mitk::Image::Pointer cuttersWithLabels = CemrgCommonUtils::Zeros(im);
    mitk::Image::Pointer originalIm = im->Clone();

    // Search for prodSeedLabels.txt - store labels, depends on the cutters.size()
    int numFiles = cutters.size();
    std::vector<int> labels(numFiles, 21); // assume all default
    QString labelsFile = qfi.absolutePath() + "/prodSeedLabels.txt";
    if(QFile::exists(labelsFile)){
        MITK_INFO << ("File found: " + labelsFile).toStdString();
        std::ifstream fl(labelsFile.toStdString());
        for (int jx = 0; jx < numFiles; jx++) {
            if(fl.eof()){
                MITK_INFO << "File finished prematurely.";
                break;
            }
            fl >> labels[jx];
        }
        fl.close();
    }

    for (int ix = 0; ix < numFiles; ix++) {
        std::ifstream fi(normals.at(ix).toStdString());
        double cline_normal[3];
        fi >> cline_normal[0];
        fi >> cline_normal[1];
        fi >> cline_normal[2];

        mitk::Surface::Pointer cutsurf = mitk::IOUtil::Load<mitk::Surface>(cutters.at(ix).toStdString());
        vtkSmartPointer<vtkPolyData> circle = cutsurf->GetVtkPolyData();

        mitk::Image::Pointer oneCutter = CreateImageCutter(circle, cline_normal, im);
        cuttersAdd = CemrgCommonUtils::AddImage(cuttersAdd, oneCutter);

        mitk::Image::Pointer oneCutterLabel = CemrgCommonUtils::ImageThreshold(oneCutter, 0, labels.at(ix), 0);
        cuttersWithLabels = CemrgCommonUtils::AddImage(cuttersWithLabels, oneCutterLabel);
    }

    mitk::Image::Pointer subIm = SubtractAndLabel(im, cuttersAdd);
    CreateCroppedAndLabelled(qfi.absolutePath(), originalIm, subIm, cuttersWithLabels, outname+"CroppedImage", outname+"Labelled");
}

void CemrgAtriaClipper::CalcParamsOfPlane(vtkSmartPointer<vtkRegularPolygonSource> plane, int ctrLineNo, int position) {

    vtkSmartPointer<vtkPolyData> line = centreLines.at(ctrLineNo)->GetOutput();
    vtkSmartPointer<vtkDoubleArray> radii = vtkDoubleArray::SafeDownCast(line->GetPointData()->GetArray("MaximumInscribedSphereRadius"));

    //Centre point
    double* clipPoint = line->GetPoint(position);

    //Radius
    double clipRadius = radii->GetValue(position) * radiusAdj;

    //Normal calculations
    double x_s = 0.0, y_s = 0.0, z_s = 0.0;
    double pointAtPosition[3], pointAtNextPosition[3];
    line->GetPoint(position + 1, pointAtNextPosition);
    line->GetPoint(position, pointAtPosition);

    x_s = pointAtNextPosition[0] - pointAtPosition[0];
    y_s = pointAtNextPosition[1] - pointAtPosition[1];
    z_s = pointAtNextPosition[2] - pointAtPosition[2];
    if (manuals[ctrLineNo] == 1) {
        double x, y, z;
        double lng = sqrt(pow(x_s, 2) + pow(y_s, 2) + pow(z_s, 2));
        x = lng * sin(normalPlAngles[ctrLineNo][0]) * cos(normalPlAngles[ctrLineNo][1]) + line->GetPoint(position)[0];
        y = lng * sin(normalPlAngles[ctrLineNo][0]) * sin(normalPlAngles[ctrLineNo][1]) + line->GetPoint(position)[1];
        z = lng * cos(normalPlAngles[ctrLineNo][0]) + line->GetPoint(position)[2];
        x_s = x - line->GetPoint(position)[0];
        y_s = y - line->GetPoint(position)[1];
        z_s = z - line->GetPoint(position)[2];
    } else {
        double y, z;
        double lng = sqrt(pow(x_s, 2) + pow(y_s, 2) + pow(z_s, 2));
        y = y_s + line->GetPoint(position)[1];
        z = z_s + line->GetPoint(position)[2];
        normalPlAngles[ctrLineNo][0] = acos((z - line->GetPoint(position)[2]) / lng);
        normalPlAngles[ctrLineNo][1] = asin((y - line->GetPoint(position)[1]) / lng * sin(normalPlAngles[ctrLineNo][0]));
    }//_if
    double clipNormal[3] = {x_s, y_s, z_s};

    //Adjust clipper
    plane->SetCenter(clipPoint);
    plane->SetRadius(clipRadius);
    plane->SetNormal(clipNormal);
}

void CemrgAtriaClipper::ResetCtrLinesClippingPlanes() {

    centreLineVeinPlanes.clear();
    centreLinePolyPlanes.clear();
    centreLinePointPlanes.clear();
}

vtkIdType CemrgAtriaClipper::CentreOfMass(mitk::Surface::Pointer surface) {

    //Polydata of surface
    vtkSmartPointer<vtkPolyData> pd = surface->GetVtkPolyData();
    CemrgCommonUtils::CalculatePolyDataNormals(pd, false);
    vtkSmartPointer<vtkFloatArray> pdNormals = vtkFloatArray::SafeDownCast(pd->GetPointData()->GetNormals());
    vtkSmartPointer<vtkCenterOfMass> centre = vtkSmartPointer<vtkCenterOfMass>::New();

    double centreInSurface[3];
    centre->SetInputData(pd);
    centre->SetUseScalarsAsWeights(false);
    centre->Update();
    centre->GetCenter(centreInSurface);

    vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
    pointLocator->SetDataSet(pd);
    pointLocator->BuildLocator();
    vtkIdType id = pointLocator->FindClosestPoint(centreInSurface); // output of method

    double * pointNormal = pdNormals->GetTuple(id);
    double * pointInSurface = pd->GetPoints()->GetPoint(id);
    double xProduct = 0;
    double mitralValvePlane[3] = {0}; //change name
    mitralValvePlane[0] = pointInSurface[0] - centreInSurface[0];
    mitralValvePlane[1] = pointInSurface[1] - centreInSurface[1];
    mitralValvePlane[2] = pointInSurface[2] - centreInSurface[2];
    for (int i = 0; i < 3; i++) {
        xProduct += (mitralValvePlane[i]) * (pointNormal[i]);
    }//_for

    //Orientation of centrelines
    if (xProduct < 0) {
        ctrlnOrientation = true;
    }//_if

    return id;
}

void CemrgAtriaClipper::VTKWriter(vtkSmartPointer<vtkPolyData> PD, QString path) {

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInputData(PD);
    writer->SetFileName(path.toStdString().c_str());
    writer->Write();
}

void CemrgAtriaClipper::SetMClipperAngles(double* value, int clippersIndex) {

    manuals[clippersIndex] = 1;
    normalPlAngles[clippersIndex][0] = value[0];
    normalPlAngles[clippersIndex][1] = value[1];
}

void CemrgAtriaClipper::SetMClipperSeeds(vtkSmartPointer<vtkPolyData> pickedCutterSeeds, int clippersIndex) {

    manuals[clippersIndex] = 2;
    centreLinePointPlanes.at(clippersIndex) = pickedCutterSeeds->GetPoints();
}

mitk::Image::Pointer CemrgAtriaClipper::FixClippingErrors(mitk::Image::Pointer image, int strel_radius){
    using ImageType = itk::Image<short, 3>;
    using IteratorType = itk::ImageRegionIterator<ImageType>;
    using ConnectedComponentImageFilterType = itk::ConnectedComponentImageFilter<ImageType, ImageType>;
    using LabelShapeKeepNObjImgFilterType = itk::LabelShapeKeepNObjectsImageFilter<ImageType>;

    // start routine
    mitk::Image::Pointer vimage = image->Clone();
    ImageType::Pointer veins = ImageType::New();

    ImageType::Pointer im = ImageType::New();
    mitk::CastToItkImage(image, im);

    // V2=imopen(V==1, strel('sphere', 3))
    IteratorType imIter(im, im->GetLargestPossibleRegion());
    imIter.GoToBegin();
    while (!imIter.IsAtEnd()){
        short value = (imIter.Get() == 1) ? 1 : 0;
        imIter.Set(value);

        ++imIter;
    }

    mitk::Image::Pointer remove_veins = mitk::Image::New();
    mitk::CastToMitkImage(im, remove_veins);

    try {
        mitk::MorphologicalOperations::StructuralElementType structuralElement = mitk::MorphologicalOperations::Ball;
        mitk::MorphologicalOperations::Opening(remove_veins, strel_radius, structuralElement);
    }
    catch (const itk::ExceptionObject &exception) {
        MITK_WARN << "Exception caught: " << exception.GetDescription();
        return NULL;
    }

    ImageType::Pointer imop = ImageType::New();
    mitk::CastToItkImage(remove_veins, imop);

    // Extract largest connected component from V2
    ConnectedComponentImageFilterType::Pointer conn1 = ConnectedComponentImageFilterType::New();
    conn1->SetInput(imop);
    conn1->Update();

    LabelShapeKeepNObjImgFilterType::Pointer keepObjs = LabelShapeKeepNObjImgFilterType::New();
    keepObjs->SetInput(conn1->GetOutput());
    keepObjs->SetBackgroundValue(0);
    keepObjs->SetNumberOfObjects(1);
    keepObjs->SetAttribute(LabelShapeKeepNObjImgFilterType::LabelObjectType::NUMBER_OF_PIXELS);
    keepObjs->Update();

    // Iterate to reattach where V==3
    // mitk::Image::Pointer vimage = mitk::IOUtil::Load<mitk::Image>(pathPVeinsCropped.toStdString());
    // ImageType::Pointer veins = ImageType::New();
    mitk::CastToItkImage(vimage, veins);

    ImageType::Pointer outIm = keepObjs->GetOutput();
    IteratorType outImIter(outIm, outIm->GetLargestPossibleRegion());
    IteratorType veinsIter(veins, veins->GetLargestPossibleRegion());

    outImIter.GoToBegin();
    veinsIter.GoToBegin();
    while (!outImIter.IsAtEnd()){
        short value = (veinsIter.Get() == 3) ? 3 : outImIter.Get();
        outImIter.Set(value);

        ++outImIter;
        ++veinsIter;
    }

    mitk::Image::Pointer outputImg = mitk::Image::New();
    mitk::CastToMitkImage(outIm, outputImg);

    return outputImg;
    // end routine
}

void CemrgAtriaClipper::FixClippingErrors(QString imagePath, int strel_radius, bool rewrite){
    QFileInfo fi(imagePath);
    mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(imagePath.toStdString());

    if (image){
        QString prefix = (rewrite) ? "" : "New_";
        QString outputPath = fi.absolutePath() + "/" + prefix + fi.baseName() + "." + fi.completeSuffix();

        mitk::Image::Pointer outputImg = CemrgAtriaClipper::FixClippingErrors(image, strel_radius);
        if (outputImg) {
            // Save image
            mitk::IOUtil::Save(outputImg, outputPath.toStdString());
        } else {
            MITK_WARN << "File not created. Errors in processing.";
        }
    } // if_image
}
