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
 * Atrial Scar Tools for MITK
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// Qmitk
#include <mitkSurface.h>
#include <mitkIOUtil.h>
#include <mitkImageCast.h>
#include <mitkImagePixelReadAccessor.h>

// VTK
#include <vtkClipPolyData.h>
#include <vtkImplicitBoolean.h>
#include <vtkPlane.h>
#include <vtkSphere.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkMath.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkPolyDataNormals.h>
#include <vtkIdList.h>

// ITK
#include <itkPoint.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>

// Qt
#include <QtDebug>
#include <QMessageBox>
#include <numeric>
#include "CemrgScar3D.h"


CemrgScar3D::CemrgScar3D() {

    this->methodType = 2;
    this->minStep = -3, this->maxStep = 3;
    this->minScalar = 1E10, this->maxScalar = -1;
    this->scalars = vtkSmartPointer<vtkFloatArray>::New();
}

mitk::Surface::Pointer CemrgScar3D::ClipMesh3D(mitk::Surface::Pointer surface, mitk::PointSet::Pointer landmarks) {

    //Retrieve mean and distance of 3 points
    double x_c = 0;
    double y_c = 0;
    double z_c = 0;
    for(int i=0; i<landmarks->GetSize(); i++) {
        x_c = x_c + landmarks->GetPoint(i).GetElement(0);
        y_c = y_c + landmarks->GetPoint(i).GetElement(1);
        z_c = z_c + landmarks->GetPoint(i).GetElement(2);
    }//_for
    x_c /= landmarks->GetSize();
    y_c /= landmarks->GetSize();
    z_c /= landmarks->GetSize();
    double distance[landmarks->GetSize()];
    for(int i=0; i<landmarks->GetSize(); i++) {
        double x_d = landmarks->GetPoint(i).GetElement(0) - x_c;
        double y_d = landmarks->GetPoint(i).GetElement(1) - y_c;
        double z_d = landmarks->GetPoint(i).GetElement(2) - z_c;
        distance[i] = sqrt(pow(x_d,2) + pow(y_d,2) + pow(z_d,2));
    }//_for
    double radius = *std::max_element(distance, distance + landmarks->GetSize());
    double centre[3] = {x_c, y_c, z_c};

    //Clipper
    vtkSmartPointer<vtkSphere> sphere = vtkSmartPointer<vtkSphere>::New();
    sphere->SetCenter(centre);
    sphere->SetRadius(radius);
    vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
    clipper->SetClipFunction(sphere);
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

    //Return the clipped mesh
    surface->SetVtkPolyData(cleaner->GetOutput());
    return surface;
}

mitk::Surface::Pointer CemrgScar3D::Scar3D(std::string directory, mitk::Image::Pointer lgeImage) {

    //Convert to itk image
    itkImageType::Pointer scarImage;
    mitk::CastToItkImage(lgeImage, scarImage);
    itkImageType::Pointer visitedImage = itkImageType::New();
    ItkDeepCopy(scarImage, visitedImage);

    //Read in the mesh
    std::string path = directory + mitk::IOUtil::GetDirectorySeparator() + "segmentation.vtk";
    mitk::Surface::Pointer surface = ReadVTKMesh(path);
    vtkSmartPointer<vtkPolyData> pd = surface->GetVtkPolyData();

    //Calculate normals
    vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
    vtkSmartPointer<vtkPolyData> tempPD = vtkSmartPointer<vtkPolyData>::New();
    tempPD->DeepCopy(pd);
    normals->ComputeCellNormalsOn();
    normals->SetInputData(tempPD);
    normals->Update();
    pd = normals->GetOutput();

    //Declarations
    vtkIdType numCellPoints;
    vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkFloatArray> cellNormals = vtkFloatArray::SafeDownCast(pd->GetCellData()->GetNormals());
    std::vector<double> allScalarsInShell;
    vtkSmartPointer<vtkFloatArray> scalarsOnlyStDev = vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray> scalarsOnlyMultiplier = vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray> scalarsOnlyIntensity = vtkSmartPointer<vtkFloatArray>::New();
    itkImageType::IndexType pixelXYZ;
    itkImageType::PointType pointXYZ;

    double pN[3], cP[3];
    double numPoints = 0;
    double cX = 0, cY = 0, cZ = 0;
    double sdev, maxSdev, sratio, maxSratio = -1e9;
    double scalar = 0;
    double mean = 0, var = 1;

    for (int i=0; i<pd->GetNumberOfCells(); i++) {

        vtkIdType neighborPoint;
        cellNormals->GetTuple(i, pN);
        cX = 0, cY = 0, cZ = 0, numPoints = 0;
        pd->GetCellPoints(i, cellPoints);
        numCellPoints = cellPoints->GetNumberOfIds();

        for (neighborPoint=0; neighborPoint<numCellPoints; ++neighborPoint) {

            //Get the neighbor point ID
            vtkIdType neighborPointID = cellPoints->GetId(neighborPoint);

            //Get the neighbor point position
            pd->GetPoint(neighborPointID, cP);

            //ITK method
            pointXYZ[0] = cP[0];
            pointXYZ[1] = cP[1];
            pointXYZ[2] = cP[2];
            scarImage->TransformPhysicalPointToIndex(pointXYZ, pixelXYZ);
            cP[0] = pixelXYZ[0];
            cP[1] = pixelXYZ[1];
            cP[2] = pixelXYZ[2];
            cX += cP[0];
            cY += cP[1];
            cZ += cP[2];
            numPoints++;

        }//_innerLoop

        //ITK method
        cX /= numPoints;
        cY /= numPoints;
        cZ /= numPoints;
        pointXYZ[0] = pN[0];
        pointXYZ[1] = pN[1];
        pointXYZ[2] = pN[2];
        scarImage->TransformPhysicalPointToIndex(pointXYZ, pixelXYZ);
        pN[0] = pixelXYZ[0];
        pN[1] = pixelXYZ[1];
        pN[2] = pixelXYZ[2];
        scalar = GetIntensityAlongNormal(scarImage, visitedImage, pN[0], pN[1], pN[2], cX, cY, cZ);

        /**
         * @brief tickbox GUI for this
         */
        int _ONLY_POSITIVE_STDEVS = 1;
        int _SCAR_AS_STANDARD_DEVIATION = 1;
        int _SCAR_MIP = 1;
        /**
         * @brief end
         */

        if (scalar > maxScalar) maxScalar = scalar;
        if (scalar < minScalar) minScalar = scalar;
        sdev = (scalar-mean) / sqrt(var);
        sratio = scalar / mean;
        if (maxSdev < sdev) maxSdev = sdev;
        if (maxSratio < sratio) maxSratio = sratio;
        if (_ONLY_POSITIVE_STDEVS == 1 && sdev < 0) sdev = 0;
        scalarsOnlyStDev->InsertTuple1(i, sdev);
        scalarsOnlyIntensity->InsertTuple1(i, scalar);
        scalarsOnlyMultiplier->InsertTuple1(i, sratio);
        //For default scalar to plot
        double scalarToPlot = (scalar-mean) / sqrt(var);

        if (scalarToPlot <= 0)
            scalarToPlot = 0;
        if (_SCAR_MIP == 1 && _SCAR_AS_STANDARD_DEVIATION == 1) {
            scalars->InsertTuple1(i, scalarToPlot);
            allScalarsInShell.push_back(scalarToPlot);
        } else {
            scalars->InsertTuple1(i, scalar);
            allScalarsInShell.push_back(scalar);
        }

    }//_for

    pd->GetCellData()->SetScalars(scalars);
    surface->SetVtkPolyData(pd);
    return surface;
}

bool CemrgScar3D::CalculateMeanStd(mitk::Image::Pointer lgeImage, mitk::Image::Pointer roiImage, double& mean, double& stdv) {

    //Access image volumes
    mitk::ImagePixelReadAccessor<float, 3> readAccess1(lgeImage);
    float* pvLGE = (float*)readAccess1.GetData();
    mitk::ImagePixelReadAccessor<float, 3> readAccess2(roiImage);
    float* pvROI = (float*)readAccess2.GetData();

    int dimsLGE = lgeImage->GetDimensions()[0] * lgeImage->GetDimensions()[1] * lgeImage->GetDimensions()[2];
    int dimsROI = roiImage->GetDimensions()[0] * roiImage->GetDimensions()[1] * roiImage->GetDimensions()[2];
    if (dimsLGE != dimsROI) {
        QMessageBox::critical(NULL, "Attention", "The mask and the image dimensions do not match!");
        return false;
    }//_wrong dimensions

    //Loop image voxels
    std::vector<float> voxelValues;
    for (int i=0; i<dimsROI; i++) {
        if (*pvROI == 1)
            voxelValues.push_back(*pvLGE);
        pvLGE++;
        pvROI++;
    }//_for

    //Calculate mean and std
    double sumDeviation = 0.0;
    double sum = std::accumulate(voxelValues.begin(), voxelValues.end(), 0.0);
    mean = sum / voxelValues.size();
    for (unsigned int i=0; i<voxelValues.size(); i++)
        sumDeviation += (voxelValues[i]-mean) * (voxelValues[i]-mean);
    stdv = std::sqrt(sumDeviation/voxelValues.size());
    return true;
}

double CemrgScar3D::Thresholding(double thresh) {

    double value;
    int ctr1 = 0; int ctr2 = 0;
    for(int i=0; i<scalars->GetNumberOfTuples(); i++) {
        value = scalars->GetValue(i);
        if (value == -1) {
            ctr1++;
            continue;
        }//_if
        if (value > thresh) ctr2++;
    }
    double percentage = (ctr2*100.0) / (scalars->GetNumberOfTuples() - ctr1);
    return percentage;
}

double CemrgScar3D::GetMinScalar() const {

    return minScalar;
}

double CemrgScar3D::GetMaxScalar() const {

    return maxScalar;
}

void CemrgScar3D::SetMinStep(int value) {

    minStep = value;
}

void CemrgScar3D::SetMaxStep(int value) {

    maxStep = value;
}

void CemrgScar3D::SetMethodType(int value) {

    methodType = value;
}

void CemrgScar3D::SetScarSegImage(const mitk::Image::Pointer image) {

    //Setup roiImage
    itkImageType::Pointer itkImage = itkImageType::New();
    mitk::CastToItkImage(image, itkImage);
    this->scarSegImage = itkImage;
}

double CemrgScar3D::GetIntensityAlongNormal(
        itkImageType::Pointer scarImage, itkImageType::Pointer visitedImage,
        double n_x, double n_y, double n_z, double centre_x, double centre_y, double centre_z) {

    //Declarations
    int a, b, c, maxX, maxY, maxZ;
    double insty = 0, x = 0, y = 0, z = 0;
    std::vector<mitk::Point3D> pointsOnAndAroundNormal;

    //Normalize
    double tempArr[3];
    tempArr[0] = n_x;
    tempArr[1] = n_y;
    tempArr[2] = n_z;
    double norm = vtkMath::Normalize(tempArr);
    n_x /= norm;
    n_y /= norm;
    n_z /= norm;

    double scar_step_min  = minStep;
    double scar_step_max  = maxStep;
    double scar_step_size = 1;

    const itkImageType::SizeType sizeOfImage = scarImage->GetLargestPossibleRegion().GetSize();
    maxX = sizeOfImage[0];
    maxY = sizeOfImage[1];
    maxZ = sizeOfImage[2];

    for (double i = scar_step_min; i <= scar_step_max; i += scar_step_size) {

        x = centre_x + (i*n_x);
        y = centre_y + (i*n_y);
        z = centre_z + (i*n_z);
        x = floor(x);
        y = floor(y);
        z = floor(z);

        for (a=-1; a<=1; a++) {
            for (b=-1; b<=1; b++) {
                for (c=-1; c<=1; c++) {
                    if (x+a>=0 && x+a<maxX && y+b>=0 && y+b<maxY && z+c>=0 && z+c<maxZ) {

                        mitk::Point3D tempPoint;
                        tempPoint.SetElement(0, x+a);
                        tempPoint.SetElement(1, y+b);
                        tempPoint.SetElement(2, z+c);
                        pointsOnAndAroundNormal.push_back(tempPoint);
                    }
                }
            }
        }
        //a=0;b=0;c=0;
        //if (x+a>=0 && x+a<maxX && y+b>=0 && y+b<maxY && z+c>=0 && z+c<maxZ) {
        //mitk::Point3D tempPoint;
        //tempPoint.SetElement(0, x+a);
        //tempPoint.SetElement(1, y+b);
        //tempPoint.SetElement(2, z+c);
        //pointsOnAndAroundNormal.push_back(tempPoint);
        //}
    }//_for

    if (methodType == 1) {

        //Statistical measure 1 returns mean
        insty = GetStatisticalMeasure(pointsOnAndAroundNormal, scarImage, visitedImage, 1);

    } else if (methodType == 2) {

        //Statistical measure 2 returns max
        insty = GetStatisticalMeasure(pointsOnAndAroundNormal, scarImage, visitedImage, 2);

    }//_if

    return insty;
}

double CemrgScar3D::GetStatisticalMeasure(
        std::vector<mitk::Point3D> pointsOnAndAroundNormal,
        itkImageType::Pointer scarImage, itkImageType::Pointer visitedImage, int measure) {

    //Declarations
    itkImageType::IndexType pixel_xyz;
    int size = pointsOnAndAroundNormal.size(), maxIndex;
    double sum = 0, max = -1, greyVal, visitedStatus, returnVal;

    //Filter out cut regions
    for (int i=0; i<size; i++) {
        pixel_xyz[0] = pointsOnAndAroundNormal.at(i).GetElement(0);
        pixel_xyz[1] = pointsOnAndAroundNormal.at(i).GetElement(1);
        pixel_xyz[2] = pointsOnAndAroundNormal.at(i).GetElement(2);
        double val = scarSegImage->GetPixel(pixel_xyz);
        if (std::abs(val - 3.0) < 1E-10)
            return -1;
    }//_for

    //Reutrn mean
    if (measure == 1) {

        for (int i=0; i<size; i++) {
            pixel_xyz[0] = pointsOnAndAroundNormal.at(i).GetElement(0);
            pixel_xyz[1] = pointsOnAndAroundNormal.at(i).GetElement(1);
            pixel_xyz[2] = pointsOnAndAroundNormal.at(i).GetElement(2);
            sum += scarImage->GetPixel(pixel_xyz);
        }
        returnVal = sum/size;
    }//_if_mean

    //Return max
    if (measure == 2) {

        for (int i=0; i<size; i++) {
            pixel_xyz[0] = pointsOnAndAroundNormal.at(i).GetElement(0);
            pixel_xyz[1] = pointsOnAndAroundNormal.at(i).GetElement(1);
            pixel_xyz[2] = pointsOnAndAroundNormal.at(i).GetElement(2);
            greyVal = scarImage->GetPixel(pixel_xyz);
            visitedStatus = visitedImage->GetPixel(pixel_xyz);
            if (greyVal > max) {// && visitedStatus < 1) {
                max = greyVal;
                maxIndex = i;
            }
        }
        if (max == -1) {
            returnVal = 0;
        } else {
            returnVal = max;
            //Now change the visited status of this max pixel
            pixel_xyz[0] = pointsOnAndAroundNormal.at(maxIndex).GetElement(0);
            pixel_xyz[1] = pointsOnAndAroundNormal.at(maxIndex).GetElement(1);
            pixel_xyz[2] = pointsOnAndAroundNormal.at(maxIndex).GetElement(2);
            visitedImage->SetPixel(pixel_xyz, 1);
        }
    }//_if_max

    //Sum along the normal (integration)
    if (measure == 3) {

        for (int i=0; i<size; i++) {
            pixel_xyz[0] = pointsOnAndAroundNormal.at(i).GetElement(0);
            pixel_xyz[1] = pointsOnAndAroundNormal.at(i).GetElement(1);
            pixel_xyz[2] = pointsOnAndAroundNormal.at(i).GetElement(2);
            greyVal = scarImage->GetPixel(pixel_xyz);
            sum += greyVal;
        }
        returnVal = sum;
    }//_if_sum

    return returnVal;
}

void CemrgScar3D::ItkDeepCopy(itkImageType::Pointer input, itkImageType::Pointer output) {

    output->SetRegions(input->GetLargestPossibleRegion());
    output->Allocate();

    itk::ImageRegionConstIterator<itkImageType> inputIterator(input, input->GetLargestPossibleRegion());
    itk::ImageRegionIterator<itkImageType> outputIterator(output, output->GetLargestPossibleRegion());

    while (!inputIterator.IsAtEnd()) {
        outputIterator.Set(0);
        ++inputIterator;
        ++outputIterator;
    }
}

mitk::Surface::Pointer CemrgScar3D::ReadVTKMesh(std::string meshPath) {

    //Load the mesh
    mitk::Surface::Pointer surface = mitk::IOUtil::LoadSurface(meshPath);
    vtkSmartPointer<vtkPolyData> pd = surface->GetVtkPolyData();

    //Prepare points for MITK visualisation
    double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;
    for (int i=0; i<pd->GetNumberOfPoints(); i++) {
        double* point = pd->GetPoint(i);
        point[0] = -point[0];
        point[1] = -point[1];
        pd->GetPoints()->SetPoint(i, point);
        //Find mins and maxs
        if (i==0) {
            Xmin = point[0];
            Xmax = point[0];
            Ymin = point[1];
            Ymax = point[1];
            Zmin = point[2];
            Zmax = point[2];
        } else {
            if (point[0]<Xmin) Xmin = point[0];
            if (point[0]>Xmax) Xmax = point[0];
            if (point[1]<Ymin) Ymin = point[1];
            if (point[1]>Ymax) Ymax = point[1];
            if (point[2]<Zmin) Zmin = point[2];
            if (point[2]>Zmax) Zmax = point[2];
        }//_if
    }//_for
    double bounds[6] = {Xmin, Xmax, Ymin, Ymax, Zmin, Zmax};
    surface->GetGeometry()->SetBounds(bounds);

    return surface;
}
