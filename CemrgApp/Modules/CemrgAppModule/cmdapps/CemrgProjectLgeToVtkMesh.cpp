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
CEMRG CMD APP TEMPLATE
This app Projects the LGE score onto a surface mesh
in the framework.
=========================================================================*/

// ITK
#include <itkPoint.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>

// VTK
#include <vtkPolyDataWriter.h>
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

// Qmitk
#include <mitkSurface.h>
#include <mitkIOUtil.h>
#include <mitkImageCast.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImage.h>
#include <mitkPointSet.h>
#include <MitkCemrgAppModuleExports.h>
#include <mitkCommandLineParser.h>
#include <mitkIOUtil.h>

// Qt
#include <QString>
#include <QFileInfo>
#include <QProcess>

// CemrgApp
#include <CemrgScar3D.h>
#include <CemrgAtriaClipper.h>
#include <CemrgCommandLine.h>

// C++ Standard
#include <algorithm>
#include <string>

typedef itk::Image<short, 3> itkImageType;
void ItkDeepCopy(itkImageType::Pointer input, itkImageType::Pointer output);
mitk::Surface::Pointer ReadVTKMesh(std::string meshPath);
double GetIntensityAlongNormal(itkImageType::Pointer scarImage, itkImageType::Pointer visitedImage,
                               double n_x, double n_y, double n_z, double centre_x, double centre_y, double centre_z, int minStep = -3, int maxStep = 3);

double GetStatisticalMeasure(itkImageType::Pointer scarSegImage, std::vector<mitk::Point3D> pointsOnAndAroundNormal,
                             itkImageType::Pointer scarImage, itkImageType::Pointer visitedImage, int measure);


int main(int argc, char* argv[]) {

    mitkCommandLineParser parser;

    // Set general information about your command-line app
    parser.setCategory("LGE projection");
    parser.setTitle("LGE score projection Command-line App");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription("Projects a lge score (file .nii) onto a vtk mesh.");

    // How should arguments be prefixed
    parser.setArgumentPrefix("--", "-");

    // Add arguments. Unless specified otherwise, each argument is optional.
    // See mitkCommandLineParser::addArgument() for more information.
    parser.addArgument(
        "input-lge", "lge", mitkCommandLineParser::String,
        "LGE Image", "Full path to the .nii file with the lge score.",
        us::Any(), false);
    parser.addArgument(
        "input-surface", "surf", mitkCommandLineParser::String,
        "Surface Image", "Full path to the .vtk file with the surface.",
        us::Any(), false);
    parser.addArgument(
        "output", "o", mitkCommandLineParser::String,
        "Output file", "Where to save the output.",
        us::Any(), false);
    parser.addArgument(
        "min-step", "minS", mitkCommandLineParser::Int,
        "number of voxels", "Number of voxels towards the interiror to project LGE. Default=1",
        1, true);
    parser.addArgument(
        "max-step", "maxS", mitkCommandLineParser::Int,
        "number of voxels", "Number of voxels towards the exterior to project LGE. Default=3",
        3, true);

    parser.addArgument( // optional
        "verbose", "v", mitkCommandLineParser::Bool,
        "Verbose Output", "Whether to produce verbose output");

    // Parse arguments.
    // This method returns a mapping of long argument names to their values.
    auto parsedArgs = parser.parseArguments(argc, argv);

    if (parsedArgs.empty())
        return EXIT_FAILURE;

    // check for mandatory arguments
    if (parsedArgs["input-lge"].Empty() || parsedArgs["input-surface"].Empty() || parsedArgs["output"].Empty()) {
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }

    // Parse, cast and set required arguments
    auto lgeFilename = us::any_cast<std::string>(parsedArgs["input-lge"]);
    auto surfFilename = us::any_cast<std::string>(parsedArgs["input-surface"]);
    auto outFilename = us::any_cast<std::string>(parsedArgs["output"]);
    // Default values for optional arguments
    auto verbose = false;

    auto minStep = -1 * 3;
    auto maxStep = 3;

    if (parsedArgs.end() != parsedArgs.find("min-step")) {
        minStep = us::any_cast<int>(parsedArgs["min-step"]);
    }

    if (parsedArgs.end() != parsedArgs.find("max-step")) {
        maxStep = us::any_cast<int>(parsedArgs["max-step"]);
    }

    //min step is negative within the code
    if (minStep > 0) {
        minStep = -1 * minStep;
    }

    //max step is positive within the code
    if (maxStep < 0) {
        maxStep = -1 * maxStep;
    }


    // Parse, cast and set optional arguments
    if (parsedArgs.end() != parsedArgs.find("verbose")) {
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);
    }

    try {
        // Code the functionality of the cmd app here.
        MITK_INFO(verbose) << "Verbose mode ON.";
        MITK_INFO << "The lge input filename:" << lgeFilename;
        MITK_INFO << "The surface input filename:" << surfFilename;
        MITK_INFO << "The output filename:" << outFilename;

        // Load the LGE image
        mitk::Image::Pointer lgeImage = mitk::IOUtil::Load<mitk::Image>(lgeFilename);

        //Convert to itk image
        itkImageType::Pointer scarImage;
        mitk::CastToItkImage(lgeImage, scarImage);
        itkImageType::Pointer visitedImage = itkImageType::New();
        ItkDeepCopy(scarImage, visitedImage);

        // Read the surface
        mitk::Surface::Pointer surface = ReadVTKMesh(surfFilename);
        vtkSmartPointer<vtkPolyData> pd = surface->GetVtkPolyData();

        //Calculate normals
        vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
        vtkSmartPointer<vtkPolyData> tempPD = vtkSmartPointer<vtkPolyData>::New();
        tempPD->DeepCopy(pd);
        normals->ComputeCellNormalsOn();
        normals->SetInputData(tempPD);
        normals->SplittingOff();
        normals->Update();
        pd = normals->GetOutput();

        // Declarations
        vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkFloatArray> cellNormals = vtkFloatArray::SafeDownCast(pd->GetCellData()->GetNormals());
        vtkSmartPointer<vtkFloatArray> scalars = vtkSmartPointer<vtkFloatArray>::New();
        itkImageType::IndexType pixelXYZ;
        itkImageType::PointType pointXYZ;
        double pN[3], cP[3];
        double maxSdev = -1e9;
        double maxSratio = -1e9;
        double mean = 0, var = 1;
        double maxScalar = -1;
        double minScalar = 1E10;

        for (int i = 0; i < pd->GetNumberOfCells(); i++) {
            cellNormals->GetTuple(i, pN);
            double cX = 0, cY = 0, cZ = 0, numPoints = 0;
            pd->GetCellPoints(i, cellPoints);
            vtkIdType numCellPoints = cellPoints->GetNumberOfIds();
            for (vtkIdType neighborPoint = 0; neighborPoint < numCellPoints; ++neighborPoint) {
                //Get the neighbor point ID
                vtkIdType neighborPointID = cellPoints->GetId(neighborPoint);

                //Get the neighbor point position
                pd->GetPoint(neighborPointID, cP);

                // ITK method
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
            cX /= numPoints;
            cY /= numPoints;
            cZ /= numPoints;
            // ITK method
            pointXYZ[0] = pN[0];
            pointXYZ[1] = pN[1];
            pointXYZ[2] = pN[2];
            scarImage->TransformPhysicalPointToIndex(pointXYZ, pixelXYZ);
            pN[0] = pixelXYZ[0];
            pN[1] = pixelXYZ[1];
            pN[2] = pixelXYZ[2];
            double scalar = GetIntensityAlongNormal(scarImage, visitedImage, pN[0], pN[1], pN[2], cX, cY, cZ, minStep, maxStep);
            if (scalar > maxScalar) maxScalar = scalar;
            if (scalar < minScalar) minScalar = scalar;
            double sdev = (scalar - mean) / sqrt(var);
            double sratio = scalar / mean;

            if (maxSdev < sdev) maxSdev = sdev;
            if (maxSratio < sratio) maxSratio = sratio;

            /**
         * @brief tickbox GUI for this
         */
            MITK_INFO(verbose) << "int _SCAR_AS_STANDARD_DEVIATION = 1";
            MITK_INFO(verbose) << "int _ONLY_POSITIVE_STDEVS = 1";
            MITK_INFO(verbose) << "int _SCAR_MIP = 1";
            /**
         * @brief end
         */

            // if (_ONLY_POSITIVE_STDEVS == 1 && sdev < 0) sdev = 0;
            if (sdev < 0) sdev = 0;

            //For default scalar to plot
            double scalarToPlot = (scalar - mean) / sqrt(var);

            if (scalarToPlot <= 0) scalarToPlot = 0;

            scalars->InsertTuple1(i, scalarToPlot); // if (_SCAR_MIP == 1 && _SCAR_AS_STANDARD_DEVIATION == 1)
            //else: scalars->InsertTuple1(i, scalar);
        }//_for
        itkImageType::Pointer scarDebugLabel = visitedImage;
        pd->GetCellData()->SetScalars(scalars);
        surface->SetVtkPolyData(pd);

        //Write the surface
        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetInputData(surface->GetVtkPolyData());
        writer->SetFileName(outFilename.c_str());
        writer->Write();

        MITK_INFO(verbose) << "Goodbye!";

    } catch (const std::exception &e) {
        MITK_ERROR << e.what();
        return EXIT_FAILURE;
    } catch (...) {
        MITK_ERROR << "Unexpected error";
        return EXIT_FAILURE;
    }
}

void ItkDeepCopy(itkImageType::Pointer input, itkImageType::Pointer output) {
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


mitk::Surface::Pointer ReadVTKMesh(std::string meshPath) {

    //Load the mesh
    mitk::Surface::Pointer surface = mitk::IOUtil::Load<mitk::Surface>(meshPath);
    vtkSmartPointer<vtkPolyData> pd = surface->GetVtkPolyData();

    //Prepare points for MITK visualisation
    double Xmin = 0, Xmax = 0, Ymin = 0, Ymax = 0, Zmin = 0, Zmax = 0;
    for (int i = 0; i < pd->GetNumberOfPoints(); i++) {
        double* point = pd->GetPoint(i);
        point[0] = -point[0];
        point[1] = -point[1];
        pd->GetPoints()->SetPoint(i, point);
        //Find mins and maxs
        if (i == 0) {
            Xmin = point[0];
            Xmax = point[0];
            Ymin = point[1];
            Ymax = point[1];
            Zmin = point[2];
            Zmax = point[2];
        } else {
            if (point[0] < Xmin) Xmin = point[0];
            if (point[0] > Xmax) Xmax = point[0];
            if (point[1] < Ymin) Ymin = point[1];
            if (point[1] > Ymax) Ymax = point[1];
            if (point[2] < Zmin) Zmin = point[2];
            if (point[2] > Zmax) Zmax = point[2];
        }//_if
    }//_for
    double bounds[6] = {Xmin, Xmax, Ymin, Ymax, Zmin, Zmax};
    surface->GetGeometry()->SetBounds(bounds);

    return surface;
}


double GetIntensityAlongNormal(itkImageType::Pointer scarImage, itkImageType::Pointer visitedImage,
                               double n_x, double n_y, double n_z, double centre_x, double centre_y, double centre_z, int minStep, int maxStep) {

    // Declarations
    int methodType = 2;
    std::vector<mitk::Point3D> pointsOnAndAroundNormal;

    // Normalize
    double tempArr[3];
    tempArr[0] = n_x;
    tempArr[1] = n_y;
    tempArr[2] = n_z;
    double norm = vtkMath::Normalize(tempArr);
    n_x /= norm;
    n_y /= norm;
    n_z /= norm;

    double scar_step_min = minStep;
    double scar_step_max = maxStep;
    double scar_step_size = 1;

    const itkImageType::SizeType sizeOfImage = scarImage->GetLargestPossibleRegion().GetSize();
    double maxX = sizeOfImage[0];
    double maxY = sizeOfImage[1];
    double maxZ = sizeOfImage[2];

    for (double i = scar_step_min; i <= scar_step_max; i += scar_step_size) {
        double x = floor(centre_x + (i * n_x));
        double y = floor(centre_y + (i * n_y));
        double z = floor(centre_z + (i * n_z));
        for (int a = -1; a <= 1; a++) {
            for (int b = -1; b <= 1; b++) {
                for (int c = -1; c <= 1; c++) {
                    if (x + a >= 0 && x + a < maxX && y + b >= 0 && y + b < maxY && z + c >= 0 && z + c < maxZ) {
                        mitk::Point3D tempPoint;
                        tempPoint.SetElement(0, x + a);
                        tempPoint.SetElement(1, y + b);
                        tempPoint.SetElement(2, z + c);
                        pointsOnAndAroundNormal.push_back(tempPoint);
                    }
                }
            }
        }
        //        a=0;b=0;c=0;
        //        if (x+a>=0 && x+a<maxX && y+b>=0 && y+b<maxY && z+c>=0 && z+c<maxZ) {
        //            mitk::Point3D tempPoint;
        //            tempPoint.SetElement(0, x+a);
        //            tempPoint.SetElement(1, y+b);
        //            tempPoint.SetElement(2, z+c);
        //            pointsOnAndAroundNormal.push_back(tempPoint);
        //        }
    }//_for

    double insty = 0;
    insty = GetStatisticalMeasure(scarImage, pointsOnAndAroundNormal, scarImage, visitedImage, methodType);

    return insty;
}


double GetStatisticalMeasure(itkImageType::Pointer scarSegImage, std::vector<mitk::Point3D> pointsOnAndAroundNormal,
                             itkImageType::Pointer scarImage, itkImageType::Pointer visitedImage, int measure) {

    //Declarations
    itkImageType::IndexType pixel_xyz;
    int size = pointsOnAndAroundNormal.size();
    double sum = 0, returnVal = 0; //, visitedStatus;

    //Filter out cut regions
    for (int i = 0; i < size; i++) {
        pixel_xyz[0] = pointsOnAndAroundNormal.at(i).GetElement(0);
        pixel_xyz[1] = pointsOnAndAroundNormal.at(i).GetElement(1);
        pixel_xyz[2] = pointsOnAndAroundNormal.at(i).GetElement(2);
        double val = scarSegImage->GetPixel(pixel_xyz);
        visitedImage->SetPixel(pixel_xyz, 1);
        if (std::abs(val - 3.0) < 1E-10)
            return -1;
    }//_for

    //Return mean
    if (measure == 1) {
        for (int i = 0; i < size; i++) {
            pixel_xyz[0] = pointsOnAndAroundNormal.at(i).GetElement(0);
            pixel_xyz[1] = pointsOnAndAroundNormal.at(i).GetElement(1);
            pixel_xyz[2] = pointsOnAndAroundNormal.at(i).GetElement(2);
            sum += scarImage->GetPixel(pixel_xyz);
        }
        returnVal = sum / size;
    }//_if_mean

    //Return max
    if (measure == 2) {
        int maxIndex = 0;
        double max = -1;
        for (int i = 0; i < size; i++) {
            pixel_xyz[0] = pointsOnAndAroundNormal.at(i).GetElement(0);
            pixel_xyz[1] = pointsOnAndAroundNormal.at(i).GetElement(1);
            pixel_xyz[2] = pointsOnAndAroundNormal.at(i).GetElement(2);
            double greyVal = scarImage->GetPixel(pixel_xyz);
            // visitedStatus = visitedImage->GetPixel(pixel_xyz);
            if (greyVal > max) { // && visitedStatus < 1) {
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
            visitedImage->SetPixel(pixel_xyz, 2);
        }
    }//_if_max

    //Sum along the normal (integration)
    if (measure == 3) {
        for (int i = 0; i < size; i++) {
            pixel_xyz[0] = pointsOnAndAroundNormal.at(i).GetElement(0);
            pixel_xyz[1] = pointsOnAndAroundNormal.at(i).GetElement(1);
            pixel_xyz[2] = pointsOnAndAroundNormal.at(i).GetElement(2);
            double greyVal = scarImage->GetPixel(pixel_xyz);
            sum += greyVal;
        }
        returnVal = sum;
    }//_if_sum

    return returnVal;
}
