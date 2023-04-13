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
This app serves as a template for the command line apps to be implemented
in the framework.
=========================================================================*/

// Qmitk
#include <mitkIOUtil.h>
#include <mitkSurface.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkCommandLineParser.h>
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
#include <itkResampleImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkGrayscaleErodeImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkImageFileWriter.h>

// Qt
#include <QtDebug>
#include <QString>
#include <QDir>
#include <QFileInfo>
#include <QProcess>
#include <QMessageBox>

// C++ Standard
#include <algorithm>
#include <string>
#include <numeric>

// CemrgApp
#include <CemrgScar3D.h>
#include <CemrgCommandLine.h>

int main(int argc, char* argv[]) {
    mitkCommandLineParser parser;

    // Set general information about your command-line app
    parser.setCategory("Post processing");
    parser.setTitle("Scar Map Projection Options Command-line App");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription("Produce Scar Map shells allowing to choose projection approach.");

    // How should arguments be prefixed
    parser.setArgumentPrefix("--", "-");

    // Add arguments. Unless specified otherwise, each argument is optional.
    // See mitkCommandLineParser::addArgument() for more information.
    // parser.addArgument(
    //   "input-path", "p", mitkCommandLineParser::String,
    //   "Input Directory Path", "Path of directory containing LGE files.",
    //   us::Any(), false);
    parser.addArgument(
        "input-lge", "i", mitkCommandLineParser::String,
        "LGE path", "Full path of LGE.nii file.",
        us::Any(), false);
    parser.addArgument( // optional
        "thresholds-method", "m", mitkCommandLineParser::Int,
        "Thresholds method", "Choose between IIR*V (1) and M + STDV*V (2). (Default=2)");
    parser.addArgument( // optional
        "single-voxel-projection", "svp", mitkCommandLineParser::Bool,
        "Single Voxel Projection", "Project LGE voxels onto Scar Map ONLY ONCE (Default=OFF)");
    parser.addArgument( // optional
        "multi-thresholds", "t", mitkCommandLineParser::Bool,
        "Multiple thresholds", "Produce the output for the scar score using multiple thresholds:\n\t  (mean+V*stdev) V = 1:0.1:5\n\t (V*IIR) V = 0.7:0.01:1.61");
    parser.addArgument( // optional
        "verbose", "v", mitkCommandLineParser::Bool,
        "Verbose Output", "Whether to produce verbose output");

    // Parse arguments.
    // This method returns a mapping of long argument names to their values.
    auto parsedArgs = parser.parseArguments(argc, argv);

    if (parsedArgs.empty())
        return EXIT_FAILURE;

    if (parsedArgs["input-lge"].Empty()) {
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }

    // Parse, cast and set required arguments
    // auto inFilename = us::any_cast<std::string>(parsedArgs["input-path"]);
    MITK_INFO << "Parsing mandatory arguments";
    auto inFilename2 = us::any_cast<std::string>(parsedArgs["input-lge"]);

    // Default values for optional arguments
    // std::string prodthresfile = "prodThresholds.txt";
    auto thresmethod = 2;
    auto singlevoxelprojection = false;
    auto multithreshold = false;
    auto verbose = false;

    // Parse, cast and set optional argument
    MITK_INFO << "Parsing optional arguments";

    if (parsedArgs.end() != parsedArgs.find("thresholds-method")) {
        thresmethod = us::any_cast<int>(parsedArgs["thresholds-method"]);
    }

    if (parsedArgs.end() != parsedArgs.find("multi-thresholds")) {
        multithreshold = us::any_cast<bool>(parsedArgs["multi-thresholds"]);
    }
    std::cout << "multithresh" << multithreshold << '\n';

    if (parsedArgs.end() != parsedArgs.find("single-voxel-projection")) {
        singlevoxelprojection = us::any_cast<bool>(parsedArgs["single-voxel-projection"]);
    }
    std::cout << "single voxel " << singlevoxelprojection << '\n';

    if (parsedArgs.end() != parsedArgs.find("verbose")) {
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);
    }
    std::cout << "verbose" << verbose << '\n';

    MITK_INFO << "Program starting...";
    try {
        // Code the functionality of the cmd app here.
        MITK_INFO(verbose) << "Verbose mode ON.";

        int method = thresmethod;
        QString methodPref = (method == 2) ? "MplusSD_" : "IIR_";
        MITK_INFO(method == 1) << "IIR METHOD";
        MITK_INFO(method == 2) << "M + SD METHOD";

        // PARSING ARGUMENTS
        QString lgename = QString::fromStdString(inFilename2);
        QString segvtk = "segmentation.vtk";
        QString outname = methodPref + "MaxScar";
        outname += singlevoxelprojection ? "-single-voxel" : "-repeated-voxels";
        outname += ".vtk";

        MITK_INFO(verbose) << "Obtaining input file path and working directory: ";

        // OBTAINING directory and lgepath variables
        QFileInfo fi(lgename);
        QString direct = fi.absolutePath();
        QString lgePath = fi.absoluteFilePath();
        QString outputFolder = direct + "/" + fi.baseName() + "_OUTPUT" + "/";

        QDir d(outputFolder);

        if (!d.exists(outputFolder)) {
            MITK_INFO(d.mkpath(outputFolder)) << "Output folder created";
        }

        typedef itk::Image<short, 3> ImageTypeCHAR;
        typedef itk::Image<short, 3> ImageTypeSHRT;

        //Scar projection
        MITK_INFO(verbose) << "Performing Scar projection.";

        int minStep = -1;
        int maxStep = 3;
        int methodType = 2;
        std::unique_ptr<CemrgScar3D> scar(new CemrgScar3D());
        scar->SetMinStep(minStep);
        scar->SetMaxStep(maxStep);
        scar->SetMethodType(methodType);

        MITK_INFO(singlevoxelprojection) << "Setting Single voxel projection";
        MITK_INFO(!singlevoxelprojection) << "Setting multiple voxels projection";
        scar->SetVoxelBasedProjection(singlevoxelprojection);

        ImageTypeCHAR::Pointer segITK = ImageTypeCHAR::New();
        ImageTypeSHRT::Pointer lgeITK = ImageTypeSHRT::New();

        mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>((direct + "/PVeinsCroppedImage.nii").toStdString()), segITK);
        mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>(lgePath.toStdString()), lgeITK);

        itk::ResampleImageFilter<ImageTypeCHAR, ImageTypeCHAR>::Pointer resampleFilter;
        resampleFilter = itk::ResampleImageFilter<ImageTypeCHAR, ImageTypeCHAR>::New();
        resampleFilter->SetInput(segITK);
        resampleFilter->SetReferenceImage(lgeITK);
        resampleFilter->SetUseReferenceImage(true);
        resampleFilter->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<ImageTypeCHAR>::New());
        resampleFilter->SetDefaultPixelValue(0);
        resampleFilter->UpdateLargestPossibleRegion();
        segITK = resampleFilter->GetOutput();
        mitk::IOUtil::Save(mitk::ImportItkImage(segITK), (direct + "/PVeinsCroppedImage.nii").toStdString());
        scar->SetScarSegImage(mitk::ImportItkImage(segITK));

        //Thresholding
        int vxls = 3;
        // int threshType = 1;

        typedef itk::Image<float, 3> ImageType;
        typedef itk::BinaryBallStructuringElement<ImageTypeCHAR::PixelType, 3> BallType;
        typedef itk::GrayscaleErodeImageFilter<ImageTypeCHAR, ImageType, BallType> ErosionFilterType;

        BallType binaryBall;
        binaryBall.SetRadius(vxls);
        binaryBall.CreateStructuringElement();
        ErosionFilterType::Pointer erosionFilter = ErosionFilterType::New();
        erosionFilter->SetInput(segITK);
        erosionFilter->SetKernel(binaryBall);
        erosionFilter->UpdateLargestPossibleRegion();
        mitk::Image::Pointer roiImage = mitk::Image::New();
        roiImage = mitk::ImportItkImage(erosionFilter->GetOutput())->Clone();

        ImageType::Pointer lgeFloat = ImageType::New();
        mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>(lgePath.toStdString()), lgeFloat);

        double mean = 0.0, stdv = 0.0;
        scar->CalculateMeanStd(mitk::ImportItkImage(lgeFloat), roiImage, mean, stdv);

        MITK_INFO(verbose) << "Performing Scar projection using " + segvtk.toStdString();

        QString prodPath = direct + "/";
        mitk::Surface::Pointer scarShell = scar->Scar3D(direct.toStdString(), mitk::ImportItkImage(lgeITK));

        MITK_INFO(verbose) << "Saving new scar map to " + outname.toStdString();

        mitk::IOUtil::Save(scarShell, (outputFolder + outname).toStdString());
        scar->SaveNormalisedScalars(mean, scarShell, outputFolder + "Normalised_" + outname);

        QFileInfo fi2(outputFolder + outname);
        QString prothresfile = fi2.baseName() + "_prodStats.txt";

        MITK_INFO(verbose) << "Writing to pordStats file" + prothresfile.toStdString();
        ofstream prodFile1;
        prodFile1.open((outputFolder + prothresfile).toStdString());
        prodFile1 << methodPref.toStdString() << std::endl;
        prodFile1 << mean << std::endl;
        prodFile1 << stdv << std::endl;

        if (multithreshold) {
            MITK_INFO << "Scores for multiple thresholds.";
            prodFile1 << "MULTIPLE SCORES:" << std::endl;
            double startVal, endVal, increment;

            startVal = (method == 2) ? 1.0 : 0.7;
            endVal = (method == 2) ? 5.0 : 1.61;
            increment = (method == 2) ? 0.1 : 0.01;

            double thisVal = startVal;
            while (thisVal <= endVal) {
                double thisthres = (method == 2) ? mean + thisVal * stdv : mean * thisVal;
                double thispercentage = scar->Thresholding(thisthres);
                prodFile1 << "V=" << thisVal << ", SCORE=" << thispercentage << std::endl;

                thisVal += increment;
            }
        }

        prodFile1.close();


        MITK_INFO(verbose) << "Goodbye!";
    } catch (const std::exception &e) {
        MITK_ERROR << e.what();
        return EXIT_FAILURE;
    } catch (...) {
        MITK_ERROR << "Unexpected error";
        return EXIT_FAILURE;
    }
}
