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
    parser.setTitle("Fix-Shells Command-line App");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription(
        "Fix Scar Map shells with incorrectly assigned zeros.");

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
    parser.addArgument(
        "output", "o", mitkCommandLineParser::String,
        "Output file", "Where to save the output.",
        us::Any(), false);
    parser.addArgument( // optional
        "segmentation-ref", "s", mitkCommandLineParser::String,
        "Segmentation Reference VTK shell", "(Not supported) Segmentation VTK to create ScarMap.");
    parser.addArgument( // optional
        "multi-thresholds", "t", mitkCommandLineParser::Bool,
        "Multiple thresholds", "Produce the output for the scar score using multiple thresholds:\n\t  (mean+V*stdev) V = 1, 2, 2.3, 3.3, 4 and 5\n\t (V*IIR) V = 0.86,0.97, 1.16, 1.2 and 1.32");
    parser.addArgument( // optional
        "verbose", "v", mitkCommandLineParser::Bool,
        "Verbose Output", "Whether to produce verbose output");

    // Parse arguments.
    // This method returns a mapping of long argument names to their values.
    auto parsedArgs = parser.parseArguments(argc, argv);

    if (parsedArgs.empty())
        return EXIT_FAILURE;

    if (//parsedArgs["input-path"].Empty() ||
        parsedArgs["input-lge"].Empty() ||
        parsedArgs["output"].Empty()) {
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }

    // Parse, cast and set required arguments
    // auto inFilename = us::any_cast<std::string>(parsedArgs["input-path"]);
    auto inFilename2 = us::any_cast<std::string>(parsedArgs["input-lge"]);
    auto outFilename = us::any_cast<std::string>(parsedArgs["output"]);

    // Default values for optional arguments
    std::string segref = "segmentation.vtk";
    auto verbose = false;
    auto multithreshold = false;

    // Parse, cast and set optional arguments
    if (parsedArgs.end() != parsedArgs.find("verbose")) {
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);
    }
    if (parsedArgs.end() != parsedArgs.find("segmentation-ref")) {
        segref = us::any_cast<std::string>(parsedArgs["segmentation-ref"]);
    }
    if (parsedArgs.end() != parsedArgs.find("multi-thresholds")) {
        multithreshold = us::any_cast<bool>(parsedArgs["multi-thresholds"]);
    }


    try {
        // Code the functionality of the cmd app here.
        MITK_INFO(verbose) << "Verbose mode ON.";

        // PARSING ARGUMENTS
        // QString direct = QString::fromStdString(inFilename);
        QString lgename = QString::fromStdString(inFilename2);
        QString segvtk = QString::fromStdString(segref);
        QString outname = QString::fromStdString(outFilename);

        if (!outname.contains(".vtk", Qt::CaseSensitive))
            outname = outname + ".vtk";

        if (!segvtk.contains(".vtk", Qt::CaseSensitive))
            segvtk = segvtk + ".vtk";


        MITK_INFO(verbose) << "Obtaining input file path and working directory: ";

        // OBTAINING directory and lgepath variables
        QFileInfo fi(lgename);
        QString direct = fi.absolutePath();
        QString lgePath = fi.absoluteFilePath();
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

        mitk::IOUtil::Save(scarShell, (prodPath + outname).toStdString());
        scar->SaveNormalisedScalars(mean, scarShell, prodPath + "Normalised_" + outname);

        QFileInfo fi2(prodPath + outname);
        QString prothresfile = fi2.baseName() + "_prodStats.txt";

        MITK_INFO(verbose) << "Writing to pordStats file" + prothresfile.toStdString();

        int method;
        double value, thres, percentage;
        double data1[5];
        ifstream prodFileRead;
        QString fileRead = prodPath + "prodThresholds.txt";
        prodFileRead.open(fileRead.toStdString());

        MITK_INFO << "READ FILE: " + fileRead.toStdString();

        if (!verbose)
            for (int i = 0; i < 5; i++)
                prodFileRead >> data1[i];
        else {
            MITK_INFO << "Data READ:";
            for (int i = 0; i < 5; i++) {
                prodFileRead >> data1[i];
                MITK_INFO << data1[i];
            }
        }

        value = data1[0];
        method = data1[1];
        thres = data1[4];

        prodFileRead.close();

        ofstream prodFile1;
        prodFile1.open((prodPath + prothresfile).toStdString());
        prodFile1 << value << std::endl;
        prodFile1 << method << std::endl;
        prodFile1 << mean << std::endl;
        prodFile1 << stdv << std::endl;
        prodFile1 << thres << std::endl;

        percentage = scar->Thresholding(thres);
        prodFile1 << "SCORE:" << percentage << "%" << std::endl;

        if (multithreshold) {
            MITK_INFO << "Scores for multiple thresholds.";
            prodFile1 << "MULTIPLE SCORES:" << std::endl;
            if (method == 2) {
                double manyvalues[6] = {1, 2, 2.3, 3.3, 4, 5};
                for (int i = 0; i < 6; i++) {
                    double thisthres = mean + manyvalues[i] * stdv;
                    double thispercentage = scar->Thresholding(thisthres);
                    prodFile1 << "V = " << manyvalues[i] <<
                        ", SCORE:" << thispercentage << "%" << std::endl;
                }
            } // mean + V*stdv
            else {
                double manyvalues[5] = {0.86, 0.97, 1.16, 1.2, 1.32};
                for (int i = 0; i < 5; i++) {
                    double thisthres = mean * manyvalues[i];
                    double thispercentage = scar->Thresholding(thisthres);
                    prodFile1 << "V = " << manyvalues[i] <<
                        ", SCORE:" << thispercentage << "%" << std::endl;
                }
            } // V*IIR
        }

        prodFile1.close();

        MITK_INFO << "Saving debug scar map labels.";
        scar->SaveScarDebugImage("Max", direct);

        MITK_INFO(verbose) << "Goodbye!";
    } catch (const std::exception &e) {
        MITK_ERROR << e.what();
        return EXIT_FAILURE;
    } catch (...) {
        MITK_ERROR << "Unexpected error";
        return EXIT_FAILURE;
    }
}
