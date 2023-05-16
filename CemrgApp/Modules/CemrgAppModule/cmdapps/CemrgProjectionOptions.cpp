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
    parser.addArgument(
        "input-segmentation", "seg", mitkCommandLineParser::String,
        "Segmentation name", "Name of segmentation file without extension (default=PVeinsCroppedImage)");
    parser.addArgument( // optional
        "output-subfolder", "o", mitkCommandLineParser::String,
        "Output subfolder name", "Name of output subfolder (Default=OUTPUT)");
    parser.addArgument( // optional
        "thresholds-method", "m", mitkCommandLineParser::Int,
        "Thresholds method", "Choose between IIR*V (1) and M + STDV*V (2). (Default=2)");
    parser.addArgument( // optional
        "single-voxel-projection", "svp", mitkCommandLineParser::Bool,
        "Single Voxel Projection", "Project LGE voxels onto Scar Map ONLY ONCE (Default=OFF)");
    parser.addArgument(     // optional
            "threshold-values", "tv", mitkCommandLineParser::String,
            "Threshold(s)", "Specific threshold (-t 1.2):\n Multithresholds (-t 0.97,1.2,1.32)\n Sweeps (-t 0.7:0.01:1.61)");
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
    std::string segFilename = "PVeinsCroppedImage";
    auto method = 2;
    std::string output_subfolder = "OUTPUT";
    auto singlevoxelprojection = false;
    auto multithreshold = false;
    std::string inThresholdString = (method == 2) ? "3.3" : "1.2";
    auto verbose = false;

    // Parse, cast and set optional argument
    MITK_INFO << "Parsing optional arguments";

    if (parsedArgs.end() != parsedArgs.find("input-segmentation")) {
        segFilename = us::any_cast<std::string>(parsedArgs["input-segmentation"]);
    }

    if (parsedArgs.end() != parsedArgs.find("output-subfolder")) {
        output_subfolder = us::any_cast<std::string>(parsedArgs["output-subfolder"]);
    }
    
    if (parsedArgs.end() != parsedArgs.find("thresholds-method")) {
        method = us::any_cast<int>(parsedArgs["thresholds-method"]);
    }

    if (parsedArgs.end() != parsedArgs.find("threshold-values")) {
        inThresholdString = us::any_cast<std::string>(parsedArgs["threshold-values"]);
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

        QString methodPref = (method == 2) ? "MplusSD_" : "IIR_";
        MITK_INFO(method == 1) << "IIR METHOD";
        MITK_INFO(method == 2) << "M + SD METHOD";

        // PARSING ARGUMENTS
        QString lgename = QString::fromStdString(inFilename2);
        QString pveinsname = QString::fromStdString(segFilename);
        QString segvtk = "segmentation.vtk";
        QString outname = methodPref + "MaxScar";
        outname += singlevoxelprojection ? "-single-voxel" : "-repeated-voxels";
        outname += ".vtk";
        pveinsname += (!pveinsname.endsWith(".nii")) ? ".nii" : "";

        MITK_INFO(verbose) << ("Segmentation name: " + pveinsname).toStdString();
        MITK_INFO(verbose) << "Obtaining input file path and working directory: ";

        // OBTAINING directory and lgepath variables
        QFileInfo fi(lgename);
        QString direct = fi.absolutePath();
        QString lgePath = fi.absoluteFilePath();
        QString outputFolder = direct + "/" + QString::fromStdString(output_subfolder) + "/";

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
        int measureType = 2; // Mean=1, Max=2, Cumulative=3, Mode=4
        std::unique_ptr<CemrgScar3D> scar(new CemrgScar3D());
        scar->SetMinStep(minStep);
        scar->SetMaxStep(maxStep);
        scar->SetMethodType(measureType);

        MITK_INFO(singlevoxelprojection) << "Setting Single voxel projection";
        MITK_INFO(!singlevoxelprojection) << "Setting multiple voxels projection";
        scar->SetVoxelBasedProjection(singlevoxelprojection);

        ImageTypeCHAR::Pointer segITK = ImageTypeCHAR::New();
        ImageTypeSHRT::Pointer lgeITK = ImageTypeSHRT::New();

        mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>((direct + "/" + pveinsname).toStdString()), segITK);
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
        mitk::IOUtil::Save(mitk::ImportItkImage(segITK), (direct + "/" + pveinsname).toStdString());
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
        scar->SaveScarDebugImage("Max_debugScar.nii", outputFolder);

        QFileInfo fi2(outputFolder + outname);
        QString prothresfile = fi2.baseName() + "_prodStats.txt";

        MITK_INFO(verbose) << "Writing to pordStats file" + prothresfile.toStdString();
        std::ofstream prodFile1;
        prodFile1.open((outputFolder + prothresfile).toStdString());
        prodFile1 << methodPref.toStdString() << std::endl;
        prodFile1 << mean << std::endl;
        prodFile1 << stdv << std::endl;

        QString thres_str = QString::fromStdString(inThresholdString);
        if (multithreshold && (thres_str.length()==3 && thres_str.contains('.'))) {
            thres_str = (method == 2) ? "1.0:0.1:5.0" : "0.7:0.01:1.61";
        }

        bool is_sweep_threshold = false;
        std::vector<double> input_threshold_vector, threshold_vector;
        QStringList thres_list = QStringList();
        if (thres_str.contains(',')){
            thres_list = thres_str.split(',');
        } else if(thres_str.contains(':')){
            thres_list = thres_str.split(':');
            is_sweep_threshold = true;
        } else{
            thres_list << thres_str;
        }

        int ix = 0;
        bool ok_to_double = true;
        while (ix<thres_list.size() && ok_to_double) {
            input_threshold_vector.push_back(thres_list.at(ix).toDouble(&ok_to_double));
            ix++;
        }

        if (!ok_to_double){
            MITK_ERROR << ("Error parsing value: " + thres_list.at(ix - 1)).toStdString();
            return EXIT_FAILURE;
        }  

        if(is_sweep_threshold) {
            double value = input_threshold_vector.at(0);
            double increment = (input_threshold_vector.size() == 2) ? 1 : input_threshold_vector.at(1);
            double endVal = input_threshold_vector.back();
            while (value < endVal){
                threshold_vector.push_back(value);
                value += increment;
            }
        } else {
            threshold_vector.assign(input_threshold_vector.begin(), input_threshold_vector.end()); 
        }

        for (auto & this_value : threshold_vector)  { 

            double this_thres = (method == 2) ? mean + this_value * stdv : mean * this_value;
            double this_percentage = scar->Thresholding(this_thres);
            prodFile1 << "V=" << this_value << ", SCORE=" << this_percentage << std::endl;
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
