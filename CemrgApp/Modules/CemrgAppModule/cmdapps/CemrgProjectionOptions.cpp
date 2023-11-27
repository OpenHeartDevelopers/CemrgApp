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
#include <QFile>
#include <QProcess>
#include <QMessageBox>
#include <QJsonObject>

// C++ Standard
#include <algorithm>
#include <string>
#include <numeric>

// CemrgApp
#include <CemrgScar3D.h>
#include <CemrgCommandLine.h>
#include <CemrgCommonUtils.h>

// helper functions
QString ParseArgumentsToOutputFolder(std::string output_subfolder, int method, bool svp, bool old, QString limits, QString thresString);
mitk::Image::Pointer Clean(mitk::Image::Pointer segmentation);
QString fromBool(bool value);
void PrintGeometry(mitk::Image::Pointer someImg, std::string name);

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
        "options-file", "opts", mitkCommandLineParser::String,
        "Option(s)", "JSON path with options. Gets copied into output-subfolder ");
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
        "roi-radius", "radius", mitkCommandLineParser::Bool,
        "ROI Radius ON", "Whether to create ROI radius");
    parser.addArgument( // optional
        "roi-legacy-projection", "old", mitkCommandLineParser::Bool,
        "Legacy (old) projection algorithm", "Legacy (old) projection algorithm");
    parser.addArgument( // optional
        "roi-limits", "limits", mitkCommandLineParser::String,
        "Limits(s)", "Limits sent as string with two values separated by a comma (default='-1,3')");
    parser.addArgument( // optional
        "alternative-geometries", "ang", mitkCommandLineParser::Bool,
        "Alternative where geometries are not copied ", "Alternative where geometries are not copied ");
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
    std::string output_subfolder = "";
    std::string options_file = "";
    auto singlevoxelprojection = false;
    std::string inThresholdString = (method == 2) ? "3.3" : "1.2";
    auto roi_radius = false;
    std::string roi_limits = "-1,3";
    auto legacy_projection = false;
    auto ang_flag = false;
    auto verbose = false;

    // Parse, cast and set optional argument
    MITK_INFO << "Parsing optional arguments";

    if (parsedArgs.end() != parsedArgs.find("input-segmentation")) {
        segFilename = us::any_cast<std::string>(parsedArgs["input-segmentation"]);
    }

    if (parsedArgs.end() != parsedArgs.find("output-subfolder")) {
        output_subfolder = us::any_cast<std::string>(parsedArgs["output-subfolder"]);
    }

    if (parsedArgs.end() != parsedArgs.find("options-file")) {
        options_file = us::any_cast<std::string>(parsedArgs["options-file"]);
    }
    
    if (parsedArgs.end() != parsedArgs.find("thresholds-method")) {
        method = us::any_cast<int>(parsedArgs["thresholds-method"]);
    }

    if (parsedArgs.end() != parsedArgs.find("threshold-values")) {
        inThresholdString = us::any_cast<std::string>(parsedArgs["threshold-values"]);
    }

    if (parsedArgs.end() != parsedArgs.find("single-voxel-projection")) {
        singlevoxelprojection = us::any_cast<bool>(parsedArgs["single-voxel-projection"]);
    }
    std::cout << "single voxel " << singlevoxelprojection << '\n';

    if (parsedArgs.end() != parsedArgs.find("roi-radius")) {
        roi_radius = us::any_cast<bool>(parsedArgs["roi-radius"]);
    }

    if (parsedArgs.end() != parsedArgs.find("roi-limits")) {
        roi_limits = us::any_cast<bool>(parsedArgs["roi-limits"]);
    }

    if (parsedArgs.end() != parsedArgs.find("roi-legacy-projection")) {
        legacy_projection = us::any_cast<bool>(parsedArgs["roi-legacy-projection"]);
    }

    if (parsedArgs.end() != parsedArgs.find("alternative-geometries")) {
        ang_flag = us::any_cast<bool>(parsedArgs["alternative-geometries"]);
    }
    if (parsedArgs.end() != parsedArgs.find("alternative-geometries")) {
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);
    }
    std::cout << "verbose" << verbose << '\n';


    MITK_INFO << "Program starting...";
    try {
        
        // Code the functionality of the cmd app here.
        MITK_INFO(verbose) << "Verbose mode ON.";

        QString lgename = QString::fromStdString(inFilename2);
        QFileInfo fi(lgename);
        QString direct = fi.absolutePath();

        // check if options_file is empty
        QString opts = QString::fromStdString(options_file);
        opts += (opts.endsWith(".json")) ? "" : ".json";
        QFileInfo fopts(opts);
        QString optsFile = "options.json";
        if (fopts.exists() && fopts.isFile()) {
            MITK_INFO << "Reading options file. THIS WILL OVERWRITE ANY OTHER ARGUMENTS PASSED IN THE COMMAND LINE.";
            QJsonObject json = CemrgCommonUtils::ReadJSONFile(fopts.absolutePath(), fopts.fileName());
            // output_dir
            if (json.contains("output_dir")){
                output_subfolder = json["output_dir"].toString().toStdString();
                if (verbose) std::cout << "JSON set output subfolder to " << output_subfolder << '\n';
            }
            // thresholds_method
            if (json.contains("thresholds_method")){
                method = json["thresholds_method"].toInt();
                if (verbose) std::cout << "JSON set thresholds method to " << method << '\n';
            }
            // threshold_values
            if (json.contains("threshold_values")){
                inThresholdString = json["threshold_values"].toString().toStdString();
                if (verbose) std::cout << "JSON set threshold values to " << inThresholdString << '\n';
            }
            // single_voxel_projection
            if (json.contains("single_voxel_projection")){
                singlevoxelprojection = json["single_voxel_projection"].toBool();
                if (verbose) std::cout << "JSON set single voxel projection to " << singlevoxelprojection << '\n';
            }
            // roi_radius
            if (json.contains("roi_radius")){
                roi_radius = json["roi_radius"].toBool();
                if (verbose) std::cout << "JSON set roi radius to " << roi_radius << '\n';
            }
            // roi_legacy_projection
            if (json.contains("roi_legacy_projection")){
                legacy_projection = json["roi_legacy_projection"].toBool();
                if (verbose) std::cout << "JSON set roi legacy projection to " << legacy_projection << '\n';
            }
            // roi_limits
            if (json.contains("roi_limits")){
                roi_limits = json["roi_limits"].toString().toStdString();
                if (verbose) std::cout << "JSON set roi limits to " << roi_limits << '\n';
            }

            // optsFile = fopts.fileName();
        }

        roi_radius = roi_radius || legacy_projection;

        QString roiLimits = QString::fromStdString(roi_limits);
        if (!roiLimits.contains(",")) {
            MITK_WARN << "ROI limits not set correctly. Using default values.";
            roiLimits = "-1,3";
        }

        QStringList limits_list = roiLimits.split(",");
        MITK_WARN(limits_list.size() != 2) << ("ROI limits have more than 2 values, using first two: " + limits_list.at(0) + ", " + limits_list.at(1)).toStdString();
        if (limits_list.size() < 2){
            MITK_WARN << "ROI limits not set correctly. Using default values.";
            limits_list.clear();
            limits_list << "-1" << "3";
        }
        
        QString thresString = QString::fromStdString(inThresholdString);

        QString outputSubfolder = ParseArgumentsToOutputFolder(output_subfolder, method, singlevoxelprojection, legacy_projection, roiLimits, thresString);
        QString lgePath = fi.absoluteFilePath();
        QString outputFolder = direct + "/" + outputSubfolder + "/";

        QDir d(outputFolder);

        if (!d.exists(outputFolder)) {
            MITK_INFO(d.mkpath(outputFolder)) << "Output folder created";
        }

        QStringList keys, values, types;
        keys << "output_dir" << "thresholds_method" << "threshold_values" << "single_voxel_projection" << "roi_radius" << "roi_legacy_projection" << "roi_limits";
        values << outputSubfolder << QString::number(method) << thresString << fromBool(singlevoxelprojection) << fromBool(roi_radius) << fromBool(legacy_projection) << roiLimits;
        types << "string" << "int" << "string" << "bool" << "bool" << "bool" << "string";

        QJsonObject json = CemrgCommonUtils::CreateJSONObject(keys, values, types);
        CemrgCommonUtils::WriteJSONFile(json, outputFolder, optsFile);

        QString methodPref = (method == 2) ? "MplusSD_" : "IIR_";
        MITK_INFO(method == 1) << "IIR METHOD";
        MITK_INFO(method == 2) << "M + SD METHOD";

        // PARSING ARGUMENTS
        QString pveinsname = QString::fromStdString(segFilename);
        QString segvtk = "segmentation.vtk";
        QString outname = methodPref + "MaxScar";
        outname += singlevoxelprojection ? "-single-voxel" : "-repeated-voxels";
        outname += ".vtk";
        pveinsname += (!pveinsname.endsWith(".nii")) ? ".nii" : "";

        MITK_INFO(verbose) << ("Segmentation name: " + pveinsname).toStdString();
        MITK_INFO(verbose) << "Obtaining input file path and working directory: ";

        //Scar projection
        MITK_INFO(verbose) << "Performing Scar projection.";
        MITK_INFO(verbose) << "Setting Scar projection parameters.";

        int small = limits_list.at(0).toInt();
        int big = limits_list.at(1).toInt();

        if (small > big) {
            int aux = small;
            small = big;
            big = aux;
        }
    
        int minStep = small;
        int maxStep = big;
        int measureType = 2; // Mean=1, Max=2, Cumulative=3, Mode=4
        std::unique_ptr<CemrgScar3D> scar(new CemrgScar3D());
        scar->SetMinStep(minStep);
        scar->SetMaxStep(maxStep);
        scar->SetMethodType(measureType);
        scar->SetRoiLegacyNormals(legacy_projection);
        scar->SetRoiRadiusOption(roi_radius);

        if (verbose){
            std::cout << "Min step: " << minStep << '\n';
            std::cout << "Max step: " << maxStep << '\n';
            std::cout << "Measure type: " << measureType << '\n';
            std::cout << "ROI legacy normals: " << legacy_projection << '\n';
            std::cout << "ROI radius: " << roi_radius << '\n';
        }

        MITK_INFO(singlevoxelprojection) << "Setting Single voxel projection";
        MITK_INFO(!singlevoxelprojection) << "Setting multiple voxels projection";
        scar->SetVoxelBasedProjection(singlevoxelprojection);

        typedef itk::Image<short, 3> ImageTypeCHAR;
        typedef itk::Image<short, 3> ImageTypeSHRT;

        ImageTypeCHAR::Pointer segITK = ImageTypeCHAR::New();
        ImageTypeSHRT::Pointer lgeITK = ImageTypeSHRT::New();

        mitk::Image::Pointer seg = Clean(mitk::IOUtil::Load<mitk::Image>((direct + "/" + pveinsname).toStdString()));
        mitk::Image::Pointer lge = mitk::IOUtil::Load<mitk::Image>(lgePath.toStdString());

        if (!ang_flag) {
            seg->SetGeometry(lge->GetGeometry());
        }
        // seg->SetGeometry(lge->GetGeometry());
        if (verbose) {
            PrintGeometry(seg, "Segmentation");
            PrintGeometry(lge, "LGE");
        }

        mitk::CastToItkImage(seg, segITK);
        mitk::CastToItkImage(lge, lgeITK);

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
        if (!ang_flag) {
            roiImage->SetGeometry(lge->GetGeometry());
        }
        if (verbose){
            PrintGeometry(roiImage, "ROI");
        }

        mitk::IOUtil::Save(roiImage, (outputFolder + "ROI.nii").toStdString());

        ImageType::Pointer lgeFloat = ImageType::New();
        mitk::CastToItkImage(lge, lgeFloat);

        double mean = 0.0, stdv = 0.0;
        scar->CalculateMeanStd(mitk::ImportItkImage(lgeFloat), roiImage, mean, stdv);

        MITK_INFO(verbose) << "Performing Scar projection using " + segvtk.toStdString();

        QString prodPath = direct + "/";
        mitk::Surface::Pointer scarShell = scar->Scar3D(direct.toStdString(), mitk::ImportItkImage(lgeITK));

        MITK_INFO(verbose) << "Saving new scar map to " + outname.toStdString();

        mitk::IOUtil::Save(scarShell, (outputFolder + outname).toStdString());
        scar->SaveNormalisedScalars(mean, scarShell, outputFolder + "Normalised_" + outname);
        scar->SaveScarDebugImage("DEBUG_" + outname , outputFolder);

        QFileInfo fi2(outputFolder + outname);
        QString prothresfile = fi2.baseName() + "_prodStats.txt";

        MITK_INFO(verbose) << "Writing to pordStats file" + prothresfile.toStdString();
        std::ofstream prodFile1;
        prodFile1.open((outputFolder + prothresfile).toStdString());
        prodFile1 << methodPref.toStdString() << std::endl;
        prodFile1 << mean << std::endl;
        prodFile1 << stdv << std::endl;

        QString thres_str = QString::fromStdString(inThresholdString);

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

QString ParseArgumentsToOutputFolder(std::string output_subfolder, int method, bool svp, bool old, QString limits, QString thresString){
    QString outputSubFolder = QString::fromStdString(output_subfolder);
   
    if (outputSubFolder.isEmpty()) {
        outputSubFolder = "OUTPUT_";
        outputSubFolder += (method == 1) ? "IIR_" : "MSD_";
        outputSubFolder += (old) ? "LEGACY_" : "";
        outputSubFolder += (svp) ? "SVP_" : "";
        outputSubFolder += "ROI";
        outputSubFolder += limits.replace(".", "dot").replace(",", "to").replace("-", "neg");
        outputSubFolder += "_THRES";

        QString auxString = thresString;
        if(auxString.contains(':')){
            QStringList thres_list = auxString.split(':');
            auxString = thres_list.at(0) + "to"+ thres_list.at(2) ;    
        } 
        outputSubFolder += auxString.replace(".", "dot").replace(",", "and").replace("-", "neg");
    }

    return outputSubFolder;
}

mitk::Image::Pointer Clean(mitk::Image::Pointer segmentation){
    using ImageType = itk::Image<float, 3>;
    using IteratorType = itk::ImageRegionIteratorWithIndex<ImageType>;

    ImageType::Pointer im = ImageType::New();
    mitk::CastToItkImage(segmentation, im);

    IteratorType imIter(im, im->GetLargestPossibleRegion());

    imIter.GoToBegin();
    while (!imIter.IsAtEnd()) {
        float value = imIter.Get();

        if (value == 2) {
            value = 3;
        } else if (value > 3){
            value = 0;
        }
        
        imIter.Set(value);
        ++imIter;
    }

    mitk::Image::Pointer outImg = mitk::Image::New();
    mitk::CastToMitkImage(im, outImg);
    return outImg;
}

QString fromBool(bool value){
    return (value) ? "true" : "false";
}

void PrintGeometry(mitk::Image::Pointer someImg, std::string name) {
    std::cout << "Image (" << name << ") Geometry: " << '\n';
    std::cout << "Origin: " << someImg->GetGeometry()->GetOrigin() << '\n';
    std::cout << "Spacing: " << someImg->GetGeometry()->GetSpacing() << '\n';
    std::cout << "Direction: " << '\n';
    std::cout << someImg->GetGeometry()->GetMatrixColumn(0) << '\n';
    std::cout << someImg->GetGeometry()->GetMatrixColumn(1) << '\n';
    std::cout << someImg->GetGeometry()->GetMatrixColumn(2) << '\n';
}
