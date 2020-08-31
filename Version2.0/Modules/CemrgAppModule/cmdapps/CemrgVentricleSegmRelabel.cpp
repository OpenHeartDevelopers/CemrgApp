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
#include <mitkImage.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkCommandLineParser.h>
#include <mitkImagePixelReadAccessor.h>

// VTK
#include <vtkImplicitBoolean.h>
#include <vtkImplicitVolume.h>
#include <vtkImplicitDataSet.h>
#include <vtkClipPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkPlane.h>
#include <vtkSphere.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkMath.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkPolyDataNormals.h>
#include <vtkIdList.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkDataSetSurfaceFilter.h>

// ITK
#include <itkPoint.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkThresholdImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkImageRegionIterator.h>

// Qt
#include <QtDebug>
#include <QString>
#include <QFileInfo>
#include <QProcess>
#include <QMessageBox>
#include <numeric>

#include <CemrgScar3D.h>
#include <CemrgCommandLine.h>

#include <algorithm>
#include <string>

int main(int argc, char* argv[]) {
    mitkCommandLineParser parser;

    // Set general information about your command-line app
    parser.setCategory("EASI processing");
    parser.setTitle("Labelled image re-labelling Command-line App");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription("Relabel segmentation (RV/LV) for Fibres.");

    // How should arguments be prefixed
    parser.setArgumentPrefix("--", "-");

    // Add arguments. Unless specified otherwise, each argument is optional.
    // See mitkCommandLineParser::addArgument() for more information.
    // parser.addArgument(
    //   "input-path", "p", mitkCommandLineParser::InputFile,
    //   "Input Directory Path", "Path of directory containing LGE files.",
    //   us::Any(), false);
    parser.addArgument(
                "input", "i", mitkCommandLineParser::InputFile,
                "Input image segmentation", "Full path of segmentation file.",
                us::Any(), false);
    parser.addArgument(
                "bloodpool", "b", mitkCommandLineParser::InputFile,
                "Input image bloodpool seg", "Full path of bloodpool file.",
                us::Any(), false);
    parser.addArgument( // optional
                "bp-label", "l", mitkCommandLineParser::String,
                "Bloodpool segmentation label", "RV segmentation label (default=30)");
    parser.addArgument( // optional
                "swap-labels", "sl", mitkCommandLineParser::String,
                "Swap levels code", "Code to swap labels (syntax: 'N,M' , default:1,5)");
    parser.addArgument(// optional
                "output", "o", mitkCommandLineParser::String,
                "Output file name", "Where to save the output (default: output.nii).");
    parser.addArgument( // optional
                "debug", "d", mitkCommandLineParser::Bool,
                "Debug Outputs", "Save intermediate steps of processing.");
    parser.addArgument( // optional
                "verbose", "v", mitkCommandLineParser::Bool,
                "Verbose Output", "Whether to produce verbose output");

    // Parse arguments.
    // This method returns a mapping of long argument names to their values.
    auto parsedArgs = parser.parseArguments(argc, argv);

    if (parsedArgs.empty())
        return EXIT_FAILURE;

    if (parsedArgs["input"].Empty() || parsedArgs["bloodpool"].Empty()) {
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }

    // Parse, cast and set required arguments
    auto inFilename = us::any_cast<std::string>(parsedArgs["input"]);
    auto bloodFilename = us::any_cast<std::string>(parsedArgs["bloodpool"]);
    // auto outFilename = us::any_cast<std::string>(parsedArgs["output"]);

    // Default values for optional arguments
    auto verbose = false;
    auto debug = false;
    std::string bp_label = "30";
    std::string swaplabels = "1,5";
    std::string outFilename = "output.nii";

    // Parse, cast and set optional arguments
    if (parsedArgs.end() != parsedArgs.find("verbose")){
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);
    }

    if (parsedArgs.end() != parsedArgs.find("bp-label")){
        swaplabels = us::any_cast<std::string>(parsedArgs["bp-label"]);
    }

    if (parsedArgs.end() != parsedArgs.find("swap-labels")){
        bp_label = us::any_cast<std::string>(parsedArgs["swap-labels"]);
    }

    if (parsedArgs.end() != parsedArgs.find("output")){
        outFilename = us::any_cast<std::string>(parsedArgs["output"]);
    }

    if (parsedArgs.end() != parsedArgs.find("debug")){
        debug = us::any_cast<bool>(parsedArgs["debug"]);
    }


    try{
        // Code the functionality of the cmd app here.
        MITK_INFO(verbose) << "Verbose mode ON.";

        // PARSING ARGUMENTS
        QString inname = QString::fromStdString(inFilename);
        QString bloodpname = QString::fromStdString(bloodFilename);
        QString outname = QString::fromStdString(outFilename);

        if (!outname.contains(".nii", Qt::CaseSensitive))
            outname = outname + ".nii";

        MITK_INFO(verbose) << "Obtaining input file path and working directory: ";

        // OBTAINING directory and inputPath variables
        QFileInfo fi(inname);
        QFileInfo bpfi(bloodpname);

        QString path = fi.absolutePath() + mitk::IOUtil::GetDirectorySeparator();
        QString inputPath = fi.absoluteFilePath();
        QString bpPath = bpfi.absoluteFilePath();
        QString outputPath = path + outname;
        QString debugPrefix = path + "DEBUG_";

        mitk::Point3D origin;

        MITK_INFO << ("INPUT SEGMENTATION: " + inputPath).toStdString();
        MITK_INFO << ("OUTPUT: " + outputPath).toStdString();
        MITK_INFO(debug) << ("DEBUG PREFIX: " + debugPrefix).toStdString();

        MITK_INFO(verbose) << "Loading Image.";
        typedef itk::Image<uint8_t,3> ImageType;
        typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleImageFilterType;
        typedef itk::NearestNeighborInterpolateImageFunction<ImageType> NNInterpolatorType;

        mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(inputPath.toStdString());
        mitk::Image::Pointer bp = mitk::IOUtil::Load<mitk::Image>(bpPath.toStdString());

        if(image && bp){

            MITK_INFO << "Resize bloodpool to image size";
            ImageType::Pointer itkInput = ImageType::New();
            ImageType::Pointer itkBp = ImageType::New();

            mitk::CastToItkImage(image, itkInput);
            mitk::CastToItkImage(bp, itkBp);

            ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
            NNInterpolatorType::Pointer nninterp = NNInterpolatorType::New();

            resampler->SetInput(itkBp);
            resampler->SetInterpolator(nninterp);
            resampler->SetOutputOrigin(itkInput->GetOrigin());
            ImageType::SizeType input_size = itkInput->GetLargestPossibleRegion().GetSize();
            ImageType::SpacingType input_spacing = itkInput->GetSpacing();

            resampler->SetSize(input_size);
            resampler->SetOutputSpacing(input_spacing);
            resampler->SetOutputDirection(itkInput->GetDirection());
            resampler->UpdateLargestPossibleRegion();

            if(debug){
                QString step1 = debugPrefix + "S1_resizedBloodPool.nii";
                mitk::IOUtil::Save(mitk::ImportItkImage(resampler->GetOutput()), step1.toStdString());
            }

            MITK_INFO << "Select ventricle with bloodpool-label";
            // threshold
            typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholdType;
            uint8_t bplabel = uint8_t(std::stoi(bp_label));

            ThresholdType::Pointer thres = ThresholdType::New();
            thres->SetInput(resampler->GetOutput());
            thres->SetLowerThreshold(bplabel);
            thres->SetUpperThreshold(bplabel);
            thres->SetInsideValue(1);
            thres->SetOutsideValue(0);
            thres->Update();

            if(debug){
                QString step2 = debugPrefix + "S2_thresholdedBpRv.nii";
                mitk::IOUtil::Save(mitk::ImportItkImage(thres->GetOutput()), step2.toStdString());
            }

            MITK_INFO << "Dilating selected blooodpool.";
            // Image dilation on thresholded image
            typedef itk::BinaryBallStructuringElement<ImageType::PixelType, 3> StrElType;
            typedef itk::BinaryDilateImageFilter<ImageType, ImageType, StrElType> ImFilterType;

            StrElType structuringElement;
            structuringElement.SetRadius(static_cast<unsigned long>(2.0));
            structuringElement.CreateStructuringElement();

            ImFilterType::Pointer dilateFilter = ImFilterType::New();
            dilateFilter->SetInput(thres->GetOutput());
            dilateFilter->SetKernel(structuringElement);
            dilateFilter->SetDilateValue(1); // same as threshold's insideValue
            // dilateFilter->UpdateLargestPossibleRegion();

            if(debug){
                QString step3 = debugPrefix + "S3_dilatedBpRv.nii";
                mitk::IOUtil::Save(mitk::ImportItkImage(dilateFilter->GetOutput()), step3.toStdString());
            }

            MITK_INFO << "Extracting bloodpool from input image";
            // Multiply input image to thresholed bloodpool (Vfragment)
            typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> ImMultiplyType;
            ImMultiplyType::Pointer imMult = ImMultiplyType::New();
            imMult->SetInput1(itkInput);
            imMult->SetInput2(dilateFilter->GetOutput());

            if (debug){
                QString step4 = debugPrefix + "S4_extractedBloodPool.nii";
                mitk::IOUtil::Save(mitk::ImportItkImage(imMult->GetOutput()), step4.toStdString());
            }

            MITK_INFO << "Relabel specific area in original.";
            // Threshold Vfragment
            // change values on input image where vfragment is 1.
            typedef itk::ImageRegionIterator<ImageType> IteratorType;
            QString swapLabelsStr = QString::fromStdString(swaplabels);
            QStringList swapLabelsList = swapLabelsStr.split(QLatin1Char(','));

            uint8_t fromLabel = uint8_t(swapLabelsList.at(0).toInt());
            uint8_t toLabel = uint8_t(swapLabelsList.at(1).toInt());

            ImageType::Pointer imageFragment = imMult->GetOutput();

            IteratorType relabelIter(itkInput, itkInput->GetLargestPossibleRegion());
            IteratorType condIter(imageFragment, imageFragment->GetLargestPossibleRegion());

            while(!relabelIter.IsAtEnd()){
                if(condIter.Get() == fromLabel){ //label LV=1, change for RV=5
                    relabelIter.Set(toLabel);
                }
                ++relabelIter;
                ++condIter;
            }

            mitk::IOUtil::Save(mitk::ImportItkImage(itkInput), outputPath.toStdString());

        }

        MITK_INFO(verbose) << "Goodbye!";
    }
    catch (const std::exception &e) {
        MITK_ERROR << e.what();
        return EXIT_FAILURE;
    }
    catch(...) {
        MITK_ERROR << "Unexpected error";
        return EXIT_FAILURE;
    }


}
