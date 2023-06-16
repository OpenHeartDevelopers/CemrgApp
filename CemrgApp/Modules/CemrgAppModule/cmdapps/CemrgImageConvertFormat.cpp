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
#include <CemrgCommonUtils.h>

int main(int argc, char* argv[]) {
    mitkCommandLineParser parser;

    // Set general information about your command-line app
    parser.setCategory("Image processing");
    parser.setTitle("Convert Image Format Command-line App");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription("Convert between MITK supported image files.");

    // How should arguments be prefixed
    parser.setArgumentPrefix("--", "-");

    // Add arguments. Unless specified otherwise, each argument is optional.
    // See mitkCommandLineParser::addArgument() for more information.
    // parser.addArgument(
    //   "input-path", "p", mitkCommandLineParser::String,
    //   "Input Directory Path", "Path of directory containing LGE files.",
    //   us::Any(), false);
    parser.addArgument(
        "input", "i", mitkCommandLineParser::String,
        "Input image", "Full path of image file.",
        us::Any(), false);
    parser.addArgument(
        "output", "o", mitkCommandLineParser::String,
        "Output filename", "Name of output file. NO EXTENSION (default=test).");
    parser.addArgument(
        "output-format", "ofmt", mitkCommandLineParser::String,
        "Output format (extension)", "Extension of output format (default=nii).");
    parser.addArgument(
        "no-resample", "nr", mitkCommandLineParser::Bool,
        "Do not resample to be isotropic", "Use flag if no resampling is required.");
    parser.addArgument(
        "no-reorient", "nrai", mitkCommandLineParser::Bool,
        "Do not reorient to RAI", "Use flag if no reorienting is required.");
    parser.addArgument(
        "binarise", "b", mitkCommandLineParser::Bool,
        "Converted image is binary", "Use flag if input image is a segmentation (ones and zeros).");
    parser.addArgument( // optional
        "verbose", "v", mitkCommandLineParser::Bool,
        "Verbose Output", "Whether to produce verbose output");

    // Parse arguments.
    // This method returns a mapping of long argument names to their values.
    auto parsedArgs = parser.parseArguments(argc, argv);

    if (parsedArgs.empty())
        return EXIT_FAILURE;

    if (parsedArgs["input"].Empty()) {
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }

    // Parse, cast and set required arguments
    auto inFilename = us::any_cast<std::string>(parsedArgs["input"]);
    // auto outFilename = us::any_cast<std::string>(parsedArgs["output"]);

    // Default values for optional arguments
    std::string outFilename = "test";
    std::string outExt = "nii";
    auto resample = true;
    auto reorient = true;
    auto binarise = false;
    auto verbose = false;

    // Parse, cast and set optional arguments
    if (parsedArgs.end() != parsedArgs.find("output"))
        outFilename = us::any_cast<std::string>(parsedArgs["output"]);

    if (parsedArgs.end() != parsedArgs.find("output-format"))
        outExt = us::any_cast<std::string>(parsedArgs["output-format"]);

    if (parsedArgs.end() != parsedArgs.find("no-resample"))
        resample = !us::any_cast<bool>(parsedArgs["no-resample"]);

    if (parsedArgs.end() != parsedArgs.find("no-reorient"))
        reorient = !us::any_cast<bool>(parsedArgs["no-reorient"]);

    if (parsedArgs.end() != parsedArgs.find("binarise"))
        binarise = us::any_cast<bool>(parsedArgs["binarise"]);

    if (parsedArgs.end() != parsedArgs.find("verbose"))
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);


    try {
        // Code the functionality of the cmd app here.
        MITK_INFO(verbose) << "Verbose mode ON. Input parameters: ";
        MITK_INFO(verbose) << inFilename;
        MITK_INFO(verbose) << outFilename;
        MITK_INFO(verbose) << outExt;
        MITK_INFO(verbose && resample) << "Resample ON";
        MITK_INFO(verbose && !resample) << "Resample OFF";
        MITK_INFO(verbose && reorient) << "Reorient ON";
        MITK_INFO(verbose && !reorient) << "Reorient OFF";
        MITK_INFO(verbose && binarise) << "binarise image ON";
        MITK_INFO(verbose && !binarise) << "binarise image OFF";

        // PARSING ARGUMENTS
        QString inname = QString::fromStdString(inFilename);
        QString outname = QString::fromStdString(outFilename);
        QString ext = QString::fromStdString(outExt);
        if(!ext.startsWith(".")){
            ext = "." + ext;
        }

        MITK_INFO(verbose) << "Obtaining input file path and working directory: ";

        // OBTAINING directory and inputPath variables
        QFileInfo fi(inname);

        QString direct = fi.absolutePath();
        QString inputPath = fi.absoluteFilePath();
        QString outputPath = direct + "/" + outname + ext;

        MITK_INFO << ("INPUT: " + inputPath).toStdString();
        MITK_INFO << ("OUTPUT: " + outputPath).toStdString();

        MITK_INFO(verbose) << "Loading Image.";
        mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(inputPath.toStdString());
        if (!image) {
            MITK_ERROR << "Problem loading image.";
            return EXIT_FAILURE;
        }

        bool success = CemrgCommonUtils::ImageConvertFormat(inputPath, outputPath, resample, reorient, binarise);

        MITK_INFO(success) << "Conversion successful";
        MITK_INFO(verbose) << "Goodbye!";
        
    } catch (const std::exception &e) {
        MITK_ERROR << e.what();
        return EXIT_FAILURE;
    } catch (...) {
        MITK_ERROR << "Unexpected error";
        return EXIT_FAILURE;
    }
}
