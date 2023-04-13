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
#include <CemrgCommonUtils.h>

int main(int argc, char* argv[]) {
    mitkCommandLineParser parser;

    // Set general information about your command-line app
    parser.setCategory("Pre processing");
    parser.setTitle("Resample/Reorient Command-line App");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription("Resample and reorient nifti files.");

    // How should arguments be prefixed
    parser.setArgumentPrefix("--", "-");

    // Add arguments. Unless specified otherwise, each argument is optional.
    parser.addArgument(
        "input", "i", mitkCommandLineParser::String,
        "NIFTI file path", "Full path of .nii file.",
        us::Any(), false);
    parser.addArgument(
        "output", "o", mitkCommandLineParser::String,
        "Output file", "Where to save the output.",
        us::Any(), false);
    parser.addArgument( // optional
        "reorient", "r", mitkCommandLineParser::Bool,
        "Reorient to RAI", "Whether the nifti should be reoriented to RAI configuration (default=true)");
    parser.addArgument( // optional
        "verbose", "v", mitkCommandLineParser::Bool,
        "Verbose Output", "Whether to produce verbose output (default=false)");

    // Parse arguments.
    // This method returns a mapping of long argument names to their values.
    auto parsedArgs = parser.parseArguments(argc, argv);

    if (parsedArgs.empty()) {
        return EXIT_FAILURE;
    }

    if (parsedArgs["input"].Empty() || parsedArgs["output"].Empty()) {
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }

    // Parse, cast and set required arguments
    auto inFilename = us::any_cast<std::string>(parsedArgs["input"]);
    auto outFilename = us::any_cast<std::string>(parsedArgs["output"]);

    // Default values for optional arguments
    auto verbose = false;
    auto reorientToRAI = true;

    // Parse, cast and set optional arguments
    if (parsedArgs.end() != parsedArgs.find("verbose")) {
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);
    }
    if (parsedArgs.end() != parsedArgs.find("reorient")) {
        reorientToRAI = us::any_cast<bool>(parsedArgs["reorient"]);
    }

    try {
        // Code the functionality of the cmd app here.
        MITK_INFO(verbose) << "Verbose mode ON.";

        // PARSING ARGUMENTS
        QString inputName = QString::fromStdString(inFilename);
        QString outname = QString::fromStdString(outFilename);

        if (!outname.contains(".nii", Qt::CaseSensitive)) {
            outname = outname + ".nii";
        }

        MITK_INFO(verbose) << "Obtaining input file path and working directory: ";

        // OBTAINING directory and image path variables
        QFileInfo fi(inputName);
        QString direct = fi.absolutePath();
        QString imagePath = fi.absoluteFilePath();

        MITK_INFO << "Resampling image to isometric.";
        MITK_INFO(verbose && reorientToRAI) << "Reorienting image to RAI.";
        mitk::Image::Pointer image = CemrgCommonUtils::IsoImageResampling(imagePath, reorientToRAI);
        MITK_INFO << "Saving...";
        mitk::IOUtil::Save(image, imagePath.toStdString());

        MITK_INFO(verbose) << "Goodbye!";
    } catch (const std::exception &e) {
        MITK_ERROR << e.what();
        return EXIT_FAILURE;
    } catch (...) {
        MITK_ERROR << "Unexpected error";
        return EXIT_FAILURE;
    }
}
