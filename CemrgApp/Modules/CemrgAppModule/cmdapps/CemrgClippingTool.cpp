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

int main(int argc, char* argv[]) {
    mitkCommandLineParser parser;

    // Set general information about your command-line app
    parser.setCategory("Scar processing");
    parser.setTitle("Clip Clipping Tool Command-line App");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription(
        "Clip Mitral Valve (or whatever) from any shape (using implicit functions).");

    // How should arguments be prefixed
    parser.setArgumentPrefix("--", "-");

    // Add arguments. Unless specified otherwise, each argument is optional.
    // See mitkCommandLineParser::addArgument() for more information.
    // parser.addArgument(
    //   "input-path", "p", mitkCommandLineParser::String,
    //   "Input Directory Path", "Path of directory containing LGE files.",
    //   us::Any(), false);
    parser.addArgument(
        "input-vtk", "i", mitkCommandLineParser::String,
        "segmentation (vtk) path", "Full path of segmentation.vtk file.",
        us::Any(), false);
    parser.addArgument(
        "output", "o", mitkCommandLineParser::String,
        "Output file", "Where to save the output.",
        us::Any(), false);
    parser.addArgument(
        "clipper", "c", mitkCommandLineParser::String,
        "Mitral valve file (.nii/.vtk)", "Image (or VTK) of mitral valve to be cut");
    parser.addArgument( // optional
        "verbose", "v", mitkCommandLineParser::Bool,
        "Verbose Output", "Whether to produce verbose output");

    // Parse arguments.
    // This method returns a mapping of long argument names to their values.
    auto parsedArgs = parser.parseArguments(argc, argv);

    if (parsedArgs.empty()) {
        return EXIT_FAILURE;
    }
    if (parsedArgs["input-vtk"].Empty() || parsedArgs["output"].Empty() || parsedArgs["clipper"].Empty()) {
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }

    // Parse, cast and set required arguments
    auto inFilename = us::any_cast<std::string>(parsedArgs["input-vtk"]);
    auto outFilename = us::any_cast<std::string>(parsedArgs["output"]);
    std::string clipperfile = us::any_cast<std::string>(parsedArgs["clipper"]);

    // Default values for optional arguments
    auto verbose = false;

    // Parse, cast and set optional arguments
    if (parsedArgs.end() != parsedArgs.find("verbose")) {
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);
    }

    try {
        // Code the functionality of the cmd app here.
        MITK_INFO(verbose) << "Verbose mode ON.";

        // PARSING ARGUMENTS
        QString inname = QString::fromStdString(inFilename);
        QString outname = QString::fromStdString(outFilename);
        QString clipname = QString::fromStdString(clipperfile);

        if (!outname.contains(".vtk", Qt::CaseSensitive))
            outname = outname + ".vtk";

        int whichImplicitFunction;
        if (clipname.contains(".nii", Qt::CaseSensitive))
            whichImplicitFunction = 1;
        else
            whichImplicitFunction = 2;

        MITK_INFO(verbose) << "Obtaining input file path and working directory: ";

        // OBTAINING directory and inputPath variables
        QFileInfo fi(inname);
        QFileInfo fi2(clipname);
        QString direct = fi.absolutePath();
        QString inputPath = fi.absoluteFilePath();
        QString clipPath = fi2.absoluteFilePath();
        QString outputPath = direct + "/" + outname;

        MITK_INFO << ("INPUT: " + inputPath).toStdString();
        MITK_INFO << ("OUTPUT: " + outputPath).toStdString();
        MITK_INFO << ("CLIPPER: " + clipPath).toStdString();

        MITK_INFO(verbose) << "Loading Shell.";
        mitk::Surface::Pointer shell = mitk::IOUtil::Load<mitk::Surface>(inputPath.toStdString());
        vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();

        MITK_INFO(verbose) << "Creating implicit function.";

        if (whichImplicitFunction == 1) {
            MITK_INFO(verbose) << "Loading Clipper image.";
            mitk::Image::Pointer ClipperImage = mitk::IOUtil::Load<mitk::Image>(clipPath.toStdString());
            vtkSmartPointer<vtkImplicitVolume> implicitFn = vtkSmartPointer<vtkImplicitVolume>::New();
            implicitFn->SetVolume(ClipperImage->GetVtkImageData());
            implicitFn->SetOutValue(0.5);
            vtkMTimeType mtime = implicitFn->GetMTime();
            MITK_INFO(verbose) << ("[...] MTime:" + QString::number(mtime)).toStdString();

            MITK_INFO(verbose) << "Creating ClipPolyData object.";
            clipper->SetClipFunction(implicitFn);
        } else {
            MITK_INFO(verbose) << "Loading Clipper surface.";
            mitk::Surface::Pointer ClipperSurface = mitk::IOUtil::Load<mitk::Surface>(clipPath.toStdString());
            vtkSmartPointer<vtkImplicitPolyDataDistance> implicitFn = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
            implicitFn->SetInput(ClipperSurface->GetVtkPolyData());
            // implicitFn->SetTolerance(0.0001);
            vtkMTimeType mtime = implicitFn->GetMTime();
            MITK_INFO(verbose) << ("[...] MTime:" + QString::number(mtime)).toStdString();

            MITK_INFO(verbose) << "Creating ClipPolyData object.";
            clipper->SetClipFunction(implicitFn);

        }

        clipper->SetInputData(shell->GetVtkPolyData());
        clipper->InsideOutOff();
        clipper->Update();

        if (verbose) {
            MITK_INFO << "[DEBUG] Preliminary output generation.";
            QString vPath = direct + "/prelim.vtk";
            shell->SetVtkPolyData(clipper->GetOutput());
            mitk::IOUtil::Save(shell, vPath.toStdString());
        }

        MITK_INFO(verbose) << "Extract and clean surface mesh.";
        vtkSmartPointer<vtkDataSetSurfaceFilter> surfer = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
        surfer->SetInputData(clipper->GetOutput());
        surfer->Update();


        MITK_INFO(verbose) << "[...] Cleaning...";
        vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
        cleaner->SetInputConnection(surfer->GetOutputPort());
        cleaner->Update();

        MITK_INFO(verbose) << "[...] Largest region...";
        vtkSmartPointer<vtkPolyDataConnectivityFilter> lrgRegion = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
        lrgRegion->SetInputConnection(cleaner->GetOutputPort());
        lrgRegion->SetExtractionModeToLargestRegion();
        lrgRegion->Update();

        MITK_INFO(verbose) << "[...] Cleaning a bit more...";
        cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
        cleaner->SetInputConnection(lrgRegion->GetOutputPort());
        cleaner->Update();

        MITK_INFO(verbose) << ("Saving to file: " + outputPath).toStdString();
        shell->SetVtkPolyData(cleaner->GetOutput());
        mitk::IOUtil::Save(shell, outputPath.toStdString());

        MITK_INFO(verbose) << "Goodbye!";
    } catch (const std::exception &e) {
        MITK_ERROR << e.what();
        return EXIT_FAILURE;
    } catch (...) {
        MITK_ERROR << "Unexpected error";
        return EXIT_FAILURE;
    }
}
