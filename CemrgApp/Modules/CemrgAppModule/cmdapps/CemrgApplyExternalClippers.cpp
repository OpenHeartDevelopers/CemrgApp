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
#include <vtkDecimatePro.h>

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
#include <QDir>
#include <QDirIterator>

// C++ Standard
#include <algorithm>
#include <string>
#include <numeric>

// CemrgApp
#include <CemrgScar3D.h>
#include <CemrgAtriaClipper.h>
#include <CemrgCommonUtils.h>

int main(int argc, char* argv[]) {
    mitkCommandLineParser parser;

    // Set general information about your command-line app
    parser.setCategory("Pre processing");
    parser.setTitle("Apply external clippers Command-line App");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription("Apply external clippers.");

    // How should arguments be prefixed
    parser.setArgumentPrefix("--", "-");

    // Add arguments. Unless specified otherwise, each argument is optional.
    parser.addArgument(
        "input", "i", mitkCommandLineParser::File,
        "Segmentation (LA-reg) NIFTI file path", "Full path of .nii file.",
        us::Any(), false);
    parser.addArgument( // optional
        "clip-mv", "mv", mitkCommandLineParser::Bool,
        "Cip Mitral Valve", "Whether to clip mitral valve (default=false)");
    parser.addArgument( // optional
        "verbose", "v", mitkCommandLineParser::Bool,
        "Verbose Output", "Whether to produce verbose output (default=false)");

    // Parse arguments.
    // This method returns a mapping of long argument names to their values.
    auto parsedArgs = parser.parseArguments(argc, argv);

    if (parsedArgs.empty()) {
        return EXIT_FAILURE;
    }

    if (parsedArgs["input"].Empty()) {
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }

    // Parse, cast and set required arguments
    auto inFilename = us::any_cast<std::string>(parsedArgs["input"]);

    // Default values for optional arguments
    auto verbose = false;
    auto clipmv = false;

    // Parse, cast and set optional arguments
    if (parsedArgs.end() != parsedArgs.find("verbose")) {
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);
    }
    if (parsedArgs.end() != parsedArgs.find("clip-mv")) {
        clipmv = us::any_cast<bool>(parsedArgs["clip-mv"]);
    }

    try {
        // Code the functionality of the cmd app here.
        MITK_INFO(verbose) << "Verbose mode ON.";

        // PARSING ARGUMENTS
        QString inputName = QString::fromStdString(inFilename);

        QFileInfo fi(inputName);
        QString base = "prodCutter";
        QString directory = fi.absolutePath();

        MITK_INFO(verbose) << directory.toStdString();

        if(!clipmv){
            MITK_INFO(verbose) << "Looking for prodCutters";
            QDirIterator qiter(directory, {base + "*.vtk", base + "*TNormals.txt"}, QDir::Files);
            QStringList cutfiles, normfiles;
            while (qiter.hasNext())  {
                QString thisFile = qiter.next();
                if (thisFile.contains("TNormals")) {
                    normfiles.push_back(thisFile);
                }
                else {
                    cutfiles.push_back(thisFile);
                }
            }
            cutfiles.sort();
            normfiles.sort();

            mitk::Surface::Pointer shell = mitk::Surface::New();
            QString segvtk = CemrgCommonUtils::GetFilePath(directory, "segmentation", ".vtk");

            if(segvtk.isEmpty()){
                MITK_ERROR << "segmentation.vtk file not found. Create it using MIRTK libraries";
                return EXIT_FAILURE;
            }

            shell = mitk::IOUtil::Load<mitk::Surface>(segvtk.toStdString());
            vtkSmartPointer<vtkDecimatePro> deci = vtkSmartPointer<vtkDecimatePro>::New();
            deci->SetInputData(shell->GetVtkPolyData());
            deci->SetTargetReduction(0.1);
            deci->PreserveTopologyOn();
            deci->Update();
            shell->SetVtkPolyData(deci->GetOutput());

            std::unique_ptr<CemrgAtriaClipper> clipper(new CemrgAtriaClipper(directory, shell));
            mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(inputName.toStdString()); // Load segmentation to transform
            clipper->ClipVeinsImgFromFileList(cutfiles, normfiles, image, "PVeins");
            int strel_radius = 3;
            bool rewriteFile = true;

            QString pveinsCroppedPath = directory + "/PVeinsCroppedImage.nii";
            CemrgAtriaClipper::FixClippingErrors(pveinsCroppedPath, strel_radius, rewriteFile);
        } else {
            MITK_INFO(verbose) << "Looking for Mitral Valve";
            QString output2 = directory + "/segmentation.vtk";
            mitk::Surface::Pointer LAShell = mitk::IOUtil::Load<mitk::Surface>(output2.toStdString());

            QString mviPath = directory + "/prodMVI.vtk";
            MITK_INFO << ("Clipping Mitral valve with file: " + mviPath).toStdString();

            mitk::Surface::Pointer mvclipper = mitk::IOUtil::Load<mitk::Surface>(mviPath.toStdString());
            CemrgCommonUtils::FlipXYPlane(mvclipper, directory, "");
            CemrgCommonUtils::ClipWithPolydata(LAShell, mvclipper, directory + "/segmentation.vtk");
        }

        
        MITK_INFO(verbose) << "Goodbye!";
    } catch (const std::exception &e) {
        MITK_ERROR << e.what();
        return EXIT_FAILURE;
    } catch (...) {
        MITK_ERROR << "Unexpected error";
        return EXIT_FAILURE;
    }
}
