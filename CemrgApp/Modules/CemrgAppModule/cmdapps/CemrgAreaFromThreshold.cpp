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
#include <QStringList>
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
#include <CemrgScarAdvanced.h>
#include <CemrgCommonUtils.h>

QString num2str(double num);

int main(int argc, char* argv[]) {
    mitkCommandLineParser parser;

    // Set general information about your command-line app
    parser.setCategory("Post processing");
    parser.setTitle("Area of Thresholded Surface Mesh Command-line App");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription("Produce surface area values and thresholded meshes.");

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
        "Surface mesh path", "Full path of mesh vtk polydata file.",
        us::Any(), false);
    parser.addArgument(
        "threshold", "t", mitkCommandLineParser::String,
        "Threshold(s)", "Multiple thresholds separated by commas:\n\t  e.g -t 1,11,13,15,17,19");
    parser.addArgument(
        "greater_than_threshold", "geq", mitkCommandLineParser::Bool,
        "Exact Threshold(s)", "Calculate values greater than threshold (Default = false)");
    parser.addArgument( // optional
        "output", "o", mitkCommandLineParser::String,
        "Output filename", "If file exists, appends at the end. (Default: output.csv)");
    parser.addArgument( // optional
        "output_mesh", "omsh", mitkCommandLineParser::String,
        "Thresholded mesh(es)", "If set, several files get output with name <omsh>_<thres>.vtk");
    parser.addArgument( // optional
        "verbose", "v", mitkCommandLineParser::Bool,
        "Verbose Output", "Whether to produce verbose output");

    // Parse arguments.
    // This method returns a mapping of long argument names to their values.
    auto parsedArgs = parser.parseArguments(argc, argv);

    if (parsedArgs.empty())
        return EXIT_FAILURE;

    if (parsedArgs["input"].Empty() || parsedArgs["threshold"].Empty()) {
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }

    // Parse, cast and set required arguments
    // auto inFilename = us::any_cast<std::string>(parsedArgs["input-path"]);
    MITK_INFO << "Parsing mandatory arguments";
    auto inFilename = us::any_cast<std::string>(parsedArgs["input"]);
    auto inThreshold = us::any_cast<std::string>(parsedArgs["threshold"]);

    // Default values for optional arguments
    auto geqThreshold = false;
    std::string outFilename = "output.csv";
    std::string omshFilename = "";
    auto verbose = false;

    // Parse, cast and set optional argument
    MITK_INFO << "Parsing optional arguments";
    if (parsedArgs.end() != parsedArgs.find("greater_than_threshold")) {
        geqThreshold = us::any_cast<bool>(parsedArgs["greater_than_threshold"]);
    }

    if (parsedArgs.end() != parsedArgs.find("output")) {
        outFilename = us::any_cast<std::string>(parsedArgs["output"]);
    }

    if (parsedArgs.end() != parsedArgs.find("output_mesh")) {
        omshFilename = us::any_cast<std::string>(parsedArgs["output_mesh"]);
    }

    if (parsedArgs.end() != parsedArgs.find("verbose")) {
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);
    }
    std::cout << "verbose" << verbose << '\n';

    MITK_INFO << "Program starting...";
    try {
        // Code the functionality of the cmd app here.
        MITK_INFO(verbose) << "Verbose mode ON.";
        // input filename (inFilename)
        QString input_file = QString::fromStdString(inFilename);
        QFileInfo fi(input_file);
        QString direct = fi.absolutePath();

        // threshold(s) (inThreshold)
        QStringList thresList = QString::fromStdString(inThreshold).split(",");
        std::vector<double> thres;
        for (int ix = 0; ix < thresList.size(); ix++) {
            bool ok;
            double th = thresList.at(ix).toDouble(&ok);
            if(ok){
                thres.push_back(th);
            } else{
                MITK_WARN << ("Threshold: " + thresList.at(ix) + " incorrect").toStdString();
            }
        }

        // output name (outFilename)
        QString output_name = direct + "/" + QString::fromStdString(outFilename);
        bool outputFileExists = QFile::exists(output_name);

        // output mesh (omshFilename)
        QString omsh_name = direct + "/" + QString::fromStdString(omshFilename);

        MITK_INFO(verbose) << "Arguments parsed";

        // read in mesh
        mitk::Surface::Pointer surface = mitk::IOUtil::Load<mitk::Surface>(input_file.toStdString());

        // csadv parameters
        mitk::Surface::Pointer tempsurf = mitk::Surface::New();
        tempsurf->SetVtkPolyData(surface->GetVtkPolyData());
        CemrgCommonUtils::SetCellDataToPointData(tempsurf);

        std::unique_ptr<CemrgScarAdvanced> csadv(new CemrgScarAdvanced());

        csadv->SetOutputFileName((direct+"corridor.csv").toStdString());
        csadv->SetOutputPath(direct.toStdString());
        csadv->SetInputData(tempsurf->GetVtkPolyData());
        csadv->SetWeightedCorridorBool(false);
        csadv->SetLeftRightPrefix("");

        std::ofstream prodFile1;
        if (outputFileExists){
            prodFile1.open(output_name.toStdString(), std::ios_base::app);
        }
        else{
            prodFile1.open(output_name.toStdString());
        }

        if(!outputFileExists){
            prodFile1 << "FOLDER,";
            for (unsigned long ix = 0; ix < thres.size(); ix++) {
                prodFile1 << thres.at(ix);
                prodFile1 << ",";
            }
            prodFile1 << '\n';
        }

        prodFile1 << direct.toStdString() << ",";
        for (unsigned long ix = 0; ix < thres.size(); ix++) {
            double th = thres.at(ix);
            double maxval = (geqThreshold) ? 500 : th;

            MITK_INFO << ("Threshold: " + QString::number(th)).toStdString();
            std::string threShellPath = csadv->ThresholdedShell(th);
            mitk::Surface::Pointer threShell = mitk::IOUtil::Load<mitk::Surface>(threShellPath);
            QFile::remove(QString::fromStdString(threShellPath));

            std::string outMeshName = (omsh_name+"_"+num2str(th)+".vtk").toStdString();
            mitk::IOUtil::Save(threShell, outMeshName);

            // area calculation and saving
            csadv->GetSurfaceAreaFromThreshold(th, maxval);
            prodFile1 << csadv->GetLargestSurfaceArea();
            prodFile1 << ",";
        }
        prodFile1 << '\n';

        MITK_INFO(verbose) << "Goodbye!";
    } catch (const std::exception &e) {
        MITK_ERROR << e.what();
        return EXIT_FAILURE;
    } catch (...) {
        MITK_ERROR << "Unexpected error";
        return EXIT_FAILURE;
    }
}

QString num2str(double num){
    QString res = QString::number(num);
    int index = res.lastIndexOf(".");
    if (index>0){
        res.replace(index, 1, "dot");
    }
    return res;
}
