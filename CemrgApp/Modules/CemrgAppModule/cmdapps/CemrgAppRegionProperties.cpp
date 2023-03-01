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
#include <mitkCommandLineParser.h>
#include <mitkImageCast.h>
#include <mitkIOUtil.h>

// ITK
#include <itkImage.h>
#include <itkPoint.h>
#include <itkImageRegionIterator.h>
#include <itkLabelGeometryImageFilter.h>

// VTK 
#include <vtkPolyData.h>
#include <vtkCenterOfMass.h>

// Qt
#include <QString>
#include <QStringList>
#include <QFileInfo>
#include <QJsonObject>
#include <QJsonValue>
#include <QJsonArray>

// C++ Standard
#include <algorithm>
#include <string>

// CemrgApp
#include <CemrgAtriaClipper.h>
#include <CemrgCommonUtils.h>

int main(int argc, char *argv[]) {
    mitkCommandLineParser parser;

    // Set general information about your command-line app
    parser.setCategory("Tests");
    parser.setTitle("CemrgApp Region Properties Command-line App");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription("Calculate region properties of a binary image or surface mesh.");

    // How should arguments be prefixed
    parser.setArgumentPrefix("--", "-");

    // Add arguments. Unless specified otherwise, each argument is optional.
    // See mitkCommandLineParser::addArgument() for more information.
    parser.addArgument(
        "input", "i", mitkCommandLineParser::String,
        "Input Image", "Any image format known to MITK.",
        us::Any(), false);
    parser.addArgument(
        "output", "o", mitkCommandLineParser::String,
        "Output filename", "Name of output text file (default='input_name_regionprops.txt').",
        us::Any(), false);
    parser.addArgument( // optional
        "vtk-surface", "vtk", mitkCommandLineParser::Bool,
        "VTK Input", "Whether input is a VTK surface mesh");
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

    // Default values for optional arguments
    std::string outFilename = "";
    auto is_vtk = false;
    auto verbose = false;

    // Parse, cast and set optional arguments
    if (parsedArgs.end() != parsedArgs.find("output")) {
        outFilename = us::any_cast<std::string>(parsedArgs["output"]);
    }
    if (parsedArgs.end() != parsedArgs.find("vtk-surface")) {
        is_vtk = us::any_cast<bool>(parsedArgs["vtk-surface"]);
    }

    if (parsedArgs.end() != parsedArgs.find("verbose")) {
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);
    }


    try {
        // Code the functionality of the cmd app here.
        MITK_INFO(verbose) << "Verbose mode ON.";

        MITK_INFO << "The input filename:" << inFilename;
        MITK_INFO << "The output filename:" << outFilename;
        MITK_INFO(is_vtk) << "Processing a VTK file";

        QFileInfo fi(QString::fromStdString(inFilename));
        QString directory = fi.absolutePath();
        QString fname = fi.baseName();
        QString extension = fi.completeSuffix();

        QString arg_in_output = QString::fromStdString(outFilename);
        QString outname = (arg_in_output == "") ? fname : arg_in_output;

        QStringList accepted_formats = is_vtk ? QStringList() << "vtk" << "vtp" << "vtu" : QStringList() << "nii" << "nrrd" << "dcm";
        QString error_msg = is_vtk ? "VTK surface" : "Image";

        if (!accepted_formats.contains(extension)){
            MITK_ERROR << ("Cannot process file <" + extension + "> as " + error_msg).toStdString();
            return EXIT_FAILURE;
        }

        QStringList keys = QStringList()<< "origin" << "bounding_box" << "centroid";
        QStringList values = QStringList()<< "0.0,0.0,0.0" << "0.0,0.0,0.0,0.0,0.0,0.0" << "0.0,0.0,0.0";
        QStringList types = QStringList()<< "array" << "array" << "array";

        QJsonObject json = CemrgCommonUtils::CreateJSONObject(keys, values, types);
        CemrgCommonUtils::WriteJSONFile(json, directory, outname);

        std::string origin_str = "";
        std::string bb_str = "";
        std::string centroid_str = "";

        double bounds[6];
        double centre[3];
        double origin[3];

        if (is_vtk) {
            
            mitk::Surface::Pointer surf = mitk::IOUtil::Load<mitk::Surface>(inFilename);
            vtkSmartPointer<vtkPolyData> pd = surf->Clone()->GetVtkPolyData();

            pd->GetBounds(bounds);

            vtkSmartPointer<vtkCenterOfMass> com = vtkSmartPointer<vtkCenterOfMass>::New();
            com->SetInputData(pd);
            com->SetUseScalarsAsWeights(false);
            com->Update();

            com->GetCenter(centre);

            for (int i = 0; i < 3; i++) {
                origin[i] = bounds[2*i];
            }
        }
        else
        {
            using ImageType = itk::Image<short, 3>;
            using LabelGeometryFilterType = itk::LabelGeometryImageFilter<ImageType>;

            mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(inFilename);
            mitk::Image::Pointer bw_image = CemrgCommonUtils::ReturnBinarised(image);

            image->GetGeometry()->GetOrigin().ToArray(origin);

            ImageType::Pointer itk_image;
            mitk::CastToItkImage(bw_image, itk_image);

            double origin_rm[3] = {0.0, 0.0, 0.0};
            itk_image->SetOrigin(origin_rm);

            LabelGeometryFilterType::Pointer regionprops = LabelGeometryFilterType::New();
            regionprops->SetInput(itk_image);
            regionprops->Update();

            LabelGeometryFilterType::LabelsType::iterator it_labels;
            MITK_INFO(verbose) << ("Number of labels: " + QString::number(regionprops->GetNumberOfLabels())).toStdString();

            for (int i = 0; i < 6; i++) {
                bounds[i] = (double)regionprops->GetBoundingBox(0)[i];
            }

            for (int i = 0; i < 3; i++) {
                centre[i] = (double)regionprops->GetCentroid(0)[i];
            }

        }
        for (int i = 0; i < (int) sizeof(origin); i++) {
            origin_str += std::to_string(origin[i]);
            origin_str += (i < (int) sizeof(origin) - 1) ? "," : "";
        }

        for (int i = 0; i < (int) sizeof(bounds); i++) {
            bb_str += std::to_string(bounds[i]);
            bb_str += (i < (int) sizeof(bounds) - 1) ? "," : "";
        }

        for (int i = 0; i < (int) sizeof(centre); i++) {
            centroid_str += std::to_string(centre[i]);
            centroid_str += (i < (int) sizeof(centre) - 1) ? "," : "";
        }

        MITK_INFO(CemrgCommonUtils::ModifyJSONFile(directory, outname, keys.at(0), QString::fromStdString(origin_str), "array")) << "File updated: origin";
        MITK_INFO(CemrgCommonUtils::ModifyJSONFile(directory, outname, keys.at(1), QString::fromStdString(bb_str), "array")) << "File updated: bb";
        MITK_INFO(CemrgCommonUtils::ModifyJSONFile(directory, outname, keys.at(2), QString::fromStdString(centroid_str), "array")) << "File updated: centroid";
        

        MITK_INFO(verbose) << "Goodbye!";
    } catch (const std::exception &e) {
        MITK_ERROR << e.what();
        return EXIT_FAILURE;
    } catch (...) {
        MITK_ERROR << "Unexpected error";
        return EXIT_FAILURE;
    }
}