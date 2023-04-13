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
// #include <mitkImageToSurfaceFilter.h>
#include "mitkManualSegmentationToSurfaceFilter.h"

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
#include <itkAddImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkLabelShapeKeepNObjectsImageFilter.h>

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

typedef itk::Image<uint16_t,3> ImageType;
typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholdType;
typedef itk::BinaryBallStructuringElement<ImageType::PixelType, 3> StrElType;
typedef itk::BinaryDilateImageFilter<ImageType, ImageType, StrElType> ImFilterType;
typedef itk::ImageRegionIterator<ImageType> IteratorType;

typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> ConnectedComponentImageFilterType;
typedef itk::LabelShapeKeepNObjectsImageFilter<ImageType> LabelShapeKeepNObjImgFilterType;

ThresholdType::Pointer thresholdImage(ImageType::Pointer input, uint16_t thresholdVal, QString debugOutput, bool debugging);
void OutputImage(ImageType::Pointer im, QString dir, QString imName);

int main(int argc, char* argv[]) {
    mitkCommandLineParser parser;
    // Set general information about your command-line app
    parser.setCategory("Image analysis");
    parser.setTitle("Colour Surface Command-line App");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription(
        "This command line app creates a surface out of a multilabeled image.");

    // How should arguments be prefixed
    parser.setArgumentPrefix("--", "-");

    // Add arguments. Unless specified otherwise, each argument is optional.
    // See mitkCommandLineParser::addArgument() for more information.
    parser.addArgument("input-segmentation", "i",
        mitkCommandLineParser::String,"Input segmentation (.nii)",
        "Input segmentation file (commonly output of CEMRGNET)", us::Any(), false
    );

    parser.addArgument("output", "o",
        mitkCommandLineParser::String, "Output name (.vtk)",
        "Output surface VTK file. "
    );

    parser.addArgument("verbose", "v",
        mitkCommandLineParser::Bool, "Verbose Output",
        "Whether to produce verbose output"
    );

    parser.addArgument("debug", "d",
        mitkCommandLineParser::Bool, "Debug Output",
        "Whether to produce debug output"
    );

    // Parse arguments.
    // This method returns a mapping of long argument names to their values.
    auto parsedArgs = parser.parseArguments(argc, argv);

    if (parsedArgs.empty())
    return EXIT_FAILURE;

    if (parsedArgs["input-segmentation"].Empty()) {
        MITK_INFO << parser.helpText();
        return EXIT_FAILURE;
    }

    // Parse, cast and set required arguments
    auto input1 = us::any_cast<std::string>(parsedArgs["input-segmentation"]);

    // Default values for optional arguments
    auto verbose = false;
    auto debug = false;
    std::string output = "labelledSurface.vtk";

    // Parse, cast and set optional arguments
    if (parsedArgs.end() != parsedArgs.find("verbose")){
        verbose = us::any_cast<bool>(parsedArgs["verbose"]);
    }
    if (parsedArgs.end() != parsedArgs.find("output")){
        output = us::any_cast<bool>(parsedArgs["output"]);
    }
    if (parsedArgs.end() != parsedArgs.find("debug")){
        debug = us::any_cast<bool>(parsedArgs["debug"]);
    }

    try{
        // Code the functionality of the cmd app here.
        MITK_INFO(verbose) << "Verbose mode ON.";

        // PARSING ARGUMENTS
        QString inname = QString::fromStdString(input1);
        QString outname = QString::fromStdString(output);

        if (!outname.contains(".vtk", Qt::CaseSensitive))
            outname = outname + ".vtk";

        MITK_INFO(verbose) << "Obtaining input file path and working directory: ";

        // OBTAINING directory and inputPath variable
        QFileInfo fi(inname);

        QString dir = fi.absolutePath();
        QString path = dir + "/";
        QString inputPath = fi.absoluteFilePath();
        QString outputPath = path + outname;

        int minStep = -1;
        int maxStep = 3;
        int methodType = 2;
        std::unique_ptr<CemrgScar3D> scar(new CemrgScar3D());
        scar->SetMinStep(minStep);
        scar->SetMaxStep(maxStep);
        scar->SetMethodType(methodType);

        if(inname.contains(".vtk", Qt::CaseSensitive)){
            MITK_INFO << "Assume steps 1-4 are finished and stored in folder.";
            QString segPath = path + "1_CleanSegmentation.nii";
            QString relabelledPath = path + "4_RelabelledVeins.nii";

            mitk::Image::Pointer segMITK = mitk::IOUtil::Load<mitk::Image>(relabelledPath.toStdString());
            // scar->SetScarSegImage(mitk::IOUtil::Load<mitk::Image>(segPath.toStdString()));
            scar->SetScarSegImage(segMITK);
            mitk::Surface::Pointer scarShell = scar->Scar3D(dir.toStdString(), segMITK);

            MITK_INFO(verbose) << ("Saving output shell to " + outputPath).toStdString();
            mitk::IOUtil::Save(scarShell, outputPath.toStdString());
        } else{
            MITK_INFO(verbose) << "Loading Image.";
            mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(inputPath.toStdString());
            if(image){
                ImageType::Pointer orgSegImage = ImageType::New();
                mitk::CastToItkImage(image, orgSegImage);

                MITK_INFO(verbose) << "Extracting clean segmentation.";
                ConnectedComponentImageFilterType::Pointer conn1 = ConnectedComponentImageFilterType::New();
                conn1->SetInput(orgSegImage);
                conn1->Update();

                LabelShapeKeepNObjImgFilterType::Pointer keepObjs = LabelShapeKeepNObjImgFilterType::New();
                keepObjs->SetInput(conn1->GetOutput());
                keepObjs->SetBackgroundValue(0);
                keepObjs->SetNumberOfObjects(1);
                keepObjs->SetAttribute(LabelShapeKeepNObjImgFilterType::LabelObjectType::NUMBER_OF_PIXELS);
                keepObjs->Update();
                ImageType::Pointer segImage = ImageType::New();
                segImage = keepObjs->GetOutput();

                OutputImage(segImage, path, "1_CleanSegmentation.nii");

                MITK_INFO(verbose) << "Relabel clean segmentation's veins";
                IteratorType segImgIter(segImage, segImage->GetLargestPossibleRegion());
                IteratorType ogImgIter(orgSegImage, orgSegImage->GetLargestPossibleRegion());
                while(!segImgIter.IsAtEnd()){
                    if(segImgIter.Get() > 0){
                        if(ogImgIter.Get() == 3){ // mitral valve
                            segImgIter.Set(30);
                        } else{
                            segImgIter.Set((int)ogImgIter.Get());
                        }
                    }
                    ++segImgIter;
                    ++ogImgIter;
                }

                OutputImage(segImage, path, "2_SegmentationWithVeins.nii");

                MITK_INFO(verbose) << "Thresholding veins from clean segmentation";
                uint16_t veinsthresh = 2;
                QString thesOutput = path + "2_ThresholdVeins.nii";
                ThresholdType::Pointer thresVeins = thresholdImage(segImage, veinsthresh, thesOutput, debug);

                // possible morphological opening step

                MITK_INFO(verbose) << "Label each detected vein";
                ConnectedComponentImageFilterType::Pointer conn2 = ConnectedComponentImageFilterType::New();
                conn2->SetInput(thresVeins->GetOutput());
                conn2->Update();

                OutputImage(conn2->GetOutput(), path, "3_ConnectedComponentsVeins.nii");

                MITK_INFO(verbose) << "Relabel original image";
                IteratorType segImgIter2(segImage, segImage->GetLargestPossibleRegion());
                IteratorType veinsLabelIter(conn2->GetOutput(), conn2->GetOutput()->GetLargestPossibleRegion());
                while(!veinsLabelIter.IsAtEnd()){
                    if(veinsLabelIter.Get() > 0){ //label LV=1, change for RV=5
                        segImgIter2.Set((int)veinsLabelIter.Get()+(int)segImgIter2.Get());
                    }
                    ++segImgIter2;
                    ++veinsLabelIter;
                }

                OutputImage(segImage, path, "4_RelabelledVeins.nii");

                MITK_INFO(verbose) << "Extract surface";
                double th   = 0.5;
                double bl   = 0.8;
                double smth = 3;
                double ds   = 0.5;

                mitk::Image::Pointer veinsRelabeledImg = mitk::Image::New();
                mitk::CastToMitkImage(segImage, veinsRelabeledImg);
                auto im2surf = mitk::ManualSegmentationToSurfaceFilter::New();

                im2surf->SetInput(veinsRelabeledImg);
                im2surf->SetThreshold(th);
                im2surf->SetUseGaussianImageSmooth(true);
                im2surf->SetSmooth(true);
                im2surf->SetMedianFilter3D(true);
                im2surf->InterpolationOn();
                im2surf->SetGaussianStandardDeviation(bl);
                im2surf->SetMedianKernelSize(smth, smth, smth);
                im2surf->SetDecimate(mitk::ImageToSurfaceFilter::QuadricDecimation);
                im2surf->SetTargetReduction(ds);
                im2surf->UpdateLargestPossibleRegion();

                // mitk::ImageToSurfaceFilter::Pointer im2surf = mitk::ImageToSurfaceFilter::New();
                // im2surf->SetInput(veinsRelabeledImg);
                // im2surf->SetThreshold(0.5);
                // im2surf->SmoothOn();
                // im2surf->SetSmoothIteration(10);
                // im2surf->SetDecimate(mitk::ImageToSurfaceFilter::DecimatePro);
                // im2surf->SetTargetReduction(0.1);
                // im2surf->Update();

                mitk::Surface::Pointer shell = im2surf->GetOutput();
                // shell->SetVtkPolyData(im2surf->GetOutput());
                mitk::IOUtil::Save(shell, (path+"segmentation.vtk").toStdString());

                MITK_INFO(verbose) << "Add scalars into surface";
                scar->SetScarSegImage(veinsRelabeledImg);
                mitk::Surface::Pointer scarShell = scar->Scar3D(dir.toStdString(), veinsRelabeledImg);

                MITK_INFO(verbose) << ("Saving output shell to " + outputPath).toStdString();
                mitk::IOUtil::Save(scarShell, outputPath.toStdString());

            }

        }

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


void OutputImage(ImageType::Pointer im, QString dir, QString imName){
    QString impath = dir + imName;
    mitk::Image::Pointer outputImg = mitk::Image::New();
    mitk::CastToMitkImage(im, outputImg);
    mitk::IOUtil::Save(outputImg, impath.toStdString());
}

ThresholdType::Pointer thresholdImage(ImageType::Pointer input, uint16_t thresholdVal, QString debugOutput, bool debugging){
    ThresholdType::Pointer thresholdOutput = ThresholdType::New();
    thresholdOutput->SetInput(input);
    thresholdOutput->SetLowerThreshold(thresholdVal);
    thresholdOutput->SetUpperThreshold(thresholdVal);
    thresholdOutput->SetInsideValue(1);
    thresholdOutput->SetOutsideValue(0);
    thresholdOutput->Update();

    if(debugging){
        mitk::IOUtil::Save(mitk::ImportItkImage(thresholdOutput->GetOutput()), debugOutput.toStdString());
    }

    return thresholdOutput;
}
