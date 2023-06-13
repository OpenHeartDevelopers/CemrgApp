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
#include <QmitkIOUtil.h>
#include <mitkCommandLineParser.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkManualSegmentationToSurfaceFilter.h>

// VTK
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataNormals.h>

// ITK
#include <itkRelabelComponentImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkLabelObject.h>
#include <itkLabelMap.h>
#include <itkLabelImageToLabelMapFilter.h>
#include <itkLabelMapToLabelImageFilter.h>
#include <itkLabelSelectionLabelMapFilter.h>

// CemrgApp
#include <CemrgMeasure.h>

int main(int argc, char* argv[]) {

    mitkCommandLineParser parser;

    // Set general information about your command-line app
    parser.setCategory("Post processing");
    parser.setTitle("Morph Analysis App");
    parser.setContributor("CEMRG, KCL");
    parser.setDescription("Atrial Size Analysis");
    parser.setArgumentPrefix("--", "-");

    // Add arguments. Unless specified otherwise, each argument is optional.
    parser.addArgument(
        "directory", "d", mitkCommandLineParser::String,
        "Directory path", "Full path of directory",
        us::Any(), false);

    // Parse arguments.
    auto parsedArgs = parser.parseArguments(argc, argv);
    if (parsedArgs.empty())
        return EXIT_FAILURE;
    auto directory = us::any_cast<std::string>(parsedArgs["directory"]);

    try {

        QString path = QString::fromStdString(directory) + "/AnalyticBloodpool.nii";
        mitk::Image::Pointer analyticImage = mitk::IOUtil::Load<mitk::Image>(path.toStdString());

        if (analyticImage) {

            //Loop through labelled image
            typedef itk::Image<short, 3> ImageType;
            typedef itk::ImageRegionIteratorWithIndex<ImageType> ItType;
            ImageType::Pointer analyticItkImage = ImageType::New();
            CastToItkImage(analyticImage, analyticItkImage);
            ItType itLbl(analyticItkImage, analyticItkImage->GetRequestedRegion());
            for (itLbl.GoToBegin(); !itLbl.IsAtEnd(); ++itLbl) {
                if ((int)itLbl.Get() == 19 || (int)itLbl.Get() == 20) {
                    itLbl.Set(0);
                }//_if
            }//_for

            //Relabel the components to separate bloodpool and appendage
            typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> ConnectedComponentImageFilterType;
            ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
            connected->SetInput(analyticItkImage);
            connected->Update();
            typedef itk::RelabelComponentImageFilter<ImageType, ImageType> RelabelFilterType;
            RelabelFilterType::Pointer relabeler = RelabelFilterType::New();
            relabeler->SetInput(connected->GetOutput());
            relabeler->Update();

            //Keep the selected labels
            typedef itk::LabelObject<short, 3> LabelObjectType;
            typedef itk::LabelMap<LabelObjectType> LabelMapType;
            typedef itk::LabelImageToLabelMapFilter< ImageType, LabelMapType > LabelImageToLabelMapFilterType;
            LabelImageToLabelMapFilterType::Pointer labelMapConverter = LabelImageToLabelMapFilterType::New();
            labelMapConverter->SetInput(relabeler->GetOutput());
            labelMapConverter->SetBackgroundValue(0);
            typedef itk::LabelSelectionLabelMapFilter<LabelMapType> SelectorType;
            SelectorType::Pointer selector = SelectorType::New();
            selector->SetInput(labelMapConverter->GetOutput());
            selector->SetLabel(2);

            //Import to MITK image
            typedef itk::LabelMapToLabelImageFilter<LabelMapType, ImageType> LabelMapToLabelImageFilterType;
            LabelMapToLabelImageFilterType::Pointer labelImageConverter = LabelMapToLabelImageFilterType::New();
            labelImageConverter->SetInput(selector->GetOutput(0));
            labelImageConverter->Update();
            mitk::Image::Pointer ap = mitk::ImportItkImage(labelImageConverter->GetOutput());
            mitk::Image::Pointer bp = mitk::IOUtil::Load<mitk::Image>(directory + "/PVeinsCroppedImage.nii");

            mitk::IOUtil::Save(ap, directory + "/AP.nii.gz");
            mitk::IOUtil::Save(bp, directory + "/BP.nii.gz");

            //Ask for user input to set the parameters
            float th = 0.5;
            float bl = 0.8;
            int smth = 1;
            float ds = 0.5;

            auto filter1 = mitk::ManualSegmentationToSurfaceFilter::New();
            filter1->SetInput(bp);
            filter1->SetThreshold(th);
            filter1->SetUseGaussianImageSmooth(true);
            filter1->SetSmooth(true);
            filter1->SetMedianFilter3D(true);
            filter1->InterpolationOn();
            filter1->SetGaussianStandardDeviation(bl);
            filter1->SetMedianKernelSize(smth, smth, smth);
            filter1->SetDecimate(mitk::ImageToSurfaceFilter::QuadricDecimation);
            filter1->SetTargetReduction(ds);
            filter1->UpdateLargestPossibleRegion();
            mitk::Surface::Pointer shell1 = filter1->GetOutput();
            vtkSmartPointer<vtkPolyData> pd1 = shell1->GetVtkPolyData();
            pd1->SetVerts(nullptr);
            pd1->SetLines(nullptr);
            vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter1 = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
            connectivityFilter1->SetInputData(pd1);
            connectivityFilter1->ColorRegionsOff();
            connectivityFilter1->SetExtractionModeToLargestRegion();
            connectivityFilter1->Update();
            vtkSmartPointer<vtkPolyDataNormals> normals1 = vtkSmartPointer<vtkPolyDataNormals>::New();
            normals1->AutoOrientNormalsOn();
            normals1->FlipNormalsOff();
            normals1->SetInputConnection(connectivityFilter1->GetOutputPort());
            normals1->Update();
            shell1->SetVtkPolyData(normals1->GetOutput());
            mitk::Surface::Pointer surfLA = shell1->Clone();

            auto filter2 = mitk::ManualSegmentationToSurfaceFilter::New();
            filter2->SetInput(ap);
            filter2->SetThreshold(th);
            filter2->SetUseGaussianImageSmooth(true);
            filter2->SetSmooth(true);
            filter2->SetMedianFilter3D(true);
            filter2->InterpolationOn();
            filter2->SetGaussianStandardDeviation(bl);
            filter2->SetMedianKernelSize(smth, smth, smth);
            filter2->SetDecimate(mitk::ImageToSurfaceFilter::QuadricDecimation);
            filter2->SetTargetReduction(ds);
            filter2->UpdateLargestPossibleRegion();
            mitk::Surface::Pointer shell2 = filter2->GetOutput();
            vtkSmartPointer<vtkPolyData> pd2 = shell2->GetVtkPolyData();
            pd2->SetVerts(nullptr);
            pd2->SetLines(nullptr);
            vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter2 = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
            connectivityFilter2->SetInputData(pd2);
            connectivityFilter2->ColorRegionsOff();
            connectivityFilter2->SetExtractionModeToLargestRegion();
            connectivityFilter2->Update();
            vtkSmartPointer<vtkPolyDataNormals> normals2 = vtkSmartPointer<vtkPolyDataNormals>::New();
            normals2->AutoOrientNormalsOn();
            normals2->FlipNormalsOff();
            normals2->SetInputConnection(connectivityFilter2->GetOutputPort());
            normals2->Update();
            shell2->SetVtkPolyData(normals2->GetOutput());
            mitk::Surface::Pointer surfAP = shell2->Clone();

            //Volume and surface calculations
            std::unique_ptr<CemrgMeasure> morphAnal = std::unique_ptr<CemrgMeasure>(new CemrgMeasure());
            double surfceLA = morphAnal->calcSurfaceMesh(surfLA);
            double volumeLA = morphAnal->calcVolumeMesh(surfLA);
            double surfceAP = morphAnal->calcSurfaceMesh(surfAP);
            double volumeAP = morphAnal->calcVolumeMesh(surfAP);
            double sphereLA = morphAnal->GetSphericity(surfLA->GetVtkPolyData());

            //Store in text file
            ofstream morphResult;
            QString morphPath = QString::fromStdString(directory) + "/morphResults_AB.txt";
            morphResult.open(morphPath.toStdString(), std::ios_base::app);
            morphResult << "SA" << " " << surfceLA << "\n";
            morphResult << "VA" << " " << volumeLA << "\n";
            morphResult << "SP" << " " << surfceAP << "\n";
            morphResult << "VP" << " " << volumeAP << "\n";
            morphResult << "SF" << " " << sphereLA << "\n";
            morphResult.close();
        }//_if
    } catch (...) {
        return -1;
    }//_try
}
