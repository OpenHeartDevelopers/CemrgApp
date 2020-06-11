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
 *
 * Simple Common Utilities
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

//ITK
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>

//VTK
#include <vtkPolyData.h>

//Qmitk
#include <mitkBoundingObjectCutter.h>
#include <mitkProgressBar.h>
#include <mitkDataNode.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkIOUtil.h>
#include <mitkDataStorage.h>
#include <mitkRenderingManager.h>

//Qt
#include <QMessageBox>
#include <QString>
#include <QFile>
#include <QTextStream>
#include "CemrgCommonUtils.h"

mitk::DataNode::Pointer CemrgCommonUtils::imageNode;
mitk::DataNode::Pointer CemrgCommonUtils::cuttingNode;
mitk::Image::Pointer CemrgCommonUtils::imageToCut;
mitk::BoundingObject::Pointer CemrgCommonUtils::cuttingCube;

mitk::Image::Pointer CemrgCommonUtils::CropImage() {

    //Test input objects
    if (imageToCut.IsNull() || cuttingCube.IsNull())
        return NULL;

    //Prepare the cutter
    mitk::BoundingObjectCutter::Pointer cutter = mitk::BoundingObjectCutter::New();
    cutter->SetBoundingObject(cuttingCube);
    cutter->SetInput(imageToCut);
    cutter->AutoOutsideValueOff();

    //Actual cutting
    try {
        cutter->Update();
        mitk::ProgressBar::GetInstance()->Progress();
    } catch (const itk::ExceptionObject& e) {
        std::string message = std::string("The Cropping filter could not process because of: \n ") + e.GetDescription();
        QMessageBox::warning(
                    NULL, "Cropping not possible!", message.c_str(),
                    QMessageBox::Ok, QMessageBox::NoButton, QMessageBox::NoButton);
        return NULL;
    }//try

    //Cutting successful
    mitk::Image::Pointer resultImage = cutter->GetOutput();
    resultImage->DisconnectPipeline();
    resultImage->SetPropertyList(imageToCut->GetPropertyList()->Clone());
    mitk::ProgressBar::GetInstance()->Progress();

    return resultImage;
}

void CemrgCommonUtils::SetImageToCut(mitk::Image::Pointer imageToCut) {

    CemrgCommonUtils::imageToCut = imageToCut;
}

void CemrgCommonUtils::SetCuttingCube(mitk::BoundingObject::Pointer cuttingCube) {

    CemrgCommonUtils::cuttingCube = cuttingCube;
}

void CemrgCommonUtils::SetImageNode(mitk::DataNode::Pointer imageNode) {

    CemrgCommonUtils::imageNode = imageNode;
}

void CemrgCommonUtils::SetCuttingNode(mitk::DataNode::Pointer cuttingNode) {

    CemrgCommonUtils::cuttingNode = cuttingNode;
}

mitk::DataNode::Pointer CemrgCommonUtils::GetImageNode() {

    return imageNode;
}

mitk::DataNode::Pointer CemrgCommonUtils::GetCuttingNode() {

    return cuttingNode;
}

mitk::Image::Pointer CemrgCommonUtils::Downsample(mitk::Image::Pointer image, int factor) {

    typedef itk::Image<short,3> ImageType;
    typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleImageFilterType;
    typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> NearestInterpolatorType;

    //Cast to ITK
    ImageType::Pointer itkImage = ImageType::New();
    mitk::CastToItkImage(image, itkImage);

    ResampleImageFilterType::Pointer downsampler = ResampleImageFilterType::New();
    downsampler->SetInput(itkImage);

    NearestInterpolatorType::Pointer interpolator = NearestInterpolatorType::New();
    downsampler->SetInterpolator(interpolator);

    downsampler->SetDefaultPixelValue(0);

    ResampleImageFilterType::SpacingType spacing = itkImage->GetSpacing();
    spacing *= (double) factor;
    downsampler->SetOutputSpacing(spacing);

    downsampler->SetOutputOrigin(itkImage->GetOrigin());
    downsampler->SetOutputDirection(itkImage->GetDirection());

    ResampleImageFilterType::SizeType size = itkImage->GetLargestPossibleRegion().GetSize();
    for (int i=0; i<3; ++i)
        size[i] /= factor;

    downsampler->SetSize(size);
    downsampler->UpdateLargestPossibleRegion();

    //Save downsampled image
    image = mitk::ImportItkImage(downsampler->GetOutput())->Clone();
    mitk::ProgressBar::GetInstance()->Progress();
    return image;
}

mitk::Surface::Pointer CemrgCommonUtils::LoadVTKMesh(std::string path) {

    try {
        //Load the mesh
        mitk::Surface::Pointer surface = mitk::IOUtil::Load<mitk::Surface>(path);
        vtkSmartPointer<vtkPolyData> pd = surface->GetVtkPolyData();

        //Prepare points for MITK visualisation
        double Xmin = 0, Xmax = 0, Ymin = 0, Ymax = 0, Zmin = 0, Zmax = 0;
        for (int i=0; i<pd->GetNumberOfPoints(); i++) {
            double* point = pd->GetPoint(i);
            point[0] = -point[0];
            point[1] = -point[1];
            pd->GetPoints()->SetPoint(i, point);
            //Find mins and maxs
            if (i==0) {
                Xmin = point[0];
                Xmax = point[0];
                Ymin = point[1];
                Ymax = point[1];
                Zmin = point[2];
                Zmax = point[2];
            } else {
                if (point[0]<Xmin) Xmin = point[0];
                if (point[0]>Xmax) Xmax = point[0];
                if (point[1]<Ymin) Ymin = point[1];
                if (point[1]>Ymax) Ymax = point[1];
                if (point[2]<Zmin) Zmin = point[2];
                if (point[2]>Zmax) Zmax = point[2];
            }//_if
        }//_for
        double bounds[6] = {Xmin, Xmax, Ymin, Ymax, Zmin, Zmax};
        surface->GetGeometry()->SetBounds(bounds);

        return surface;

    } catch (...) {
        return mitk::Surface::New();
    }//_catch
}

QString CemrgCommonUtils::m3dlibParamFileGenerator(QString dir, QString filename, QString thicknessCalc){
    QString path2file = dir + mitk::IOUtil::GetDirectorySeparator() + filename;
    QFile fi(path2file);
    if(thicknessCalc.compare("0", Qt::CaseSensitive)!=0 && thicknessCalc.compare("1", Qt::CaseSensitive)!=0){
        MITK_INFO << "Thickness calculation set to default (OFF)";
        thicknessCalc = "0";
    }

    if (fi.open(QFile::WriteOnly | QFile::Truncate)) {
        QTextStream out(&fi);
        out << "[segmentation]" << "\n\n";
        out << "seg_dir" << "=" << "./example" << "\n";
        out << "seg_name" << "=" << "converted.inr" << "\n";
        out << "mesh_from_segmentation" << "=" << "1" << "\n\n";

        out << "[meshing]" << "\n\n";
        out << "readTheMesh" << "=" << "0" << "\n";
        out << "mesh_dir" << "=" << "." << "\n";
        out << "mesh_name" << "=" << "mesh" << "\n\n";

        out << "facet_angle" << "=" << "30"<< "\n";
        out << "facet_size" << "=" << "5.0" << "\n";
        out << "facet_distance" << "=" << "4"<< "\n";
        out << "cell_rad_edge_ratio" << "=" << "2.0" << "\n";
        out << "cell_size" << "=" << "1.0" << "\n\n";

        out << "rescaleFactor" << "=" << "1.0  # rescaling for carp and vtk output" << "\n\n";

        out << "[laplacesolver]" << "\n\n";
        out << "abs_toll" << "=" << "1e-6 # Also for evaluating the thickness" << "\n";
        out << "rel_toll" << "=" << "1e-6" << "\n";
        out << "itr_max" << "=" << "500" << "\n";
        out << "dimKrilovSp" << "=" << "150" << "\n";
        out << "verbose" << "=" << "0" << "\n\n";

        out << "[output]" << "\n\n";
        out << "outdir" << "=" << "." << "\n";
        out << "name" << "=" << "imgmesh" << "\n\n";

        out << "out_medit" << "=" << "0" << "\n";
        out << "out_carp" << "=" << "1" << "\n";
        out << "out_carp_binary" << "=" << "0" << "\n";
        out << "out_vtk" << "=" << "1" << "\n";
        out << "out_vtk_binary" << "=" << "0" << "\n";
        out << "out_potential" << "=" << "0" << "\n";
        out << "debug_output" << "=" << "0" << "\n";
        out << "debug_frequency" << "=" << "10000" << "\n\n";

        out << "[others]" << "\n\n";
        out << "eval_thickness" << "=" << thicknessCalc << "\n";
        out << "thickalgo" << "=" << "1" << "\n"; //#1: Martin Bishop Algorithm; 2: Cesare Corrado Algorithm
        out << "swapregions" << "=" << "1" << "\n";
        out << "verbose" << "=" << "0" << "\n";

        return path2file;
    } else {
        MITK_WARN << ("File " + path2file + "not created.").toStdString();
        return "ERROR_IN_PROCESSING";
    }

}

mitk::DataNode::Pointer CemrgCommonUtils::AddToStorage(
        mitk::BaseData* data, std::string nodeName, mitk::DataStorage::Pointer ds) {

    if (!data) return mitk::DataNode::New();

    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(data);
    node->SetName(nodeName);
    ds->Add(node);

    mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(ds);
    return node;
}
