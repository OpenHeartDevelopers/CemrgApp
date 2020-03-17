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
 * Simple Image Utilities for MITK
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

//ITK
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>

//Qmitk
#include <mitkBoundingObjectCutter.h>
#include <mitkProgressBar.h>
#include <mitkDataNode.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>

//Qt
#include <QMessageBox>

#include "CemrgImageUtils.h"


mitk::DataNode::Pointer CemrgImageUtils::imageNode;
mitk::DataNode::Pointer CemrgImageUtils::cuttingNode;
mitk::Image::Pointer CemrgImageUtils::imageToCut;
mitk::BoundingObject::Pointer CemrgImageUtils::cuttingCube;

mitk::Image::Pointer CemrgImageUtils::CropImage() {

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

void CemrgImageUtils::SetImageToCut(mitk::Image::Pointer imageToCut) {

    CemrgImageUtils::imageToCut = imageToCut;
}

void CemrgImageUtils::SetCuttingCube(mitk::BoundingObject::Pointer cuttingCube) {

    CemrgImageUtils::cuttingCube = cuttingCube;
}

void CemrgImageUtils::SetImageNode(mitk::DataNode::Pointer imageNode) {

    CemrgImageUtils::imageNode = imageNode;
}

void CemrgImageUtils::SetCuttingNode(mitk::DataNode::Pointer cuttingNode) {

    CemrgImageUtils::cuttingNode = cuttingNode;
}

mitk::DataNode::Pointer CemrgImageUtils::GetImageNode() {

    return imageNode;
}

mitk::DataNode::Pointer CemrgImageUtils::GetCuttingNode() {

    return cuttingNode;
}

mitk::Image::Pointer CemrgImageUtils::Downsample(mitk::Image::Pointer image, int factor) {

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
