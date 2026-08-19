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

// ITK
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkResampleImageFilter.h>
#include <itkOrientImageFilter.h>
#include <itkPasteImageFilter.h>

// VTK
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkCell.h>
#include <vtkImageMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor2D.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkCenterOfMass.h>
#include <vtkImageActor.h>
#include <vtkImageMapper3D.h>
#include <vtkExtractVOI.h>
#include <vtkPlane.h>
#include <vtkProperty.h>
#include <vtkCutter.h>
#include <vtkCamera.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkColorTransferFunction.h>
#include <vtkSphere.h>
#include <vtkSphereSource.h>
#include <vtkClipPolyData.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkCellDataToPointData.h>
#include <vtkPointDataToCellData.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkFillHolesFilter.h>

// Qmitk
#include <mitkBoundingObjectCutter.h>
#include <mitkProgressBar.h>
#include <mitkDataNode.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkIOUtil.h>
#include <mitkDataStorage.h>
#include <mitkImagePixelReadAccessor.h>
#include <mitkImageToSurfaceFilter.h>
#include <mitkManualSegmentationToSurfaceFilter.h>
#include <mitkRenderingManager.h>

// Qt
#include <QMessageBox>
#include <QString>
#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QTextStream>
#include <QJsonParseError>
#include <QJsonDocument>
#include <QJsonArray>

#include "CemrgCommonUtils.h"


mitk::DataNode::Pointer CemrgCommonUtils::imageNode;
mitk::DataNode::Pointer CemrgCommonUtils::cuttingNode;
mitk::Image::Pointer CemrgCommonUtils::imageToCut;
mitk::BoundingObject::Pointer CemrgCommonUtils::cuttingCube;

mitk::Image::Pointer CemrgCommonUtils::CropImage() {

    //Test input objects
    if (imageToCut.IsNull() || cuttingCube.IsNull()) {
        return nullptr;
    }

    //Prepare the cutter
    mitk::BoundingObjectCutter::Pointer cutter = mitk::BoundingObjectCutter::New();
    cutter->SetBoundingObject(cuttingCube);
    cutter->SetInput(imageToCut);
    cutter->AutoOutsideValueOff();

    //Actual cutting
    try {
        cutter->Update();
    } catch (const itk::ExceptionObject& e) {
        std::string message = std::string("The Cropping filter could not process because of: \n ") + e.GetDescription();
        QMessageBox::warning(
            NULL, "Cropping not possible!", message.c_str(),
            QMessageBox::Ok, QMessageBox::NoButton, QMessageBox::NoButton);
        return nullptr;
    }//try

    //Cutting successful
    mitk::Image::Pointer resultImage = cutter->GetOutput();
    resultImage->DisconnectPipeline();
    resultImage->SetPropertyList(imageToCut->GetPropertyList()->Clone());

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

    typedef itk::Image<short, 3> ImageType;
    typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleImageFilterType;
    typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> NearestInterpolatorType;

    //Cast to ITK
    ImageType::Pointer itkImage = ImageType::New();
    mitk::CastToItkImage(image, itkImage);

    //Downsampler
    ResampleImageFilterType::Pointer downsampler = ResampleImageFilterType::New();
    downsampler->SetInput(itkImage);
    NearestInterpolatorType::Pointer interpolator = NearestInterpolatorType::New();
    downsampler->SetInterpolator(interpolator);
    downsampler->SetDefaultPixelValue(0);
    ResampleImageFilterType::SpacingType spacing = itkImage->GetSpacing();
    spacing *= (double)factor;
    downsampler->SetOutputSpacing(spacing);
    downsampler->SetOutputOrigin(itkImage->GetOrigin());
    downsampler->SetOutputDirection(itkImage->GetDirection());
    ResampleImageFilterType::SizeType size = itkImage->GetLargestPossibleRegion().GetSize();
    for (int i = 0; i < 3; ++i)
        size[i] /= factor;
    downsampler->SetSize(size);
    downsampler->UpdateLargestPossibleRegion();

    //Save downsampled image
    image = mitk::ImportItkImage(downsampler->GetOutput())->Clone();
    return image;
}

mitk::Image::Pointer CemrgCommonUtils::IsoImageResampleReorient(mitk::Image::Pointer image, bool resample, bool reorientToRAI, bool isBinary) {

    MITK_INFO(resample) << "Resampling image to be isometric.";
    MITK_INFO(reorientToRAI) << "Doing a reorientation to RAI.";

    typedef itk::Image<short, 3> ImageType;
    typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleImageFilterType;
    ImageType::Pointer itkInputImage = ImageType::New();
    ImageType::Pointer resampleOutput, outputImage;
    mitk::CastToItkImage(image, itkInputImage);

    if (resample) {
        ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();

        if(isBinary){
            typedef itk::BSplineInterpolateImageFunction<ImageType, double, double> BSplineInterpolatorType;
            BSplineInterpolatorType::Pointer bsplineInterp = BSplineInterpolatorType::New();
            bsplineInterp->SetSplineOrder(3);
            resampler->SetInterpolator(bsplineInterp);
        } else{
            typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> NearestInterpolatorType;
            NearestInterpolatorType::Pointer nnInterp = NearestInterpolatorType::New();
            resampler->SetInterpolator(nnInterp);
        }

        resampler->SetInput(itkInputImage);
        resampler->SetOutputOrigin(itkInputImage->GetOrigin());
        ImageType::SizeType input_size = itkInputImage->GetLargestPossibleRegion().GetSize();
        ImageType::SpacingType input_spacing = itkInputImage->GetSpacing();
        ImageType::SizeType output_size;
        ImageType::SpacingType output_spacing;
        output_size[0] = input_size[0] * (input_spacing[0] / 1.0);
        output_size[1] = input_size[1] * (input_spacing[1] / 1.0);
        output_size[2] = input_size[2] * (input_spacing[2] / 1.0);
        output_spacing[0] = 1.0;
        output_spacing[1] = 1.0;
        output_spacing[2] = 1.0;
        resampler->SetSize(output_size);
        resampler->SetOutputSpacing(output_spacing);
        resampler->SetOutputDirection(itkInputImage->GetDirection());
        resampler->UpdateLargestPossibleRegion();
        resampleOutput = resampler->GetOutput();

    } else {
        resampleOutput = itkInputImage;
    }//_if

    if (reorientToRAI) {

        typedef itk::OrientImageFilter<ImageType, ImageType> OrientImageFilterType;
        OrientImageFilterType::Pointer orienter = OrientImageFilterType::New();
        orienter->UseImageDirectionOn();
        orienter->SetDesiredCoordinateOrientationToAxial(); // RAI
        orienter->SetInput(resampleOutput);
        orienter->Update();
        outputImage = orienter->GetOutput();

    } else {
        outputImage = resampleOutput;
    }//_if

    image = mitk::ImportItkImage(outputImage)->Clone();
    return image;
}

mitk::Image::Pointer CemrgCommonUtils::IsoImageResampleReorient(QString imPath, bool resample,  bool reorientToRAI, bool isBinary) {

    return CemrgCommonUtils::IsoImageResampleReorient(mitk::IOUtil::Load<mitk::Image>(imPath.toStdString()), resample, reorientToRAI, isBinary);
}

void CemrgCommonUtils::Binarise(mitk::Image::Pointer image, float background){
    using ImageType = itk::Image<float, 3>;
    using IteratorType = itk::ImageRegionIteratorWithIndex<ImageType>;

    ImageType::Pointer im = ImageType::New();
    mitk::CastToItkImage(image, im);

    IteratorType imIter(im, im->GetLargestPossibleRegion());

    imIter.GoToBegin();
    while(!imIter.IsAtEnd()){
        float value = (imIter.Get() > background) ? 1 : 0;
        imIter.Set(value);

        ++imIter;
    }

    image = mitk::ImportItkImage(im)->Clone();
}

mitk::Image::Pointer CemrgCommonUtils::ReturnBinarised(mitk::Image::Pointer image, float background){
    using ImageType = itk::Image<float, 3>;
    using IteratorType = itk::ImageRegionIteratorWithIndex<ImageType>;

    ImageType::Pointer im = ImageType::New();
    mitk::CastToItkImage(image, im);

    IteratorType imIter(im, im->GetLargestPossibleRegion());

    imIter.GoToBegin();
    while(!imIter.IsAtEnd()){
        float value = (imIter.Get() > background) ? 1 : 0;
        imIter.Set(value);

        ++imIter;
    }

    mitk::Image::Pointer outImg = mitk::Image::New();
    mitk::CastToMitkImage(im, outImg);
    return outImg;
}

bool CemrgCommonUtils::ConvertToNifti(mitk::BaseData::Pointer oneNode, QString path2file, bool resample, bool reorient) {

    bool successful = false;

    if (oneNode) {

        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(oneNode.GetPointer());
        if (image) { //Test if this data item is an image
            image = CemrgCommonUtils::IsoImageResampleReorient(image, resample, reorient);
            mitk::IOUtil::Save(image, path2file.toStdString());
            successful = true;
        } else {
            MITK_INFO << "[...] Problem casting node data to image";
        }//_if

    } else {
        MITK_INFO << "[...] Problem with node";
    }//_if

    return successful;
}

mitk::Image::Pointer CemrgCommonUtils::PadImageWithConstant(mitk::Image::Pointer image, int vxlsToExtend, short constant){
    using ImageType = itk::Image<short,3>;
    MITK_WARN(constant != 0) << "Constant != 0 not supported yet";

    ImageType::Pointer inputImg = ImageType::New();
    mitk::CastToItkImage(image, inputImg);

    ImageType::Pointer outputImg = ImageType::New();
    ImageType::IndexType start;
    start[0] = 0; start[1] = 0; start[2] = 0;

    ImageType::SizeType size;
    size[0] = image->GetDimension(0) + 2*vxlsToExtend;
    size[1] = image->GetDimension(1) + 2*vxlsToExtend;
    size[2] = image->GetDimension(2) + 2*vxlsToExtend;

    ImageType::RegionType region;
    region.SetSize(size);
    region.SetIndex(start);

    outputImg->SetRegions(region);
    outputImg->Allocate();

    using IteratorType = itk::ImageRegionIterator<ImageType>;
    IteratorType imIter(outputImg, outputImg->GetLargestPossibleRegion());
    imIter.GoToBegin();
    while(!imIter.IsAtEnd()){
        imIter.Set(0);
        ++imIter;
    }

    ImageType::IndexType indexForOutput;
    indexForOutput[0] = vxlsToExtend;
    indexForOutput[1] = vxlsToExtend;
    indexForOutput[2] = vxlsToExtend;

    using PasteImage = itk::PasteImageFilter<ImageType, ImageType>;
    PasteImage::Pointer paste = PasteImage::New();
    paste->SetSourceImage(inputImg);
    paste->SetSourceRegion(inputImg->GetLargestPossibleRegion());
    paste->SetDestinationImage(outputImg);
    paste->SetDestinationIndex(indexForOutput);

    image = mitk::ImportItkImage(paste->GetOutput())->Clone();

    return image;
}

void CemrgCommonUtils::SavePadImageWithConstant(QString inputPath, QString outputPath, int vxlsToExtend, short constant){
    QString out = (outputPath.isEmpty()) ? inputPath : outputPath;
    mitk::Image::Pointer inImg = mitk::IOUtil::Load<mitk::Image>(inputPath.toStdString());
    mitk::Image::Pointer outImg = CemrgCommonUtils::PadImageWithConstant(inImg, vxlsToExtend, constant);
    mitk::IOUtil::Save(outImg, out.toStdString());
}

bool CemrgCommonUtils::ImageConvertFormat(QString pathToImage, QString pathToOutput, bool optResample, bool optReorient, bool optImgBinary){
    mitk::Image::Pointer image = CemrgCommonUtils::IsoImageResampleReorient(pathToImage, optResample, optReorient, optImgBinary);
    if(optImgBinary){
        image = CemrgCommonUtils::ReturnBinarised(image);
    }

    mitk::IOUtil::Save(image, pathToOutput.toStdString());

    return (QFile::exists(pathToOutput));

}

void CemrgCommonUtils::SetSegmentationEdgesToZero(mitk::Image::Pointer image, QString outPath){
    using ImageType = itk::Image<short,3>;

    ImageType::Pointer im = ImageType::New();
    mitk::CastToItkImage(image, im);

    ImageType::SizeType size;
    size[0] = image->GetDimension(0);
    size[1] = image->GetDimension(1);
    size[2] = image->GetDimension(2);

    ImageType::IndexType pixelIndexStart, pixelIndexEnd;

    for (unsigned int ix = 0; ix < 3; ix++) {
        for (unsigned int jx = 0; jx < size[(ix + 1) % 3]; jx++) {
            for (unsigned int kx = 0; kx < size[(ix + 2) % 3]; kx++) {
                pixelIndexStart[ix] = 0;
                pixelIndexStart[(ix + 1) % 3] = jx;
                pixelIndexStart[(ix + 2) % 3] = kx;

                pixelIndexEnd[ix] = size[ix] - 1;
                pixelIndexEnd[(ix + 1) % 3] = jx;
                pixelIndexEnd[(ix + 2) % 3] = kx;

                im->SetPixel(pixelIndexStart, 0);
                im->SetPixel(pixelIndexEnd, 0);
            }
        }
    }

    mitk::Image::Pointer outImg = mitk::ImportItkImage(im)->Clone();
    if (!outPath.isEmpty()) {
        mitk::IOUtil::Save(outImg, outPath.toStdString());
    }
}

void CemrgCommonUtils::RoundPixelValues(QString pathToImage, QString outputPath) {
    QFileInfo fi(pathToImage);
    if (fi.exists()) {
        using ImageType = itk::Image<double, 3>;
        using IteratorType = itk::ImageRegionIteratorWithIndex<ImageType>;

        ImageType::Pointer im = ImageType::New();
        mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>(pathToImage.toStdString()), im);

        IteratorType imIter(im, im->GetLargestPossibleRegion());

        imIter.GoToBegin();
        while (!imIter.IsAtEnd()) {
            double pixelValue = imIter.Get();
            imIter.Set(std::round(pixelValue));

            ++imIter;
        }

        QString writingPath = (outputPath.isEmpty()) ? pathToImage : outputPath;

        mitk::Image::Pointer outputImg = mitk::Image::New();
        mitk::CastToMitkImage(im, outputImg);

        MITK_INFO(outputPath.isEmpty()) << ("Overwriting: " + pathToImage).toStdString();
        mitk::IOUtil::Save(outputImg, writingPath.toStdString());

    } else {
        MITK_WARN << ("Path: " + pathToImage + " does not exist.").toStdString();
    }
}


mitk::Surface::Pointer CemrgCommonUtils::LoadVTKMesh(std::string path) {

    try {
        //Load the mesh
        mitk::Surface::Pointer surface = mitk::IOUtil::Load<mitk::Surface>(path);
        vtkSmartPointer<vtkPolyData> pd = surface->GetVtkPolyData();

        //Prepare points for MITK visualisation
        double Xmin = 0, Xmax = 0, Ymin = 0, Ymax = 0, Zmin = 0, Zmax = 0;
        for (int i = 0; i < pd->GetNumberOfPoints(); i++) {
            double* point = pd->GetPoint(i);
            point[0] = -point[0];
            point[1] = -point[1];
            pd->GetPoints()->SetPoint(i, point);
            //Find mins and maxs
            if (i == 0) {
                Xmin = point[0];
                Xmax = point[0];
                Ymin = point[1];
                Ymax = point[1];
                Zmin = point[2];
                Zmax = point[2];
            } else {
                if (point[0] < Xmin) Xmin = point[0];
                if (point[0] > Xmax) Xmax = point[0];
                if (point[1] < Ymin) Ymin = point[1];
                if (point[1] > Ymax) Ymax = point[1];
                if (point[2] < Zmin) Zmin = point[2];
                if (point[2] > Zmax) Zmax = point[2];
            }//_if
        }//_for
        double bounds[6] = {Xmin, Xmax, Ymin, Ymax, Zmin, Zmax};
        surface->GetGeometry()->SetBounds(bounds);

        return surface;

    } catch (...) {
        return mitk::Surface::New();
    }//_catch
}

mitk::Surface::Pointer CemrgCommonUtils::ExtractSurfaceFromSegmentation(mitk::Image::Pointer image, double thresh, double blur, double smooth, double decimation) {
    auto im2surf = mitk::ManualSegmentationToSurfaceFilter::New();

    im2surf->SetInput(image);
    im2surf->SetThreshold(thresh);
    im2surf->SetUseGaussianImageSmooth(true);
    im2surf->SetSmooth(true);
    im2surf->SetMedianFilter3D(true);
    im2surf->InterpolationOn();
    im2surf->SetGaussianStandardDeviation(blur);
    im2surf->SetMedianKernelSize(smooth, smooth, smooth);
    im2surf->SetDecimate(mitk::ImageToSurfaceFilter::QuadricDecimation);
    im2surf->SetTargetReduction(decimation);
    im2surf->UpdateLargestPossibleRegion();

    mitk::Surface::Pointer shell = im2surf->GetOutput();
    return shell;
}

void CemrgCommonUtils::SetCellDataToPointData(mitk::Surface::Pointer surface, QString outputPath, QString fieldname){
    vtkSmartPointer<vtkCellDataToPointData> cell_to_point = vtkSmartPointer<vtkCellDataToPointData>::New();
    cell_to_point->SetInputData(surface->GetVtkPolyData());
    cell_to_point->PassCellDataOn();
    cell_to_point->SetContributingCellOption(0); // All=0, Patch=1, DataSetMax=2
    cell_to_point->Update();
    surface->SetVtkPolyData(cell_to_point->GetPolyDataOutput());
    surface->GetVtkPolyData()->GetPointData()->GetScalars()->SetName(fieldname.toStdString().c_str());

    if(!outputPath.isEmpty()){
        mitk::IOUtil::Save(surface, outputPath.toStdString());
    }
}

void CemrgCommonUtils::SetPointDataToCellData(mitk::Surface::Pointer surface, bool categories, QString outputPath){
    vtkSmartPointer<vtkPointDataToCellData> point_to_cell = vtkSmartPointer<vtkPointDataToCellData>::New();
    point_to_cell->SetInputData(surface->GetVtkPolyData());
    point_to_cell->PassPointDataOn();
    point_to_cell->SetCategoricalData(categories);
    point_to_cell->Update();
    surface->SetVtkPolyData(point_to_cell->GetPolyDataOutput());

    if(!outputPath.isEmpty()){
        mitk::IOUtil::Save(surface, outputPath.toStdString());
    }
}

mitk::Surface::Pointer CemrgCommonUtils::ClipWithSphere(mitk::Surface::Pointer surface, double x_c, double y_c, double z_c, double radius, QString saveToPath){
    double centre[3] = {x_c, y_c, z_c};
    //Clipper
    vtkSmartPointer<vtkSphere> sphere = vtkSmartPointer<vtkSphere>::New();
    sphere->SetCenter(centre);
    sphere->SetRadius(radius);
    vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
    clipper->SetClipFunction(sphere);
    clipper->SetInputData(surface->GetVtkPolyData());
    clipper->InsideOutOff();
    clipper->Update();

    if (!saveToPath.isEmpty()) {
        MITK_INFO << ("Saving clipper sphere to: " + saveToPath).toStdString();
        mitk::Surface::Pointer outSphere = mitk::Surface::New();
        vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetCenter(centre[0], centre[1], centre[2]);
        sphereSource->SetRadius(radius);
        sphereSource->SetPhiResolution(40);
        sphereSource->SetThetaResolution(40);
        sphereSource->Update();

        outSphere->SetVtkPolyData(sphereSource->GetOutput());
        mitk::IOUtil::Save(outSphere, saveToPath.toStdString());
    }

    //Extract and clean surface mesh
    vtkSmartPointer<vtkDataSetSurfaceFilter> surfer = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfer->SetInputData(clipper->GetOutput());
    surfer->Update();

    vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->SetInputConnection(surfer->GetOutputPort());
    cleaner->Update();

    vtkSmartPointer<vtkPolyDataConnectivityFilter> lrgRegion = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    lrgRegion->SetInputConnection(cleaner->GetOutputPort());
    lrgRegion->SetExtractionModeToLargestRegion();
    lrgRegion->Update();

    cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->SetInputConnection(lrgRegion->GetOutputPort());
    cleaner->Update();

    //Return the clipped mesh
    surface->SetVtkPolyData(cleaner->GetOutput());

    return surface;
}

void CemrgCommonUtils::FlipXYPlane(mitk::Surface::Pointer surf, QString dir, QString vtkname) {

    //Prepare points for MITK visualisation - (CemrgCommonUtils::LoadVTKMesh)
    vtkSmartPointer<vtkPolyData> pd = surf->GetVtkPolyData();
    for (int ix = 0; ix < pd->GetNumberOfPoints(); ix++) {
        double* point = pd->GetPoint(ix);
        point[0] = -point[0];
        point[1] = -point[1];
        pd->GetPoints()->SetPoint(ix, point);
    }

    if (!vtkname.isEmpty()) {
        vtkname += (!vtkname.contains(".vtk")) ? ".vtk" : "";
        QString path = dir + "/" + vtkname;
        mitk::IOUtil::Save(surf, path.toStdString());
    }
}

QString CemrgCommonUtils::M3dlibParamFileGenerator(QString dir, QString filename, QString thicknessCalc) {

    QString path2file = dir + "/" + filename;
    QFile fi(path2file);

    if (thicknessCalc.compare("0", Qt::CaseSensitive) != 0 && thicknessCalc.compare("1", Qt::CaseSensitive) != 0) {
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

        out << "facet_angle" << "=" << "30" << "\n";
        out << "facet_size" << "=" << "5.0" << "\n";
        out << "facet_distance" << "=" << "4" << "\n";
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

QString CemrgCommonUtils::OpenCarpParamFileGenerator(QString dir, QString filename, QString meshname, QString zeroBoundaryName, QString oneBoundaryName){
    QString path2file = dir + "/" + filename;
    QDir home(dir);
    QFile fi(path2file);

    if (fi.open(QFile::WriteOnly | QFile::Truncate)) {
        QTextStream out(&fi);

        out << "meshname " << "= " << meshname << "\n";
        out << "experiment " << "= " << " 2" << "\n";
        out << "bidomain " << "= " << " 1" << "\n";

        out << "num_imp_regions " << "= " << " 1" << "\n";
        out << "imp_region[0].im " << "= " << " MBRDR" << "\n";
        out << "num_gregions " << "= " << " 1" << "\n";

        out << "gregion[0].g_et " << "= " << " 1" << "\n";
        out << "gregion[0].g_el " << "= " << " 1" << "\n";
        out << "gregion[0].g_en " << "= " << " 1" << "\n";
        out << "gregion[0].g_il " << "= " << " 1" << "\n";
        out << "gregion[0].g_it " << "= " << " 1" << "\n";
        out << "gregion[0].g_in " << "= " << " 1" << "\n";
        out << "gregion[0].num_IDs " << "= " << " 1" << "\n";
        out << "gregion[0].ID[0] " << "= " << " 2" << "\n";

        out << "ellip_use_pt " << "= " << " 1" << "\n";
        out << "parab_use_pt " << "= " << " 1" << "\n";

        int ix = 0;
        bool onlyOneStimulusCheck = zeroBoundaryName.isEmpty();
        if(onlyOneStimulusCheck){
            out << "num_stim " << "= " << " 1" << "\n";

        } else {
            out << "num_stim " << "= " << " 2" << "\n";
            out << "stimulus["+ QString::number(ix)+"].stimtype " << "= " << " 3" << "\n";
            out << "stimulus["+ QString::number(ix)+"].vtx_file " << "= " << zeroBoundaryName << "\n";
            ix++;
        }

        out << "stimulus["+ QString::number(ix)+"].stimtype " << "= " << " 2" << "\n";
        out << "stimulus["+ QString::number(ix)+"].duration " << "= " << " 1" << "\n";
        out << "stimulus["+ QString::number(ix)+"].strength " << "= " << " 1.0" << "\n";

        if(onlyOneStimulusCheck){
            out << "stimulus["+ QString::number(ix)+"].vtx_fcn  " << "= " << " 1" << "\n";
        }
        out << "stimulus["+ QString::number(ix)+"].vtx_file " << "= " << oneBoundaryName << "\n";

        return path2file;

    } else {
        MITK_WARN << ("File " + path2file + "not created.").toStdString();
        return "ERROR_IN_PROCESSING";
    }
}

bool CemrgCommonUtils::ConvertToCarto(
        std::string vtkPath, std::vector<double> thresholds, double meanBP, double stdvBP, int methodType, bool discreteScheme) {

    //Read vtk from the file
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(vtkPath.c_str());
    reader->Update();
    vtkSmartPointer<vtkPolyData> pd = reader->GetOutput();

    //Output path
    QString qoutputPath = QString::fromStdString(vtkPath);
    std::string outputPath = qoutputPath.left(qoutputPath.lastIndexOf(QChar('.'))).toStdString();
    outputPath = outputPath + "-carto.vtk";

    //File
    std::ofstream cartoFile;
    cartoFile.open(outputPath);

    //Header
    cartoFile << "# vtk DataFile Version 3.0\n";
    cartoFile << "PatientData Anon Anon 00000000\n";
    cartoFile << "ASCII\n";
    cartoFile << "DATASET POLYDATA\n";

    //Points
    cartoFile << "POINTS\t" << pd->GetNumberOfPoints() << "\tfloat\n";
    for (int i = 0; i < pd->GetNumberOfPoints(); i++) {
        double* point = pd->GetPoint(i);
        cartoFile << point[0] << " " << point[1] << " " << point[2] << "\n";
    }
    cartoFile << "\n";

    //Cells
    cartoFile << "POLYGONS\t";
    cartoFile << pd->GetNumberOfCells() << "\t";
    cartoFile << pd->GetNumberOfCells() * 4 << "\n";
    for (int i = 0; i < pd->GetNumberOfCells(); i++) {
        vtkCell* cell = pd->GetCell(i);
        vtkIdList* list = cell->GetPointIds();
        cartoFile << "3";
        for (int j = 0; j < list->GetNumberOfIds(); j++)
            cartoFile << " " << list->GetId(j);
        cartoFile << "\n";
    }

    //Point data
    vtkSmartPointer<vtkFloatArray> pointData = vtkSmartPointer<vtkFloatArray>::New();
    try {
        if (pd->GetCellData()->GetScalars() != NULL) {

            vtkSmartPointer<vtkCellDataToPointData> cellToPoint = vtkSmartPointer<vtkCellDataToPointData>::New();
            cellToPoint->SetInputData(pd);
            cellToPoint->PassCellDataOn();
            cellToPoint->Update();
            pointData = vtkFloatArray::SafeDownCast(cellToPoint->GetPolyDataOutput()->GetPointData()->GetScalars());

        } else
            throw;
    } catch (...) {
        MITK_ERROR << "Storing point data failed! Check your input";
        return false;
    }//_try

    float min = pointData->GetRange()[0];
    float max = pointData->GetRange()[1];

    MITK_INFO << "Storing point data, number of tuples: " << pointData->GetNumberOfTuples();
    MITK_INFO << "Storing point data, number of components: " << pointData->GetNumberOfComponents();

    if (pointData->GetNumberOfTuples() != 0) {

        cartoFile << "\nPOINT_DATA\t";
        cartoFile << pointData->GetNumberOfTuples() << "\n";

        if (pointData->GetNumberOfComponents() == 1) {

            cartoFile << "SCALARS scalars float\n";
            cartoFile << "LOOKUP_TABLE lookup_table\n";
            for (int i = 0; i < pointData->GetNumberOfTuples(); i++) {

                //Get scalar raw value
                double value = static_cast<double>(pointData->GetTuple1(i));

                //Colouring
                if (discreteScheme) {
                    if (methodType == 1) {
                        if (value < (meanBP * thresholds.at(0))) value = 0.0;
                        else if (thresholds.size() == 2 && value < (meanBP * thresholds.at(1))) value = 0.5;
                        else value = 1.0;
                    } else {
                        if (value < (meanBP + thresholds.at(0) * stdvBP)) value = 0.0;
                        else if (thresholds.size() == 2 && value < (meanBP + thresholds.at(1) * stdvBP)) value = 0.5;
                        else value = 1.0;
                    }//_if
                } else {
                    value = (value - min) / (max - min);
                }//_if

                std::stringstream stream;
                stream << std::fixed << std::setprecision(2) << value;
                cartoFile << stream.str() << "\n";

            }//_for
            cartoFile << "\n";

        } else {

            for (int i = 0; pointData->GetNumberOfComponents(); i++) {

                cartoFile << "SCALARS " << "scalars" << i << " float\n";
                cartoFile << "LOOKUP_TABLE lookup_table\n";
                for (int j = 0; j < pointData->GetNumberOfTuples(); j++)
                    cartoFile << pointData->GetTuple(j)[i] << " ";
                cartoFile << "\n";

            }//_for
        }//_if
    }//_point_data

    MITK_INFO << "Storing lookup table, min/max scalar values: " << min << " " << max;

    //LUT
    int numCols = discreteScheme ? 3 : 256;
    cartoFile << "LOOKUP_TABLE lookup_table " << numCols << "\n";
    vtkSmartPointer<vtkColorTransferFunction> lut = vtkSmartPointer<vtkColorTransferFunction>::New();
    lut->SetColorSpaceToRGB();
    lut->AddRGBPoint(0.0, 0.04, 0.21, 0.25);
    lut->AddRGBPoint((numCols - 1.0) / 2.0, 0.94, 0.47, 0.12);
    lut->AddRGBPoint((numCols - 1.0), 0.90, 0.11, 0.14);
    lut->SetScaleToLinear();
    for (int i = 0; i < numCols; i++) {
        cartoFile << lut->GetColor(i)[0] << " ";
        cartoFile << lut->GetColor(i)[1] << " ";
        cartoFile << lut->GetColor(i)[2] << " ";
        cartoFile << "1.0" << "\n";
    }//_for

    cartoFile.close();
    return true;
}

void CemrgCommonUtils::MotionTrackingReport(QString directory, int timePoints) {

    for (int tS = 0; tS < timePoints; tS++) {

        //Image
        QString path = directory + "/dcm-" + QString::number(tS) + ".nii";
        mitk::Image::Pointer img3D = mitk::IOUtil::Load<mitk::Image>(path.toStdString());
        int* extent = img3D->GetVtkImageData()->GetExtent();
        mitk::Vector3D spacing = img3D->GetGeometry()->GetSpacing();
        vtkSmartPointer<vtkMatrix4x4> direction = img3D->GetGeometry()->GetVtkMatrix();

        //Mesh
        path = directory + "/Model-" + QString::number(tS) + ".vtk";
        mitk::Surface::Pointer sur3D = CemrgCommonUtils::LoadVTKMesh(path.toStdString());

        //Window
        int xSliceMin = extent[0];
        int xSliceMax = extent[1];
        int ySliceMin = extent[2];
        int ySliceMax = extent[3];
        //int zSliceMin = extent[4];
        int zSliceMax = extent[5];
        double xmins[8] = {0.00, 0.25, 0.50, 0.75, 0.00, 0.25, 0.50, 0.75};
        double xmaxs[8] = {0.25, 0.50, 0.75, 1.00, 0.25, 0.50, 0.75, 1.00};
        double ymins[8] = {0.00, 0.00, 0.00, 0.00, 0.50, 0.50, 0.50, 0.50};
        double ymaxs[8] = {0.50, 0.50, 0.50, 0.50, 1.00, 1.00, 1.00, 1.00};
        vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->SetAlphaBitPlanes(1);
        renderWindow->SetSize(500, 500);

        for (int view = 0; view < 8; view++) {

            //Setup views
            int zSlice = zSliceMax - view * floor(zSliceMax / 8);
            double zPlane = zSlice * spacing[2];

            //Image mapper
            vtkSmartPointer<vtkExtractVOI> extractSlice = vtkSmartPointer<vtkExtractVOI>::New();
            extractSlice->SetInputData(img3D->GetVtkImageData());
            extractSlice->SetVOI(xSliceMin, xSliceMax, ySliceMin, ySliceMax, zSlice, zSlice);
            extractSlice->Update();
            vtkSmartPointer<vtkImageData> slice = extractSlice->GetOutput();
            vtkSmartPointer<vtkImageActor> imgActor = vtkSmartPointer<vtkImageActor>::New();
            imgActor->GetMapper()->SetInputData(slice);

            //Mesh mapper
            vtkSmartPointer<vtkPolyData> pd = sur3D->GetVtkPolyData();
            vtkSmartPointer<vtkTransform> scaling = vtkSmartPointer<vtkTransform>::New();
            scaling->Scale(spacing[0], spacing[1], spacing[2]);
            vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
            transform->SetMatrix(direction);
            transform->Inverse();
            transform->PostMultiply();
            transform->Concatenate(scaling);
            vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
            transformFilter->SetInputData(pd);
            transformFilter->SetTransform(transform);
            transformFilter->Update();
            pd = transformFilter->GetOutput();
            vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
            plane->SetOrigin(0, 0, zPlane);
            plane->SetNormal(0, 0, 1);
            vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
            cutter->SetCutFunction(plane);
            cutter->SetInputData(pd);
            cutter->Update();
            vtkSmartPointer<vtkPolyDataMapper> mapMesh = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapMesh->SetInputConnection(cutter->GetOutputPort());
            mapMesh->SetScalarModeToUsePointData();
            mapMesh->SetScalarVisibility(1);
            mapMesh->SetScalarRange(1, 4);
            vtkSmartPointer<vtkActor> mshActor = vtkSmartPointer<vtkActor>::New();
            mshActor->GetProperty()->SetRepresentationToPoints();
            mshActor->GetProperty()->SetPointSize(3);
            mshActor->SetMapper(mapMesh);

            //Image renderer
            vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
            renderer->AddActor(imgActor);
            renderer->AddActor(mshActor);
            renderer->ResetCamera();
            renderer->GetActiveCamera()->ParallelProjectionOn();
            renderer->GetActiveCamera()->SetParallelScale(.5 * imgActor->GetBounds()[1]);
            renderWindow->AddRenderer(renderer);
            renderer->SetViewport(xmins[view], ymins[view], xmaxs[view], ymaxs[view]);
            renderWindow->Render();

        }//_for

        //Screenshot
        vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
        windowToImageFilter->SetInput(renderWindow);
        windowToImageFilter->SetInputBufferTypeToRGBA();
        windowToImageFilter->ReadFrontBufferOff();
        windowToImageFilter->FixBoundaryOn();
        windowToImageFilter->Update();
        vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
        writer->SetFileName((directory.toStdString() + "/dcm-" + QString::number(tS).toStdString() + ".png").c_str());
        writer->SetInputConnection(windowToImageFilter->GetOutputPort());
        writer->Write();
    }
}

void CemrgCommonUtils::CalculatePolyDataNormals(vtkSmartPointer<vtkPolyData>& pd, bool celldata) {

    vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
    vtkSmartPointer<vtkPolyData> tempPD = vtkSmartPointer<vtkPolyData>::New();
    tempPD->DeepCopy(pd);
    if (celldata) {
        normals->ComputeCellNormalsOn();
    } else { // pointdata
        normals->ComputePointNormalsOn();
    }
    normals->SetInputData(tempPD);
    normals->SplittingOff();
    normals->Update();
    pd = normals->GetOutput();
}

mitk::DataNode::Pointer CemrgCommonUtils::AddToStorage(
    mitk::BaseData* data, std::string nodeName, mitk::DataStorage::Pointer ds, bool init) {

    if (!data)
        return mitk::DataNode::New();

    //DS node creation
    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData(data);
    node->SetName(nodeName);
    ds->Add(node);

    if (init)
        mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(ds);

    return node;
}

QJsonObject CemrgCommonUtils::ReadJSONFile(QString dir, QString fname) {
    fname += (!fname.endsWith(".json")) ? ".json" : "";
    QFile file(dir + "/" + fname);

    if (!file.open(QIODevice::ReadOnly)) {
        qWarning("Failed to open file");
        return QJsonObject();
    }

    QByteArray jsonData = file.readAll();

    QJsonParseError error;
    QJsonDocument jsonDoc = QJsonDocument::fromJson(jsonData, &error);

    MITK_WARN((jsonDoc.isNull())) << ("Failed to parse JSON: " + error.errorString()).toStdString();

    QJsonObject jsonObj = jsonDoc.object();

    file.close();

    return jsonObj;
}

bool CemrgCommonUtils::WriteJSONFile(QJsonObject json, QString dir, QString fname) {
    // Create the JSON document
    QJsonDocument jsonDoc(json);
    bool success = true;

    // Open the file for writing
    fname += (!fname.endsWith(".json")) ? ".json" : "";
    QFile file(dir + "/" + fname);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
        qWarning("Failed to open file");
        success = false;
    }

    file.write(jsonDoc.toJson());
    file.close();

    return success;
}

bool CemrgCommonUtils::ModifyJSONFile(QString dir, QString fname, QString key, QString value, QString type) {
    bool success = true;
    QJsonObject json = ReadJSONFile(dir, fname);

    if (json.empty()) {
        return false;
    }

    if (value.isEmpty()) {
        json.remove(key);
    } else {
        // Determine the data type of the value
        QJsonValue jsonValue;
        QJsonArray json_array;
        if (type == "string") {
            jsonValue = QJsonValue::fromVariant(value);
        }
        else if (type == "int") {
            jsonValue = QJsonValue::fromVariant(value.toInt());
        }
        else if (type == "float" || type == "double") {
            jsonValue = QJsonValue::fromVariant(value.toDouble());
        }
        else if (type == "bool") {
            jsonValue = QJsonValue::fromVariant(value.toLower() == "true");
        }
        else if (type == "array") {
            // split value
            QStringList array_values = value.split(",");
            for (int ix = 0; ix < array_values.size(); ix++) {
                // QJsonValue json_value = array_values.at(ix).toDouble();
                json_array.push_back(array_values.at(ix).toDouble());
            }
        }
        else {
            // If the data type is not recognized, use a null value
            jsonValue = QJsonValue();
        }

        // Add the key-value pair to the object
        json[key] = (type == "array") ? json_array : jsonValue;
    }

    success = WriteJSONFile(json, dir, fname);
    return success;
}

QJsonObject CemrgCommonUtils::CreateJSONObject(QStringList keys_list, QStringList values_list, QStringList types_list) {
    QJsonObject jsonObj;

    for (int i = 0; i < keys_list.size(); i++) {
        QString key = keys_list.at(i);
        QString value = values_list.at(i);
        QString type = types_list.at(i);

        // Determine the data type of the value
        QJsonValue jsonValue;
        QJsonArray json_array;
        if (type == "string") {
            jsonValue = QJsonValue::fromVariant(value);
        }
        else if (type == "int") {
            jsonValue = QJsonValue::fromVariant(value.toInt());
        }
        else if (type == "float" || type == "double") {
            jsonValue = QJsonValue::fromVariant(value.toDouble());
        }
        else if (type == "bool") {
            jsonValue = QJsonValue::fromVariant(value.toLower() == "true");
        }
        else  if (type == "array") {
            // split value
            QStringList array_values = value.split(",");
            for (int ix = 0; ix < array_values.size(); ix++) {
                // QJsonValue json_value = array_values.at(ix).toDouble();
                json_array.push_back(array_values.at(ix).toDouble());
            }
        }
        else {
            // If the data type is not recognized, use a null value
            jsonValue = QJsonValue();
        }

        // Add the key-value pair to the object
        jsonObj[key] = (type == "array") ? json_array : jsonValue;
    }

    return jsonObj;
}

mitk::Image::Pointer CemrgCommonUtils::ImageFromSurfaceMesh(mitk::Surface::Pointer surf, double origin[3], double spacing[3]){
    vtkSmartPointer<vtkPolyData> pd = surf->GetVtkPolyData();
    double bounds[6];
    pd->GetBounds(bounds);

    int dimensions[3];
    for (int ix = 0; ix < 3; ix++) {
        dimensions[ix] = static_cast<int>(std::ceil((bounds[ix * 2 + 1] - bounds[ix * 2]) / spacing[ix]));
    }

    for (int jx = 0; jx < 3; jx++) {
        origin[jx] = bounds[2*jx] + spacing[jx]/2;
    }
    std::cout << "o = (" << origin[0] << ", " << origin[1] << ", " << origin[2] << ")"<< '\n';

    //Prepare empty image
    vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();
    whiteImage->SetOrigin(origin);
    whiteImage->SetSpacing(spacing);
    whiteImage->SetDimensions(dimensions);
    whiteImage->SetExtent(0, dimensions[0] - 1, 0, dimensions[1] - 1, 0, dimensions[2] - 1);
    whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

    unsigned char inval = 1;
    unsigned char otval = 0;

    vtkIdType count = whiteImage->GetNumberOfPoints();
    for (vtkIdType i = 0; i < count; ++i){
        whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
    }

    //Image Stencil
    vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
    pol2stenc->SetTolerance(0.5);
    pol2stenc->SetInputData(pd);
    pol2stenc->SetOutputOrigin(origin);
    pol2stenc->SetOutputSpacing(spacing);
    pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
    pol2stenc->Update();

    vtkSmartPointer<vtkImageStencil> imgstenc = vtkSmartPointer<vtkImageStencil>::New();
    imgstenc->SetInputData(whiteImage);
    imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
    imgstenc->ReverseStencilOff();
    imgstenc->SetBackgroundValue(otval);
    imgstenc->Update();

    //VTK to ITK conversion
    mitk::Image::Pointer cutImg = mitk::Image::New();
    cutImg->Initialize(imgstenc->GetOutput());
    cutImg->SetVolume(imgstenc->GetOutput()->GetScalarPointer());

    return cutImg;

}

void CemrgCommonUtils::SaveImageFromSurfaceMesh(QString surfPath, double origin[3], double spacing[3], QString outputPath){
    QString out;

    if(outputPath.isEmpty()){
        QFileInfo fi(surfPath);
        out = fi.absolutePath() + "/" + fi.baseName() + ".nii";
    } else{
        out = outputPath;
    }
    mitk::Surface::Pointer surf = mitk::IOUtil::Load<mitk::Surface>(surfPath.toStdString());
    mitk::Image::Pointer im = CemrgCommonUtils::ImageFromSurfaceMesh(surf, origin, spacing);

    mitk::IOUtil::Save(im, out.toStdString());
}

void CemrgCommonUtils::FillHoles(mitk::Surface::Pointer surf, QString dir, QString vtkname){
    vtkSmartPointer<vtkPolyData> pd = surf->GetVtkPolyData();
    vtkSmartPointer<vtkFillHolesFilter> fillholes = vtkSmartPointer<vtkFillHolesFilter>::New();
    fillholes->SetInputData(pd);
    fillholes->Update();

    surf->SetVtkPolyData(fillholes->GetOutput());

    if (!dir.isEmpty() && !vtkname.isEmpty()) {
        vtkname += (!vtkname.contains(".vtk")) ? ".vtk" : "";
        QString outPath = dir + "/" + vtkname;
        mitk::IOUtil::Save(surf, outPath.toStdString());
    }
}

double CemrgCommonUtils::GetSphereParametersFromLandmarks(mitk::PointSet::Pointer landmarks, double * centre){
    //Retrieve mean and distance of 3 points
    double x_c = 0;
    double y_c = 0;
    double z_c = 0;
    for(int i=0; i<landmarks->GetSize(); i++) {
        x_c = x_c + landmarks->GetPoint(i).GetElement(0);
        y_c = y_c + landmarks->GetPoint(i).GetElement(1);
        z_c = z_c + landmarks->GetPoint(i).GetElement(2);
    }//_for
    x_c /= landmarks->GetSize();
    y_c /= landmarks->GetSize();
    z_c /= landmarks->GetSize();
    double * distance = new double [landmarks->GetSize()];
    for(int i=0; i<landmarks->GetSize(); i++) {
        double x_d = landmarks->GetPoint(i).GetElement(0) - x_c;
        double y_d = landmarks->GetPoint(i).GetElement(1) - y_c;
        double z_d = landmarks->GetPoint(i).GetElement(2) - z_c;
        distance[i] = sqrt(pow(x_d,2) + pow(y_d,2) + pow(z_d,2));
    }//_for
    double radius = *std::max_element(distance, distance + landmarks->GetSize());
    centre[0] = x_c;
    centre[1] = y_c;
    centre[2] = z_c;

    return radius;
}

//UTILities for CARP - operations with .elem and .pts files
void CemrgCommonUtils::OriginalCoordinates(QString imagePath, QString pointPath, QString outputPath, double scaling) {
    if (QFileInfo::exists(imagePath) && QFileInfo::exists(pointPath)) {
        typedef itk::Image<uint8_t, 3> ImageType;
        mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(imagePath.toStdString());
        ImageType::Pointer itkInput = ImageType::New();
        ImageType::PointType origin;
        mitk::CastToItkImage(image, itkInput);
        origin = itkInput->GetOrigin();

        std::ifstream pointFileRead;

        int nPts;
        pointFileRead.open(pointPath.toStdString());
        pointFileRead >> nPts;

        double x, y, z;

        std::ofstream outputFileWrite;
        outputFileWrite.open(outputPath.toStdString());
        outputFileWrite << nPts << std::endl;

        for (int iPt = 0; iPt < nPts; iPt++) {
            double xt, yt, zt;
            pointFileRead >> x;
            pointFileRead >> y;
            pointFileRead >> z;
            if (pointFileRead.eof()) {
                MITK_WARN << " WARNING!: File ended prematurely " << endl;
                break;
            }

            xt = x + (origin[0] * scaling);
            yt = y + (origin[1] * scaling);
            zt = z + (origin[2] * scaling);

            outputFileWrite << std::fixed << xt << " ";
            outputFileWrite << std::fixed << yt << " ";
            outputFileWrite << std::fixed << zt << std::endl;
        }
        pointFileRead.close();
        outputFileWrite.close();
        MITK_INFO << ("Saved to file: " + outputPath).toStdString();

    } else {
        MITK_ERROR(QFileInfo::exists(imagePath)) << ("Could not read file" + imagePath).toStdString();
        MITK_ERROR(QFileInfo::exists(pointPath)) << ("Could not read file" + pointPath).toStdString();
    }

}

void CemrgCommonUtils::CalculateCentreOfGravity(QString pointPath, QString elemPath, QString outputPath) {
    if (QFileInfo::exists(elemPath) && QFileInfo::exists(pointPath)) {
        FILE* pointFileRead = fopen(pointPath.toStdString().c_str(), "r");
        FILE* elemFileRead = fopen(elemPath.toStdString().c_str(), "r");
        int nPts, nElem;

        int fscanOut1 = fscanf(pointFileRead, "%d\n", &nPts);
        MITK_INFO(fscanOut1 == 1) << "File size read correctly";

        double* pts_array = (double*)malloc(nPts * 3 * sizeof(double));
        if (pts_array == NULL) {
            MITK_ERROR << "pts_array malloc FAIL";
            fclose(pointFileRead);
            fclose(elemFileRead);
            return;
        }

        MITK_INFO << "Beginning input .pts file";
        for (int i = 0; i < nPts; i++) {
            double* loc = pts_array + 3 * i;
            fscanOut1 = fscanf(pointFileRead, "%lf %lf %lf\n", loc, loc + 1, loc + 2);
            MITK_INFO(fscanOut1 != 3) << ("Error reading file at line: " + QString::number(i)).toStdString();
        }
        MITK_INFO << "Completed input .pts file";
        fclose(pointFileRead);

        fscanOut1 = fscanf(elemFileRead, "%d\n", &nElem);
        MITK_INFO(fscanOut1 == 1) << "File size read correctly";

        std::ofstream outputFileWrite;
        outputFileWrite.open(outputPath.toStdString());
        outputFileWrite << nElem << " 3" << std::endl;

        MITK_INFO << "Beginning input .elem file, simultaneous output";
        for (int i = 0; i < nElem; i++) {
            // 		int* loc = elem_array + 4*i;
            int p1, p2, p3, p4, region;
            fscanOut1 = fscanf(elemFileRead, "Tt %d %d %d %d %d\n", &p1, &p2, &p3, &p4, &region);
            MITK_INFO(fscanOut1 != 5) << ("Error reading file at line: " + QString::number(i)).toStdString();

            // Calculate and output cog
            double x = 0.0, y = 0.0, z = 0.0;

            double *loc = pts_array + 3 * p1;
            x += loc[0];
            y += loc[1];
            z += loc[2];

            loc = pts_array + 3 * p2;
            x += loc[0];
            y += loc[1];
            z += loc[2];

            loc = pts_array + 3 * p3;
            x += loc[0];
            y += loc[1];
            z += loc[2];

            loc = pts_array + 3 * p4;
            x += loc[0];
            y += loc[1];
            z += loc[2];

            x /= 4.0 * 1000;
            y /= 4.0 * 1000;
            z /= 4.0 * 1000;

            outputFileWrite << std::fixed << std::setprecision(6) << x << std::endl;
            outputFileWrite << std::fixed << std::setprecision(6) << y << std::endl;
            outputFileWrite << std::fixed << std::setprecision(6) << z << std::endl;

        }
        MITK_INFO << "Completed input .elem file";

        fclose(elemFileRead);
        outputFileWrite.close();
        MITK_INFO << "Completed input .elem file";


    } else {
        MITK_ERROR(QFileInfo::exists(elemPath)) << ("Could not read file" + elemPath).toStdString();
        MITK_ERROR(QFileInfo::exists(pointPath)) << ("Could not read file" + pointPath).toStdString();
    }
}

void CemrgCommonUtils::RegionMapping(QString bpPath, QString pointPath, QString elemPath, QString outputPath) {
    if (QFileInfo::exists(bpPath) && QFileInfo::exists(pointPath) && QFileInfo::exists(elemPath)) {
        typedef itk::Image<uint8_t, 3> ImageType;
        typedef itk::Index<3> IndexType;
        // ScarImage image(bpPath.toStdString());
        mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(bpPath.toStdString());
        ImageType::Pointer itkInput = ImageType::New();
        ImageType::SpacingType spacing;
        ImageType::RegionType region;
        ImageType::PointType origin;
        ImageType::SizeType size;

        mitk::CastToItkImage(image, itkInput);
        origin = itkInput->GetOrigin();
        spacing = itkInput->GetSpacing();
        region = itkInput->GetLargestPossibleRegion();
        size = region.GetSize();

        double Tx[4][4];
        int max_i = size[0] - 1;
        int max_j = size[1] - 1;
        int max_k = size[2] - 1;
        int min_i = 0;
        int min_j = 0;
        int min_k = 0;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                Tx[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }

        std::ofstream outputFileWrite;
        std::ifstream cogFileRead, elemFileRead;
        cogFileRead.open(pointPath.toStdString());

        int nElemCOG, dim, count;
        double x, y, z;

        cogFileRead >> nElemCOG;
        cogFileRead >> dim;


        outputFileWrite.open(outputPath.toStdString());
        outputFileWrite << nElemCOG << std::endl;

        MITK_INFO << ("Number of elements (COG file):" + QString::number(nElemCOG)).toStdString();
        MITK_INFO << ("Dimension= " + QString::number(dim)).toStdString();

        elemFileRead.open(elemPath.toStdString());
        int nElem;

        elemFileRead >> nElem;
        if (nElem != nElemCOG) {
            MITK_ERROR << "Number of elements in files are not consistent.";
        }

        count = 0;
        char type[2];
        int nodes[4], imregion, newRegion;
        int newRegionCount = 0;

        for (int iElem = 0; iElem < nElemCOG; iElem++) {
            cogFileRead >> x;
            cogFileRead >> y;
            cogFileRead >> z;
            if (cogFileRead.eof()) {
                MITK_WARN << "File ended prematurely";
                break;
            }

            elemFileRead >> type;
            elemFileRead >> nodes[0];
            elemFileRead >> nodes[1];
            elemFileRead >> nodes[2];
            elemFileRead >> nodes[3];
            elemFileRead >> imregion;

            // checking point belonging to imregion (cm2carp/carp_scar_map::inScar())
            double xt, yt, zt;
            xt = Tx[0][0] * x + Tx[0][1] * y + Tx[0][2] * z + Tx[0][3] - origin[0] + spacing[0] / 2;
            yt = Tx[1][0] * x + Tx[1][1] * y + Tx[1][2] * z + Tx[1][3] - origin[1] + spacing[0] / 2;
            zt = Tx[2][0] * x + Tx[2][1] * y + Tx[2][2] * z + Tx[2][3] - origin[2] + spacing[0] / 2;

            int i, j, k;
            i = static_cast<int>(xt / spacing[0]);
            j = static_cast<int>(yt / spacing[1]);
            k = static_cast<int>(zt / spacing[2]);

            if (i < 0 || j < 0 || k < 0) {
                min_i = std::min(i, min_i);
                min_j = std::min(j, min_j);
                min_k = std::min(k, min_k);
                newRegion = 0;
            } else if ((unsigned)i > size[0] - 1 || (unsigned)j > size[1] - 1 || (unsigned)k > size[2] - 1) {
                max_i = std::max(i, max_i);
                max_j = std::max(j, max_j);
                max_k = std::max(k, max_k);
                newRegion = 0;

            } else {
                IndexType index = {{i, j, k}};
                newRegion = itkInput->GetPixel(index);
            }

            if (newRegion != 0) {
                imregion = newRegion;
                newRegionCount++;
            }

            outputFileWrite << type << " ";
            outputFileWrite << nodes[0] << " ";
            outputFileWrite << nodes[1] << " ";
            outputFileWrite << nodes[2] << " ";
            outputFileWrite << nodes[3] << " ";
            outputFileWrite << imregion << std::endl;

            count++;
        }

        outputFileWrite.close();

        MITK_INFO << ("Number of element COG read: " + QString::number(count)).toStdString();
        MITK_INFO << ("Number of new regions determined: " + QString::number(newRegionCount)).toStdString();

        if (min_i < 0 || min_j < 0 || min_k < 0) {
            std::cerr << "WARNING: The elemCOG file falls outside the image bounds! Code assumes that no scar lies in this region." << std::endl;
            std::cerr << "If scar does lie in this region, then you need to pad the image at the start by (in pixel space):" << std::endl;
            std::cerr << "[" << -min_i << ", " << -min_j << ", " << -min_k << "]" << std::endl;
            std::cerr << "And add the following transformation to the TransMatFile (in geometric space):" << std::endl;
            std::cerr << "[" << -min_i * spacing[0] << ", " << -min_j * spacing[1] << ", " << -min_k * spacing[2] << "]" << std::endl;
        }
        if (max_i > int(size[0] - 1) || max_j > int(size[1] - 1) || max_k > int(size[2] - 1)) {
            std::cerr << "WARNING: The elemCOG file falls outside the image bounds! Code assumes that no scar lies in this region." << std::endl;
            std::cerr << "If scar does lie in this region, then you need to pad the image at the end by (in pixel space):" << std::endl;
            std::cerr << "[" << max_i - (size[0] - 1) << ", " << max_j - (size[1] - 1) << ", " << max_k - (size[2] - 1) << "]" << std::endl;
            std::cerr << "No need to change TransMatFile" << std::endl;
        }

    } else {
        MITK_ERROR(QFileInfo::exists(bpPath)) << ("File does not exist: " + bpPath).toStdString();
        MITK_ERROR(QFileInfo::exists(pointPath)) << ("File does not exist: " + pointPath).toStdString();
    }
}

void CemrgCommonUtils::NormaliseFibreFiles(QString fibresPath, QString outputPath) {
    MITK_INFO << "Normalise fibres file";
    std::ifstream ffibres(fibresPath.toStdString());
    std::ofstream fo(outputPath.toStdString());

    int numVect;
    ffibres >> numVect;
    MITK_INFO << ("Number of vectors per line in file: " + QString::number(numVect)).toStdString();

    double x, y, z;
    double norm;
    // prime read
    ffibres >> x;
    ffibres >> y;
    ffibres >> z;
    fo << numVect << std::endl;
    while (!ffibres.eof()) {
        for (int i = 0; i < numVect; i++) {
            norm = sqrt(x * x + y * y + z * z);
            if (norm > 0) {
                fo << std::fixed << std::setprecision(8) << x / norm << " " << y / norm << " " << z / norm << " ";
            } else {
                fo << std::fixed << std::setprecision(8) << x << " " << y << " " << z << " ";
            }
            ffibres >> x;
            ffibres >> y;
            ffibres >> z;
        }
        fo << std::endl;
    }
    ffibres.close();
    fo.close();
}

void CemrgCommonUtils::CarpToVtk(QString elemPath, QString ptsPath, QString outputPath, bool saveRegionlabels) {
    std::ofstream VTKFile;
    std::ifstream ptsFileRead, elemFileRead;
    short int precision = 12;
    short int numColsLookupTable = 1;

    VTKFile.open(outputPath.toStdString());
    MITK_INFO << "Writing vtk file header.";
    VTKFile << "# vtk DataFile Version 4.0" << std::endl;
    VTKFile << "vtk output" << std::endl;
    VTKFile << "ASCII" << std::endl;
    VTKFile << "DATASET UNSTRUCTURED_GRID" << std::endl;

    int nElem, nPts;
    ptsFileRead.open(ptsPath.toStdString());
    ptsFileRead >> nPts;

    MITK_INFO << "Setting geometry - Points";
    VTKFile << "POINTS " << nPts << " float" << std::endl;
    double x, y, z;
    for (int ix = 0; ix < nPts; ix++) {
        ptsFileRead >> x;
        ptsFileRead >> y;
        ptsFileRead >> z;

        VTKFile << std::setprecision(precision) << x << " " << y << " " << z << std::endl;
    }
    ptsFileRead.close();

    MITK_INFO << "Setting geometry - Tetrahedral elements";
    elemFileRead.open(elemPath.toStdString());
    elemFileRead >> nElem;
    std::string type;
    int p0, p1, p2, p3;
    std::vector<double> regionVector(nElem);
    VTKFile << "CELLS " << nElem << " " << (4 + 1) * nElem << std::endl;
    for (int ix = 0; ix < nElem; ix++) {
        elemFileRead >> type;
        elemFileRead >> p0;
        elemFileRead >> p1;
        elemFileRead >> p2;
        if (type.compare("Tr") == 0) {
            elemFileRead >> p3;
        }
        elemFileRead >> regionVector[ix];

        if (type.compare("Tr") == 0) {
            VTKFile << "3 " << p0 << " " << p1 << " " << p2 << std::endl;
        } else {
            VTKFile << "4 " << p0 << " " << p1 << " " << p2 << " " << p3 << std::endl;
        }
    }

    VTKFile << "CELL_TYPES " << nElem << " ";
    for (int ix = 0; ix < nElem; ix++) {
        VTKFile << "10"; // type for tetrahedral mesh
        if (((1 + ix) % numColsLookupTable) && (ix < (nElem - 1))) {
            VTKFile << " ";
        } else {
            VTKFile << std::endl;
        }
    }
    elemFileRead.close();
    VTKFile.close();
    if (saveRegionlabels) {
        AppendScalarFieldToVtk(outputPath, "region_labels", "CELL", regionVector);
    }
}

void CemrgCommonUtils::RectifyFileValues(QString pathToFile, double minVal, double maxVal) {
    QFileInfo fi(pathToFile);
    QString copyName = fi.absolutePath() + "/" + fi.baseName() + "_copy." + fi.completeSuffix();
    MITK_INFO << "Copying path name";
    QFile::rename(pathToFile, copyName);

    std::ifstream readInFile(copyName.toStdString());
    std::ofstream writeOutFile(pathToFile.toStdString());

    double valueOut;
    int precision = 16;
    int count = 0;
    writeOutFile << std::scientific;
    MITK_INFO << "Rectifying file.";
    while (!readInFile.eof()) {
        double valueIn = -1;
        readInFile >> valueIn;
        if (valueIn > maxVal) {
            valueOut = maxVal;
        } else if (valueIn < minVal) {
            valueOut = minVal;
        } else {
            valueOut = valueIn;
        }
        if (valueIn != -1) {
            writeOutFile << std::setprecision(precision) << valueOut << std::endl;
            count++;
        }
    }
    MITK_INFO << ("Finished rectifying file with :" + QString::number(count) + " points.").toStdString();
    readInFile.close();
    writeOutFile.close();

    QFile::remove(copyName);
}

int CemrgCommonUtils::GetTotalFromCarpFile(QString pathToFile, bool totalAtTop) {
    int total = -1;
    std::ifstream fi(pathToFile.toStdString());
    if (totalAtTop) {
        fi >> total;
    } else {
        int count;
        double value;
        count = 0;
        while (!fi.eof()) {
            value = -1;
            fi >> value;
            count++;
        }
        total = count;
        total -= (value == -1) ? 1 : 0; // checks if last line in file is empty
    }
    fi.close();
    return total;
}

std::vector<double> CemrgCommonUtils::ReadScalarField(QString pathToFile) {
    std::ifstream fi(pathToFile.toStdString());
    int n = CemrgCommonUtils::GetTotalFromCarpFile(pathToFile, false);
    std::vector<double> field(n, 0.0);

    for (int ix = 0; ix < n; ix++) {
        if (fi.eof()) {
            MITK_INFO << "File finished prematurely.";
            break;
        }
        fi >> field[ix];
    }
    fi.close();

    return field;
}

std::vector<double> CemrgCommonUtils::ReadVectorField(QString pathToFile, bool totalAtTop) {
    std::ifstream fi(pathToFile.toStdString());
    int n = CemrgCommonUtils::GetTotalFromCarpFile(pathToFile, totalAtTop);
    std::vector<double> field((n*3), 0.0);

    for (int ix = 0; ix < n; ix++) {
        if (fi.eof()) {
            MITK_INFO << "File finished prematurely.";
            break;
        }
        fi >> field[3*ix + 0];
        fi >> field[3*ix + 1];
        fi >> field[3*ix + 2];
    }
    fi.close();

    return field;
}

void CemrgCommonUtils::AppendScalarFieldToVtk(QString vtkPath, QString fieldName, QString typeData, std::vector<double> field, bool setHeader) {
    std::ofstream VTKFile;
    short int precision = 12;
    short int numColsLookupTable = 1;

    VTKFile.open(vtkPath.toStdString(), std::ios_base::app);
    int fieldSize = field.size();

    if (setHeader) {
        MITK_INFO << "Setting POINT_DATA header.";
        VTKFile << typeData.toStdString() << "_DATA " << fieldSize << std::endl;
    }

    MITK_INFO << ("Appending scalar field <<" + fieldName + ">> to VTK file.").toStdString();
    VTKFile << "SCALARS " << fieldName.toStdString() << " FLOAT " << numColsLookupTable << " " << std::endl;
    VTKFile << "LOOKUP_TABLE default " << std::endl;

    for (int ix = 0; ix < fieldSize; ix++) {
        VTKFile << std::setprecision(precision) << field.at(ix);
        if (((1 + ix) % numColsLookupTable) && (ix < fieldSize - 1)) {
            VTKFile << " ";
        } else {
            VTKFile << std::endl;
        }
    }
    VTKFile.close();
}

void CemrgCommonUtils::AppendVectorFieldToVtk(QString vtkPath, QString fieldName, QString dataType, std::vector<double> field, bool setHeader) {
    std::ofstream VTKFile;
    short int precision = 12;
    // short int numColsLookupTable=1;

    VTKFile.open(vtkPath.toStdString(), std::ios_base::app);
    int nElem = field.size() / 3;

    if (setHeader) {
        MITK_INFO << "Setting CELL_DATA header.";
        VTKFile << dataType.toStdString() << "_DATA " << nElem << std::endl;
    }

    MITK_INFO << ("Appending vector field <<" + fieldName + ">> to VTK file.").toStdString();
    VTKFile << "VECTORS " << fieldName.toStdString() << " FLOAT " << std::endl;

    for (int ix = 0; ix < nElem; ix++) {
        double x, y, z;
        x = field[ix + 0 * nElem];
        y = field[ix + 1 * nElem];
        z = field[ix + 2 * nElem];

        VTKFile << std::setprecision(precision) << x << " " << y << " " << z << std::endl;
    }

    VTKFile.close();
}

void CemrgCommonUtils::VtkScalarToFile(QString vtkPath, QString outPath, QString fieldName, bool isElem){
    mitk::Surface::Pointer surface = mitk::IOUtil::Load<mitk::Surface>(vtkPath.toStdString());
    vtkFloatArray *scalars = vtkFloatArray::New();
    vtkIdType numObjects;

    std::cout << "fieldName: " << fieldName.toStdString() << '\n';

    if (isElem){
        // surface->GetVtkPolyData()->GetCellData()->SetActiveScalars(fieldName.toStdString().c_str());
        scalars = vtkFloatArray::SafeDownCast(surface->GetVtkPolyData()->GetCellData()->GetScalars());
        numObjects = surface->GetVtkPolyData()->GetNumberOfCells();
    } else{
        // surface->GetVtkPolyData()->GetCellData()->SetActiveScalars(fieldName.toStdString().c_str());
        scalars = vtkFloatArray::SafeDownCast(surface->GetVtkPolyData()->GetCellData()->GetScalars());
        numObjects = surface->GetVtkPolyData()->GetNumberOfCells();
    }

    std::ofstream fo(outPath.toStdString());


    for (vtkIdType ix=0;ix<numObjects;ix++) {
        double s = scalars->GetTuple1(ix);
        fo << std::setprecision(12) << s;
        if(ix<numObjects-1){
            fo << std::endl;
        }
    }
}

void CemrgCommonUtils::VtkPointScalarToFile(QString vtkPath, QString outPath, QString fieldName){
    VtkScalarToFile(vtkPath, outPath, fieldName, false);
}

void CemrgCommonUtils::VtkCellScalarToFile(QString vtkPath, QString outPath, QString fieldName){
    VtkScalarToFile(vtkPath, outPath, fieldName, true);
}
