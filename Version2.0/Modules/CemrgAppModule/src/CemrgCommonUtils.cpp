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
#include <itkBSplineInterpolateImageFunction.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkResampleImageFilter.h>
#include <itkOrientImageFilter.h>

//VTK
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

//Qmitk
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

//Qt
#include <QMessageBox>
#include <QString>
#include <QFile>
#include <QFileInfo>
#include <QTextStream>
#include "CemrgCommonUtils.h"


mitk::DataNode::Pointer CemrgCommonUtils::imageNode;
mitk::DataNode::Pointer CemrgCommonUtils::cuttingNode;
mitk::Image::Pointer CemrgCommonUtils::imageToCut;
mitk::BoundingObject::Pointer CemrgCommonUtils::cuttingCube;

mitk::Image::Pointer CemrgCommonUtils::CropImage() {

    //Test input objects
    if (imageToCut.IsNull() || cuttingCube.IsNull()) {
        return NULL;
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
        return NULL;
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

    typedef itk::Image<short,3> ImageType;
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
    return image;
}

mitk::Image::Pointer CemrgCommonUtils::IsoImageResampleReorient(mitk::Image::Pointer image, bool resample, bool reorientToRAI) {

    MITK_INFO(resample) << "Resampling image to be isometric.";
    MITK_INFO(reorientToRAI) << "Doing a reorientation to RAI.";

    typedef itk::Image<short,3> ImageType;
    typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleImageFilterType;
    typedef itk::BSplineInterpolateImageFunction<ImageType, double, double> BSplineInterpolatorType;
    ImageType::Pointer itkInputImage = ImageType::New();
    ImageType::Pointer resampleOutput = ImageType::New();
    ImageType::Pointer outputImage = ImageType::New();
    mitk::CastToItkImage(image, itkInputImage);

    if (resample) {

        ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
        BSplineInterpolatorType::Pointer binterp = BSplineInterpolatorType::New();
        binterp->SetSplineOrder(3);
        resampler->SetInterpolator(binterp);
        resampler->SetInput(itkInputImage);
        resampler->SetOutputOrigin(itkInputImage->GetOrigin());
        ImageType::SizeType input_size = itkInputImage->GetLargestPossibleRegion().GetSize();
        ImageType::SpacingType input_spacing = itkInputImage->GetSpacing();
        ImageType::SizeType output_size;
        ImageType::SpacingType output_spacing;
        output_size[0] = input_size[0] * (input_spacing[0] / 1.0);
        output_size[1] = input_size[1] * (input_spacing[1] / 1.0);
        output_size[2] = input_size[2] * (input_spacing[2] / 1.0);
        output_spacing [0] = 1.0;
        output_spacing [1] = 1.0;
        output_spacing [2] = 1.0;
        resampler->SetSize(output_size);
        resampler->SetOutputSpacing(output_spacing);
        resampler->SetOutputDirection(itkInputImage->GetDirection());
        resampler->UpdateLargestPossibleRegion();
        resampleOutput = resampler->GetOutput();

    } else {
        resampleOutput = itkInputImage;
    }//_if

    if (reorientToRAI) {

        typedef itk::OrientImageFilter<ImageType,ImageType> OrientImageFilterType;
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

mitk::Image::Pointer CemrgCommonUtils::IsoImageResampleReorient(QString imPath, bool resample,  bool reorientToRAI) {

    return CemrgCommonUtils::IsoImageResampleReorient(mitk::IOUtil::Load<mitk::Image>(imPath.toStdString()), resample, reorientToRAI);
}

bool CemrgCommonUtils::ConvertToNifti(mitk::BaseData::Pointer oneNode, QString path2file, bool resample, bool reorient) {

    bool successful = false;

    if (oneNode) {

        mitk::Image::Pointer image = dynamic_cast<mitk::Image*>(oneNode.GetPointer());
        if (image) { //Test if this data item is an image
            image = CemrgCommonUtils::IsoImageResampleReorient(image, resample, reorient);
            mitk::IOUtil::Save(image, path2file.toStdString());
            successful = true;
        } else{
            MITK_INFO << "[...] Problem casting node data to image";
        }//_if

    } else{
        MITK_INFO << "[...] Problem with node";
    }//_if

    return successful;
}

void CemrgCommonUtils::RoundPixelValues(QString pathToImage, QString outputPath){
    QFileInfo fi(pathToImage);
    if(fi.exists()){
        using ImageType = itk::Image<double, 3>;
        using IteratorType = itk::ImageRegionIteratorWithIndex<ImageType>;

        ImageType::Pointer im = ImageType::New();
        mitk::CastToItkImage(mitk::IOUtil::Load<mitk::Image>(pathToImage.toStdString()), im);

        IteratorType imIter(im, im->GetLargestPossibleRegion());

        imIter.GoToBegin();
        while(!imIter.IsAtEnd()){
            double pixelValue = imIter.Get();
            imIter.Set(std::round(pixelValue));

            ++imIter;
        }

        QString writingPath = (outputPath.isEmpty()) ? pathToImage : outputPath;

        mitk::Image::Pointer outputImg = mitk::Image::New();
        mitk::CastToMitkImage(im, outputImg);

        MITK_INFO(outputPath.isEmpty()) << ("Overwriting: " + pathToImage).toStdString();
        mitk::IOUtil::Save(outputImg, writingPath.toStdString());

    } else{
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

mitk::Surface::Pointer CemrgCommonUtils::ExtractSurfaceFromSegmentation(mitk::Image::Pointer image, double thresh, double blur, double smooth, double decimation){
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

    if(!saveToPath.isEmpty()){
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

void CemrgCommonUtils::FlipXYPlane(mitk::Surface::Pointer surf, QString dir, QString vtkname){

    //Prepare points for MITK visualisation - (CemrgCommonUtils::LoadVTKMesh)
    vtkSmartPointer<vtkPolyData> pd = surf->GetVtkPolyData();
    for (int ix=0; ix<pd->GetNumberOfPoints(); ix++) {
        double* point = pd->GetPoint(ix);
        point[0] = -point[0];
        point[1] = -point[1];
        pd->GetPoints()->SetPoint(ix, point);
    }

    if(!vtkname.isEmpty()){
        vtkname += (!vtkname.contains(".vtk")) ? ".vtk" : "";
        QString path = dir + mitk::IOUtil::GetDirectorySeparator()+vtkname;
        mitk::IOUtil::Save(surf, path.toStdString());
    }
}

QString CemrgCommonUtils::M3dlibParamFileGenerator(QString dir, QString filename, QString thicknessCalc) {

    QString path2file = dir + mitk::IOUtil::GetDirectorySeparator() + filename;
    QFile fi(path2file);

    if (thicknessCalc.compare("0", Qt::CaseSensitive)!=0 && thicknessCalc.compare("1", Qt::CaseSensitive)!=0) {
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

void CemrgCommonUtils::ConvertToCarto(
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
    ofstream cartoFile;
    cartoFile.open(outputPath);

    //Header
    cartoFile << "# vtk DataFile Version 3.0\n";
    cartoFile << "PatientData Anon Anon 00000000\n";
    cartoFile << "ASCII\n";
    cartoFile << "DATASET POLYDATA\n";

    //Points
    cartoFile << "POINTS\t" << pd->GetNumberOfPoints() << "\tfloat\n";
    for (int i=0; i<pd->GetNumberOfPoints(); i++) {
        double* point = pd->GetPoint(i);
        cartoFile << point[0] << " " << point[1] << " " << point[2] << "\n";
    }
    cartoFile << "\n";

    //Cells
    cartoFile << "POLYGONS\t";
    cartoFile << pd->GetNumberOfCells() << "\t";
    cartoFile << pd->GetNumberOfCells()*4 << "\n";
    for (int i=0; i<pd->GetNumberOfCells(); i++) {
        vtkCell* cell = pd->GetCell(i);
        vtkIdList* list = cell->GetPointIds();
        cartoFile << "3";
        for (int j=0; j<list->GetNumberOfIds(); j++)
            cartoFile << " " << list->GetId(j);
        cartoFile << "\n";
    }

    //Point data
    vtkSmartPointer<vtkFloatArray> pointData = vtkSmartPointer<vtkFloatArray>::New();
    try {
        pointData = vtkFloatArray::SafeDownCast(pd->GetPointData()->GetScalars());
    } catch (...) {
        MITK_WARN << "Storing point data failed! Check your input";
        return;
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
            for (int i=0; i<pointData->GetNumberOfTuples(); i++) {

                //Get scalar raw value
                double value = static_cast<double>(pointData->GetTuple1(i));

                //Colouring
                if (discreteScheme) {
                    if (methodType == 1) {
                        if (value < (meanBP * thresholds.at(0))) value = 0.0;
                        else if (thresholds.size() == 2 && value < (meanBP * thresholds.at(1))) value = 0.5;
                        else value = 1.0;
                    } else {
                        if (value < (meanBP + thresholds.at(0)*stdvBP)) value = 0.0;
                        else if (thresholds.size() == 2 && value < (meanBP + thresholds.at(1)*stdvBP)) value = 0.5;
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

            for (int i=0; pointData->GetNumberOfComponents(); i++) {

                cartoFile << "SCALARS " << "scalars" << i << " float\n";
                cartoFile << "LOOKUP_TABLE lookup_table\n";
                for (int j=0; j<pointData->GetNumberOfTuples(); j++)
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
    lut->AddRGBPoint((numCols-1.0)/2.0, 0.94, 0.47, 0.12);
    lut->AddRGBPoint((numCols-1.0), 0.90, 0.11, 0.14);
    lut->SetScaleToLinear();
    for (int i=0; i<numCols; i++) {
        cartoFile << lut->GetColor(i)[0] << " ";
        cartoFile << lut->GetColor(i)[1] << " ";
        cartoFile << lut->GetColor(i)[2] << " ";
        cartoFile << "1.0" <<"\n";
    }//_for

    cartoFile.close();
}

void CemrgCommonUtils::MotionTrackingReport(QString directory, int timePoints) {

    for (int tS=0; tS<timePoints; tS++) {

        //Image
        QString path = directory + mitk::IOUtil::GetDirectorySeparator() + "dcm-" + QString::number(tS) + ".nii";
        mitk::Image::Pointer img3D = mitk::IOUtil::Load<mitk::Image>(path.toStdString());
        int* extent = img3D->GetVtkImageData()->GetExtent();
        mitk::Vector3D spacing = img3D->GetGeometry()->GetSpacing();
        vtkSmartPointer<vtkMatrix4x4> direction = img3D->GetGeometry()->GetVtkMatrix();

        //Mesh
        path = directory + mitk::IOUtil::GetDirectorySeparator() + "Model-" + QString::number(tS) + ".vtk";
        mitk::Surface::Pointer sur3D = CemrgCommonUtils::LoadVTKMesh(path.toStdString());

        //Window
        int xSliceMin = extent[0];
        int xSliceMax = extent[1];
        int ySliceMin = extent[2];
        int ySliceMax = extent[3];
        //int zSliceMin = extent[4];
        int zSliceMax = extent[5];
        double xmins[8] = {0.00,0.25,0.50,0.75,0.00,0.25,0.50,0.75};
        double xmaxs[8] = {0.25,0.50,0.75,1.00,0.25,0.50,0.75,1.00};
        double ymins[8] = {0.00,0.00,0.00,0.00,0.50,0.50,0.50,0.50};
        double ymaxs[8] = {0.50,0.50,0.50,0.50,1.00,1.00,1.00,1.00};
        vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->SetAlphaBitPlanes(1);
        renderWindow->SetSize(500,500);

        for (int view=0; view<8; view++) {

            //Setup views
            int zSlice = zSliceMax - view * floor(zSliceMax/8);
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
            plane->SetOrigin(0,0,zPlane);
            plane->SetNormal(0,0,1);
            vtkSmartPointer<vtkCutter> cutter = vtkSmartPointer<vtkCutter>::New();
            cutter->SetCutFunction(plane);
            cutter->SetInputData(pd);
            cutter->Update();
            vtkSmartPointer<vtkPolyDataMapper> mapMesh = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapMesh->SetInputConnection(cutter->GetOutputPort());
            mapMesh->SetScalarModeToUsePointData();
            mapMesh->SetScalarVisibility(1);
            mapMesh->SetScalarRange(1,4);
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
            renderer->GetActiveCamera()->SetParallelScale(.5*imgActor->GetBounds()[1]);
            renderWindow->AddRenderer(renderer);
            renderer->SetViewport(xmins[view],ymins[view],xmaxs[view],ymaxs[view]);
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
        writer->SetFileName((directory.toStdString() + mitk::IOUtil::GetDirectorySeparator() + "dcm-" + QString::number(tS).toStdString() + ".png").c_str());
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
    } else{ // pointdata
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
