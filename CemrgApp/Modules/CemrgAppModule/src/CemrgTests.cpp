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
 * CEMRG TESTS
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#include <mitkIOUtil.h>
#include <mitkPointSet.h>
#include <mitkDataNode.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkUnstructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetSurfaceFilter.h>
#include <CemrgMeasure.h>
#include <CemrgStrains.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkCubeSource.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkStructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPointData.h>
#include <vtkMetaImageWriter.h>
#include <vtkPointLocator.h>
#include <mitkBoundingObjectCutter.h>
#include <itkConformalFlatteningMeshFilter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkExtractUnstructuredGrid.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <vtkPolyDataWriter.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkPolygon.h>
#include <vtkDelaunay2D.h>
#include <vtkLookupTable.h>
#include <vtkCellData.h>


void autoNIIconvert() {

    std::string frames;
    std::ifstream file0("/home/or15/Desktop/MRI/process/frames.txt");
    if (file0.is_open())
        getline(file0, frames);
    file0.close();

    for (int i = 1; i <= QString::fromStdString(frames).toInt(); i++) {

        unsigned int slices = 0;
        mitk::Image::Pointer img3D = mitk::Image::New();

        std::string line;
        std::ifstream file1("/home/or15/Desktop/MRI/process/" + QString::number(i).toStdString() + "/path.txt");
        if (file1.is_open()) {
            getline(file1, line);
            while (getline(file1, line)) {
                if (line == "path.txt")
                    break;
                slices++;
            }
        }
        file1.close();
        std::ifstream file2("/home/or15/Desktop/MRI/process/" + QString::number(i).toStdString() + "/path.txt");
        if (file2.is_open()) {
            bool first = true;
            int ctr = 0;
            getline(file2, line);
            while (getline(file2, line)) {
                if (line == "path.txt")
                    break;
                mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>("/home/or15/Desktop/MRI/process/" + QString::number(i).toStdString() + "/" + line);
                if (first) {
                    mitk::ImageDescriptor::Pointer dsc = image->GetImageDescriptor();
                    img3D->Initialize(dsc->GetChannelDescriptor(0).GetPixelType(), *image->GetGeometry(), slices, 1);
                    first = false;
                }
                img3D->SetVolume(mitk::ImageReadAccessor(image).GetData(), 0, ctr);
                ctr++;
            }
        }
        file2.close();

        mitk::IOUtil::Save(img3D, "/home/or15/Desktop/MRI/process/dcm-" + QString::number(i - 1).toStdString() + ".nii");
    }//for
}

void CreateSyntImage() {
    int ctr = -5;

    //Prepare empty image
    vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();
    double spacing[3] = {0.83, 0.83, 0.83};
    whiteImage->SetSpacing(spacing);
    int dimensions[3] = {25, 25, 25};
    whiteImage->SetDimensions(dimensions);
    int extents[3] = {dimensions[0] - 1, dimensions[1] - 1, dimensions[2] - 1};
    whiteImage->SetExtent(0, extents[0], 0, extents[1], 0, extents[2]);
    double origin[3] = {0, 0, 0};
    whiteImage->SetOrigin(origin);
    whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

    //Colouring
    for (vtkIdType i = 0; i < whiteImage->GetNumberOfPoints(); ++i)
        whiteImage->GetPointData()->GetScalars()->SetTuple1(i, 0);

    //Tracking test
    ctr += 5;
    for (int x = 10; x < 15; x++) {
        for (int y = 10 + ctr; y < 15 + ctr; y++) {
            for (int z = 10 + ctr; z < 15 + ctr; z++) {
                vtkIdType idx = x + dimensions[0] * (y + dimensions[1] * z);
                whiteImage->GetPointData()->GetScalars()->SetTuple1(idx, 50);
            }
        }
    }

    //    for(int x=10; x<15; x++) {
    //        for(int y=10; y<15; y++) {
    //            for(int z=10; z<15; z++) {
    //                vtkIdType idx = x + dimensions[0] * (y + dimensions[1] * z);
    //                whiteImage->GetPointData()->GetScalars()->SetTuple1(idx, 50);
    //            }
    //        }
    //    }
    //    for(int x=9; x<16; x++) {
    //        for(int y=9; y<16; y++) {
    //            for(int z=9; z<16; z++) {
    //                vtkIdType idx = x + dimensions[0] * (y + dimensions[1] * z);
    //                if (whiteImage->GetPointData()->GetScalars()->GetTuple1(idx) != 50)
    //                    whiteImage->GetPointData()->GetScalars()->SetTuple1(idx, 100);
    //                if (x==12 && y==12 && z==9)
    //                    whiteImage->GetPointData()->GetScalars()->SetTuple1(idx, 200);
    //            }
    //        }
    //    }

    mitk::Image::Pointer image = mitk::Image::New();
    image->Initialize(whiteImage);
    image->SetVolume(whiteImage->GetScalarPointer());
    mitk::IOUtil::Save(image, "/home/or15/Desktop/test.nii");
}

void ColourMesh() {

    for (int m = 0; m < 10; m++) {
        mitk::Surface::Pointer mesh;
        mesh = mitk::IOUtil::Load<mitk::Surface>("/home/or15/Downloads/SqzMeshes/Mesh" + QString::number(m).toStdString() + ".vtk");
        vtkSmartPointer<vtkUnsignedCharArray> segmentColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
        segmentColors->SetNumberOfComponents(3);
        segmentColors->SetNumberOfTuples(mesh->GetVtkPolyData()->GetNumberOfCells());
        double min = mesh->GetVtkPolyData()->GetScalarRange()[0];
        double max = mesh->GetVtkPolyData()->GetScalarRange()[1];

        //Create the color map
        vtkSmartPointer<vtkLookupTable> colorLookupTable = vtkSmartPointer<vtkLookupTable>::New();
        colorLookupTable->SetTableRange(min, max);
        colorLookupTable->Build();
        vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
        colors->SetNumberOfComponents(3);
        colors->SetName("Colors");
        colors->SetNumberOfTuples(mesh->GetVtkPolyData()->GetNumberOfCells());
        for (int i = 0; i < mesh->GetVtkPolyData()->GetNumberOfCells(); i++) {
            double sqz = mesh->GetVtkPolyData()->GetCellData()->GetScalars()->GetTuple1(i);
            double dcolor[3];
            colorLookupTable->GetColor(sqz, dcolor);
            unsigned char color[3];
            for (unsigned int j = 0; j < 3; j++)
                color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
            colors->InsertTuple3(i, color[0], color[1], color[2]);
        }//_for

        mesh->GetVtkPolyData()->GetCellData()->SetScalars(colors);
        mitk::IOUtil::Save(mesh, "/home/or15/Downloads/SqzMeshes/MeshColoured" + QString::number(m).toStdString() + ".vtk");
    }//_for
}

void binariseImage() {

    typedef itk::Image<double, 3> ImageType;
    mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>("/home/or15/Downloads/wallThicknessNRRD/leftAtrialWall900101_UPsampled.nrrd");
    ImageType::Pointer imageITK = ImageType::New();
    CastToItkImage(image, imageITK);

    typedef itk::ImageRegionIteratorWithIndex<ImageType> ItType;
    ItType it(imageITK, imageITK->GetRequestedRegion());
    for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
        if (it.Get() > 0.5 && it.Get() != 0)// && it.Get() != 3 && it.Get() != 2)
            it.Set(1);
        else if (it.Get() <= 0.5 && it.Get() != 0)// && it.Get() != 3 && it.Get() != 2)
            it.Set(0);
    }//_for

    mitk::IOUtil::Save(mitk::ImportItkImage(imageITK), "/home/or15/Downloads/wallThicknessNRRD/leftAtrialWall900101_UPsampled_mask.nrrd");
}

#include <itkSimpleContourExtractorImageFilter.h>

void MeshError() {

    for (int i = 1; i <= 6; i++) {
        //if (i==4) continue;
        QString pathImage = "/home/or15/Desktop/Proj/RZ/StrainsWork/MeshError/N" + QString::number(i) + "/segmentation.nii";
        typedef itk::Image<unsigned char, 3> ImageType;
        typedef itk::ImageRegionIteratorWithIndex<ImageType> ItType;
        ImageType::Pointer segItkImage = ImageType::New();
        CastToItkImage(mitk::IOUtil::Load<mitk::Image>(pathImage.toStdString()), segItkImage);
        typedef itk::SimpleContourExtractorImageFilter<ImageType, ImageType> SimpleContourExtractorImageFilterType;
        SimpleContourExtractorImageFilterType::Pointer contourFilter = SimpleContourExtractorImageFilterType::New();
        contourFilter->SetInput(segItkImage);
        contourFilter->SetInputBackgroundValue(0);
        contourFilter->SetInputForegroundValue(1);
        contourFilter->Update();
        ImageType::Pointer ctrItkImage = contourFilter->GetOutput();
        QString pathMesh = "/home/or15/Desktop/Proj/RZ/StrainsWork/MeshError/N" + QString::number(i) + "/segmentation.vtk";
        mitk::Surface::Pointer surface = mitk::IOUtil::Load<mitk::Surface>(pathMesh.toStdString());
        vtkSmartPointer<vtkPolyData> pd = surface->GetVtkPolyData();
        for (int j = 0; j < pd->GetNumberOfPoints(); j++) {
            double* point = pd->GetPoint(j);
            point[0] = -point[0];
            point[1] = -point[1];
            pd->GetPoints()->SetPoint(j, point);
        }//_for
        vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
        pointLocator->SetDataSet(pd);
        pointLocator->BuildLocator();

        std::vector<double> distances;
        ItType it(ctrItkImage, ctrItkImage->GetRequestedRegion());
        for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
            if ((int)it.Get() != 0) {
                double imgPoint[3];
                ImageType::PointType point;
                segItkImage->TransformIndexToPhysicalPoint(it.GetIndex(), point);
                imgPoint[0] = point[0];
                imgPoint[1] = point[1];
                imgPoint[2] = point[2];
                vtkIdType id = pointLocator->FindClosestPoint(imgPoint);
                double* mshPoint;
                mshPoint = pd->GetPoint(id);
                double x_d = imgPoint[0] - mshPoint[0];
                double y_d = imgPoint[1] - mshPoint[1];
                double z_d = imgPoint[2] - mshPoint[2];
                distances.push_back(sqrt(pow(x_d, 2) + pow(y_d, 2) + pow(z_d, 2)));
            }//_if
        }//_for

        std::ofstream file;
        QString pathOutput = "/home/or15/Desktop/Proj/RZ/StrainsWork/MeshError/N" + QString::number(i) + "/distances.csv";
        file.open(pathOutput.toStdString());
        for (size_t j = 0; j < distances.size(); j++) {
            file << distances.at(j);
            if (j == distances.size() - 1) file << endl;
            else file << ",";
        }
        qDebug() << "Case" << i;
    }//_i_18
}

void polygonCutter() {

    //Create the geometry of a point (the coordinate)
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    const float p0[3] = {29, -151, -149};
    const float p1[3] = {59, -142, -142};
    const float p2[3] = {29, -134, -149};
    const float p3[3] = {29, -136, -160};
    const float p4[3] = {39, -146, -168};
    const float p5[3] = {29, -154, -160};

    points->InsertNextPoint(p0);
    points->InsertNextPoint(p1);
    points->InsertNextPoint(p2);
    points->InsertNextPoint(p3);
    points->InsertNextPoint(p4);
    points->InsertNextPoint(p5);

    // Create the polygon
    vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
    polygon->GetPointIds()->SetNumberOfIds(6); //make a quad
    polygon->GetPointIds()->SetId(0, 0);
    polygon->GetPointIds()->SetId(1, 1);
    polygon->GetPointIds()->SetId(2, 2);
    polygon->GetPointIds()->SetId(3, 3);
    polygon->GetPointIds()->SetId(4, 4);
    polygon->GetPointIds()->SetId(5, 5);
    vtkSmartPointer<vtkCellArray> polygons = vtkSmartPointer<vtkCellArray>::New();
    polygons->InsertNextCell(polygon);
    vtkSmartPointer<vtkPolyData> polygonPolyData = vtkSmartPointer<vtkPolyData>::New();
    polygonPolyData->SetPoints(points);
    polygonPolyData->SetPolys(polygons);

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInputData(polygonPolyData);
    writer->SetFileName("/home/or15/Desktop/test.vtk");
    writer->Write();
}

/**** The MeshUtil class was removed from the MITK ****
void flattening() {

    mitk::Surface::Pointer surface = mitk::IOUtil::Load<mitk::Surface>("/home/or15/Downloads/LGECART30MIN-Scar.vtk");

    typedef itk::Mesh<double, 3> MeshType;
    typedef itk::ConformalFlatteningMeshFilter<MeshType, MeshType>  FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(MeshUtil<MeshType>::MeshFromSurface(surface));
    filter->SetPolarCellIdentifier(0);
    filter->MapToSphere();
    //filter->MapToPlane();
    //filter->SetScale(scale);
    try {
        filter->Update();
    } catch (itk::ExceptionObject&) {
    }
    MeshType::Pointer mesh = filter->GetOutput();

    surface->SetVtkPolyData(MeshUtil<MeshType>::MeshToPolyData(mesh));
    mitk::IOUtil::Save(surface, "/home/or15/Downloads/test.vtk");
}
*/

void RRcalcsAuto() {

    std::vector<double> values;
    std::unique_ptr<CemrgMeasure> rr(new CemrgMeasure());
    CemrgMeasure::Points points;
    std::string lineD;

    //1
    std::ifstream dirc("/home/or15/Desktop/Proj/RZ/StrainsWork/TrackTSFFD/Transforms/dirPaths.txt");
    std::vector<int> DS = {2, 5, 6, 8, 9, 10, 12, 13, 16, 17};

    if (dirc.is_open()) {
        for (int ds = 0; ds < 10; ds++) {
            dirc.clear();
            dirc.seekg(0, ios::beg);
            while (getline(dirc, lineD)) {

                //2
                std::string prePath;
                auto const posSL2 = lineD.find_last_of('/');
                auto const posSL1 = lineD.find_last_of('/', posSL2 - 10);
                prePath = std::string(lineD.substr(0, posSL1 + 1).c_str()) + "GridSearchBE6SW/Dataset" + std::to_string(DS.at(ds)) + lineD.substr(posSL2);

                //3
                for (int i = 0; i < 1; i++) {
                    for (int j = 0; j < 101; j++) {

                        //4
                        std::stringstream stream;
                        //stream << "/Perm-RG" << std::fixed << std::setprecision(1) << j/10.0;
                        stream << "/Perm-BE" << std::fixed << std::setprecision(2) << i / 100.0 << "-SW" << std::fixed << std::setprecision(2) << j / 100.0;
                        std::string permStr = stream.str();

                        for (int frame = 0; frame < 10; frame++) {
                            QString path = QString::fromStdString(prePath) + QString::fromStdString(permStr);
                            qDebug() << path;
                            points = rr->Deconvert(path, frame);
                            if (QString::fromStdString(prePath).endsWith("PP", Qt::CaseInsensitive))
                                values.push_back(rr->CalcPerimeter(points));
                            else if (QString::fromStdString(prePath).endsWith("MAA", Qt::CaseInsensitive))
                                values.push_back(rr->CalcArea(points));
                            else
                                values.push_back(rr->CalcDistance(points));
                        }
                        std::ofstream file;
                        QString path = QString::fromStdString(prePath) + QString::fromStdString(permStr) + ".csv";
                        file.open(path.toStdString());
                        for (size_t z = 0; z < values.size(); z++) {
                            file << values.at(z);
                            if (z == values.size() - 1) file << endl;
                            else file << ",";
                        }
                        values.clear();
                        points.clear();
                        file.close();
                    }
                }
            }
        }
    }
}

void RRcnvrtAuto() {

    //Converts user MPS to input.vtk for later calcs
    std::string lineF, lineD;
    std::ifstream file("/home/or15/Desktop/Proj/Tom/Points/Tools/filPaths.txt");
    std::ifstream dirc("/home/or15/Desktop/Proj/Tom/Points/Tools/dirPaths.txt");
    if (file.is_open() && dirc.is_open()) {
        while (getline(file, lineF)) {

            getline(dirc, lineD);
            mitk::PointSet::Pointer MIPS = mitk::IOUtil::Load<mitk::PointSet>(lineF);
            mitk::DataNode::Pointer node = mitk::DataNode::New();
            node->SetData(MIPS);
            std::unique_ptr<CemrgMeasure> rr(new CemrgMeasure());
            rr->Convert(QString::fromStdString(lineD), node);
        }
    }//_if
    file.close();
    dirc.close();
}

/**
 * @brief MeshReader try flipping mesh to overlap segmentation
 */
void MeshReader() {

    std::string inputFilename = "/home/or15/Desktop/Marina/BiV.vtk";
    std::string otputFilename = "/home/or15/Desktop/Marina/BiV_flipped.vtk";
    //inputFilename = "/home/or15/Desktop/Proj/Sander/ControlCMR/VTKs/transformed-5.vtk";
    //otputFilename = "/home/or15/Desktop/Proj/Sander/ControlCMR/Results/testest5.vtk";

    mitk::BaseData::Pointer meshData = mitk::IOUtil::Load(inputFilename).at(0);
    mitk::UnstructuredGrid::Pointer mitkVtkGrid = dynamic_cast<mitk::UnstructuredGrid*>(meshData.GetPointer());
    vtkSmartPointer<vtkUnstructuredGrid> vtkGrid = mitkVtkGrid->GetVtkUnstructuredGrid();
    for (vtkIdType i = 0; i < vtkGrid->GetNumberOfPoints(); i++) {
        double* point = vtkGrid->GetPoint(i);
        point[0] = -point[0];
        point[1] = -point[1];
        //point[2] = -point[2];
        vtkGrid->GetPoints()->SetPoint(i, point);
    }//_for
    vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetInputData(vtkGrid);
    writer->SetFileName(otputFilename.c_str());
    writer->Write();
    /**
    inputFilename = "/home/or15/Desktop/Proj/Sander/ControlCMR/VTKs/transformed-10.vtk";
    otputFilename = "/home/or15/Desktop/Proj/Sander/ControlCMR/Results/testest10.vtk";
    meshData = mitk::IOUtil::Load(inputFilename).at(0);
    mitkVtkGrid = dynamic_cast<mitk::UnstructuredGrid*>(meshData.GetPointer());
    vtkGrid = mitkVtkGrid->GetVtkUnstructuredGrid();
    for (vtkIdType i=0; i<vtkGrid->GetNumberOfPoints(); i++) {
        double* point = vtkGrid->GetPoint(i);
        point[0] = -point[0];
        point[1] = -point[1];
        //point[2] = -point[2];
        vtkGrid->GetPoints()->SetPoint(i, point);
    }//_for
    writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetInputData(vtkGrid);
    writer->SetFileName(otputFilename.c_str());
    writer->Write();

    inputFilename = "/home/or15/Desktop/Proj/Sander/ControlCMR/VTKs/transformed-20.vtk";
    otputFilename = "/home/or15/Desktop/Proj/Sander/ControlCMR/Results/testest20.vtk";
    meshData = mitk::IOUtil::Load(inputFilename).at(0);
    mitkVtkGrid = dynamic_cast<mitk::UnstructuredGrid*>(meshData.GetPointer());
    vtkGrid = mitkVtkGrid->GetVtkUnstructuredGrid();
    for (vtkIdType i=0; i<vtkGrid->GetNumberOfPoints(); i++) {
        double* point = vtkGrid->GetPoint(i);
        point[0] = -point[0];
        point[1] = -point[1];
        //point[2] = -point[2];
        vtkGrid->GetPoints()->SetPoint(i, point);
    }//_for
    writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetInputData(vtkGrid);
    writer->SetFileName(otputFilename.c_str());
    writer->Write();

    inputFilename = "/home/or15/Desktop/Proj/Sander/ControlCMR/VTKs/transformed-25.vtk";
    otputFilename = "/home/or15/Desktop/Proj/Sander/ControlCMR/Results/testest25.vtk";
    meshData = mitk::IOUtil::Load(inputFilename).at(0);
    mitkVtkGrid = dynamic_cast<mitk::UnstructuredGrid*>(meshData.GetPointer());
    vtkGrid = mitkVtkGrid->GetVtkUnstructuredGrid();
    for (vtkIdType i=0; i<vtkGrid->GetNumberOfPoints(); i++) {
        double* point = vtkGrid->GetPoint(i);
        point[0] = -point[0];
        point[1] = -point[1];
        //point[2] = -point[2];
        vtkGrid->GetPoints()->SetPoint(i, point);
    }//_for
    writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetInputData(vtkGrid);
    writer->SetFileName(otputFilename.c_str());
    writer->Write();

    inputFilename = "/home/or15/Desktop/Proj/Sander/ControlCMR/VTKs/transformed-30.vtk";
    otputFilename = "/home/or15/Desktop/Proj/Sander/ControlCMR/Results/testest30.vtk";
    meshData = mitk::IOUtil::Load(inputFilename).at(0);
    mitkVtkGrid = dynamic_cast<mitk::UnstructuredGrid*>(meshData.GetPointer());
    vtkGrid = mitkVtkGrid->GetVtkUnstructuredGrid();
    for (vtkIdType i=0; i<vtkGrid->GetNumberOfPoints(); i++) {
        double* point = vtkGrid->GetPoint(i);
        point[0] = -point[0];
        point[1] = -point[1];
        //point[2] = -point[2];
        vtkGrid->GetPoints()->SetPoint(i, point);
    }//_for
    writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetInputData(vtkGrid);
    writer->SetFileName(otputFilename.c_str());
    writer->Write();

    inputFilename = "/home/or15/Desktop/Proj/Sander/ControlCMR/VTKs/transformed-40.vtk";
    otputFilename = "/home/or15/Desktop/Proj/Sander/ControlCMR/Results/testest40.vtk";
    meshData = mitk::IOUtil::Load(inputFilename).at(0);
    mitkVtkGrid = dynamic_cast<mitk::UnstructuredGrid*>(meshData.GetPointer());
    vtkGrid = mitkVtkGrid->GetVtkUnstructuredGrid();
    for (vtkIdType i=0; i<vtkGrid->GetNumberOfPoints(); i++) {
        double* point = vtkGrid->GetPoint(i);
        point[0] = -point[0];
        point[1] = -point[1];
        //point[2] = -point[2];
        vtkGrid->GetPoints()->SetPoint(i, point);
    }//_for
    writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetInputData(vtkGrid);
    writer->SetFileName(otputFilename.c_str());
    writer->Write();

    inputFilename = "/home/or15/Desktop/Proj/Sander/ControlCMR/VTKs/transformed-45.vtk";
    otputFilename = "/home/or15/Desktop/Proj/Sander/ControlCMR/Results/testest45.vtk";
    meshData = mitk::IOUtil::Load(inputFilename).at(0);
    mitkVtkGrid = dynamic_cast<mitk::UnstructuredGrid*>(meshData.GetPointer());
    vtkGrid = mitkVtkGrid->GetVtkUnstructuredGrid();
    for (vtkIdType i=0; i<vtkGrid->GetNumberOfPoints(); i++) {
        double* point = vtkGrid->GetPoint(i);
        point[0] = -point[0];
        point[1] = -point[1];
        //point[2] = -point[2];
        vtkGrid->GetPoints()->SetPoint(i, point);
    }//_for
    writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetInputData(vtkGrid);
    writer->SetFileName(otputFilename.c_str());
    writer->Write();
    **/
}

void StrainsTests() {

    /**
      * Segmentation overlap test
      **/
    if (false) {
        std::string inputFilename = "/home/or15/Desktop/Proj/JB/OrodNormals/n1/Manuals/transformed-9.vtk";
        std::string outputFilename = "/home/or15/Desktop/Proj/JB/OrodNormals/n1/Manuals/aut-segmentation.mhd";
        // Get all data from the file
        vtkSmartPointer<vtkGenericDataObjectReader> reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
        reader->SetFileName(inputFilename.c_str());
        reader->Update();
        vtkPolyData* output = reader->GetPolyDataOutput();
        for (int i = 0; i < output->GetNumberOfPoints(); i++) {
            double* point = output->GetPoint(i);
            point[0] = -point[0];
            point[1] = -point[1];
            output->GetPoints()->SetPoint(i, point);
        }//_for

        // Write out
        vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();
        double bounds[6];
        output->GetBounds(bounds);

        double voxel = 1;
        double spacing[3]; // desired volume spacing
        spacing[0] = voxel;
        spacing[1] = voxel;
        spacing[2] = voxel;
        whiteImage->SetSpacing(spacing);

        // compute dimensions
        int dim[3];
        for (int i = 0; i < 3; i++) {
            dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]));
        }
        whiteImage->SetDimensions(dim);
        whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

        double origin[3];
        origin[0] = bounds[0] + spacing[0] / 2;
        origin[1] = bounds[2] + spacing[1] / 2;
        origin[2] = bounds[4] + spacing[2] / 2;
        whiteImage->SetOrigin(origin);
        whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

        // fill the image with foreground voxels:
        unsigned char inval = 1;
        unsigned char outval = 0;
        vtkIdType count = whiteImage->GetNumberOfPoints();
        for (vtkIdType i = 0; i < count; ++i) {
            whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
        }

        // polygonal data --> image stencil:
        vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
        pol2stenc->SetInputData(output);
        pol2stenc->SetOutputOrigin(origin);
        pol2stenc->SetOutputSpacing(spacing);
        pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
        pol2stenc->Update();

        // cut the corresponding white image and set the background:
        vtkSmartPointer<vtkImageStencil> imgstenc = vtkSmartPointer<vtkImageStencil>::New();
        imgstenc->SetInputData(whiteImage);
        imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
        imgstenc->ReverseStencilOff();
        imgstenc->SetBackgroundValue(outval);
        imgstenc->Update();

        vtkSmartPointer<vtkMetaImageWriter> writer = vtkSmartPointer<vtkMetaImageWriter>::New();
        writer->SetFileName(outputFilename.c_str());
        writer->SetInputData(imgstenc->GetOutput());
        writer->Write();
    }
    /**
      * Cutting manual test
      **/
    if (false) {

        mitk::Image::Pointer manIm = mitk::IOUtil::Load<mitk::Image>("/home/or15/Desktop/Proj/JB/OrodNormals/n1/Manuals/man-segmentation.nii");
        mitk::Image::Pointer autIm = mitk::IOUtil::Load<mitk::Image>("/home/or15/Desktop/Proj/JB/OrodNormals/n1/Manuals/aut-segmentation.nii");

        mitk::BoundingObject::Pointer cuttingCube;
        QList<mitk::DataNode::Pointer> nodes;// = this->GetDataManagerSelection();
        cuttingCube = dynamic_cast<mitk::BoundingObject*>(nodes.at(0)->GetData());
        mitk::BoundingObjectCutter::Pointer cutter1 = mitk::BoundingObjectCutter::New();
        cutter1->SetBoundingObject(cuttingCube);
        cutter1->SetInput(manIm);
        cutter1->Update();
        mitk::BoundingObjectCutter::Pointer cutter2 = mitk::BoundingObjectCutter::New();
        cutter2->SetBoundingObject(cuttingCube);
        cutter2->SetInput(autIm);
        cutter2->Update();

        mitk::IOUtil::Save(cutter1->GetOutput(), "/home/or15/Desktop/Proj/JB/OrodNormals/n1/Manuals/resultMan.nii");
        mitk::IOUtil::Save(cutter2->GetOutput(), "/home/or15/Desktop/Proj/JB/OrodNormals/n1/Manuals/resultAut.nii");
    }
    /**
      * Tom Tec tests
      **/
    if (false) {
        for (int i = 0; i < 11; i++) {
            mitk::BaseData::Pointer meshData = mitk::IOUtil::Load(
                "/home/or15/Desktop/Proj/JB/TomTec/transformed-" +
                QString::number(i).toStdString() + ".vtk").at(0);
            mitk::UnstructuredGrid::Pointer mitkVtkGrid = dynamic_cast<mitk::UnstructuredGrid*>(meshData.GetPointer());
            vtkSmartPointer<vtkUnstructuredGrid> vtkGrid = mitkVtkGrid->GetVtkUnstructuredGrid();

            vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
            surfaceFilter->SetInputData(vtkGrid);
            surfaceFilter->Update();
            vtkSmartPointer<vtkPolyData> polydata = surfaceFilter->GetOutput();
            mitk::Surface::Pointer converted = mitk::Surface::New();
            converted->SetVtkPolyData(polydata);
            mitk::IOUtil::Save(converted,
                "/home/or15/Desktop/Proj/JB/TomTec/transformed-" +
                QString::number(i).toStdString() + ".vtk");
        }
    }
    /**
      * Papillary tests
      **/
    if (false) {
        QList<mitk::DataNode::Pointer> nodes;// = this->GetDataManagerSelection();
        mitk::PointSet::Pointer s0 = dynamic_cast<mitk::PointSet*>(nodes.at(0)->GetData());
        mitk::PointSet::Pointer s1 = dynamic_cast<mitk::PointSet*>(nodes.at(1)->GetData());
        mitk::PointSet::Pointer s2 = dynamic_cast<mitk::PointSet*>(nodes.at(2)->GetData());
        mitk::PointSet::Pointer s3 = dynamic_cast<mitk::PointSet*>(nodes.at(3)->GetData());
        mitk::PointSet::Pointer s4 = dynamic_cast<mitk::PointSet*>(nodes.at(4)->GetData());
        mitk::PointSet::Pointer s5 = dynamic_cast<mitk::PointSet*>(nodes.at(5)->GetData());
        mitk::PointSet::Pointer s6 = dynamic_cast<mitk::PointSet*>(nodes.at(6)->GetData());
        mitk::PointSet::Pointer s7 = dynamic_cast<mitk::PointSet*>(nodes.at(7)->GetData());
        mitk::PointSet::Pointer s8 = dynamic_cast<mitk::PointSet*>(nodes.at(8)->GetData());
        mitk::PointSet::Pointer s9 = dynamic_cast<mitk::PointSet*>(nodes.at(9)->GetData());
        std::tuple<double, double, double> set0(s0->GetPoint(0).GetElement(0), s0->GetPoint(0).GetElement(1), s0->GetPoint(0).GetElement(2));
        std::tuple<double, double, double> set1(s1->GetPoint(0).GetElement(0), s1->GetPoint(0).GetElement(1), s1->GetPoint(0).GetElement(2));
        std::tuple<double, double, double> set2(s2->GetPoint(0).GetElement(0), s2->GetPoint(0).GetElement(1), s2->GetPoint(0).GetElement(2));
        std::tuple<double, double, double> set3(s3->GetPoint(0).GetElement(0), s3->GetPoint(0).GetElement(1), s3->GetPoint(0).GetElement(2));
        std::tuple<double, double, double> set4(s4->GetPoint(0).GetElement(0), s4->GetPoint(0).GetElement(1), s4->GetPoint(0).GetElement(2));
        std::tuple<double, double, double> set5(s5->GetPoint(0).GetElement(0), s5->GetPoint(0).GetElement(1), s5->GetPoint(0).GetElement(2));
        std::tuple<double, double, double> set6(s6->GetPoint(0).GetElement(0), s6->GetPoint(0).GetElement(1), s6->GetPoint(0).GetElement(2));
        std::tuple<double, double, double> set7(s7->GetPoint(0).GetElement(0), s7->GetPoint(0).GetElement(1), s7->GetPoint(0).GetElement(2));
        std::tuple<double, double, double> set8(s8->GetPoint(0).GetElement(0), s8->GetPoint(0).GetElement(1), s8->GetPoint(0).GetElement(2));
        std::tuple<double, double, double> set9(s9->GetPoint(0).GetElement(0), s9->GetPoint(0).GetElement(1), s9->GetPoint(0).GetElement(2));
        // QString path = "/home/or15/Desktop/Proj/JB/OrodNormals/n1/TrackingPoints";
        // CemrgMeasure* measure = new CemrgMeasure();
        // std::tuple<double,double,double> tra0 = measure->Deconvert(path,0).at(0);
        // std::tuple<double,double,double> tra1 = measure->Deconvert(path,1).at(0);
        // std::tuple<double,double,double> tra2 = measure->Deconvert(path,2).at(0);
        // std::tuple<double,double,double> tra3 = measure->Deconvert(path,3).at(0);
        // std::tuple<double,double,double> tra4 = measure->Deconvert(path,4).at(0);
        // std::tuple<double,double,double> tra5 = measure->Deconvert(path,5).at(0);
        // std::tuple<double,double,double> tra6 = measure->Deconvert(path,6).at(0);
        // std::tuple<double,double,double> tra7 = measure->Deconvert(path,7).at(0);
        // std::tuple<double,double,double> tra8 = measure->Deconvert(path,8).at(0);
        // std::tuple<double,double,double> tra9 = measure->Deconvert(path,9).at(0);
        // qDebug() << measure->CalcDist3D(set0,tra0);
        // qDebug() << measure->CalcDist3D(set1,tra1);
        // qDebug() << measure->CalcDist3D(set2,tra2);
        // qDebug() << measure->CalcDist3D(set3,tra3);
        // qDebug() << measure->CalcDist3D(set4,tra4);
        // qDebug() << measure->CalcDist3D(set5,tra5);
        // qDebug() << measure->CalcDist3D(set6,tra6);
        // qDebug() << measure->CalcDist3D(set7,tra7);
        // qDebug() << measure->CalcDist3D(set8,tra8);
        // qDebug() << measure->CalcDist3D(set9,tra9);
    }
    /**
      * Manual segmentation test
      **/
    if (false) {
        //sanity check
        vtkSmartPointer<vtkCubeSource> mCT = vtkSmartPointer<vtkCubeSource>::New(); mCT->Update();
        vtkSmartPointer<vtkCubeSource> aCT = vtkSmartPointer<vtkCubeSource>::New(); aCT->Update();
        for (int i = 0; i < 4; i++) {
            double* point = mCT->GetOutput()->GetPoint(i);
            point[0] = point[0] * .5;
            point[1] = point[2] * .5;
            point[2] = point[1] * .5;
            aCT->GetOutput()->GetPoints()->SetPoint(i, point);
        }
        mitk::Surface::Pointer mST = mitk::Surface::New(); mST->SetVtkPolyData(mCT->GetOutput());
        mitk::Surface::Pointer aST = mitk::Surface::New(); aST->SetVtkPolyData(aCT->GetOutput());
        mitk::DataNode::Pointer node1 = mitk::DataNode::New(); node1->SetData(mST); //this->GetDataStorage()->Add(node1);
        mitk::DataNode::Pointer node2 = mitk::DataNode::New(); node2->SetData(aST); //this->GetDataStorage()->Add(node2);
        mitk::Surface::Pointer mS = mST;//mitk::IOUtil::Load<mitk::Surface>("/home/or15/Desktop/Proj/JB/OrodNormals/n1/Manuals/man-segmentation.vtk");
        mitk::Surface::Pointer aS = aST;//mitk::IOUtil::Load<mitk::Surface>("/home/or15/Desktop/Proj/JB/OrodNormals/n1/Manuals/aut-segmentation.vtk");
        //return;

        vtkSmartPointer<vtkPolyData> mPD = mS->GetVtkPolyData();
        vtkSmartPointer<vtkPolyData> aPD = aS->GetVtkPolyData();
        vtkSmartPointer<vtkBooleanOperationPolyDataFilter> booleanOperation = vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
        booleanOperation->SetOperationToDifference();
        booleanOperation->SetInputData(0, mPD);
        booleanOperation->SetInputData(1, aPD);
        booleanOperation->Update();
        vtkSmartPointer<vtkPolyData> rPD = booleanOperation->GetOutput();
        mitk::Surface::Pointer rS = mitk::Surface::New();
        rS->SetVtkPolyData(rPD);
        mitk::IOUtil::Save(rS, "/home/or15/Desktop/Proj/JB/OrodNormals/n1/Manuals/res-segmentation.vtk");
        qDebug() << "Filter Done!";

        //double areaM = 0;
        //CemrgStrains* str = new CemrgStrains();
        //for (vtkIdType cellID = 0; cellID < mPD->GetNumberOfCells(); cellID++)
        //    areaM += str->GetCellArea(mPD, cellID);
        //double areaR = 0;
        //for (vtkIdType cellID = 0; cellID < rPD->GetNumberOfCells(); cellID++)
        //    areaR += str->GetCellArea(rPD, cellID);
        //qDebug() << "Difference in %:" << (areaR*100.0)/areaM;
    }
}

void curvesCalculator() {

    for (int i = 0; i < 10; i++) {

        QString directory = "/home/or15/Work/Strain/LR/case0" + QString::number(i);
        mitk::DataNode::Pointer lmNode = mitk::IOUtil::Load<mitk::DataNode>((directory + "/PointSet.mps").toStdString());

        int segRatios[3] = {40, 40, 20};
        std::unique_ptr<CemrgStrains> strain;
        strain = std::unique_ptr<CemrgStrains>(new CemrgStrains(directory, 0));
        strain->ReferenceAHA(lmNode, segRatios, false);
        std::vector<std::vector<double>> plotValueVectorsSQZ;
        std::vector<std::vector<double>> plotValueVectorsCRC;
        std::vector<std::vector<double>> plotValueVectorsLNG;

        for (int j = 0; j < 10; j++) {
            plotValueVectorsSQZ.push_back(strain->CalculateSqzPlot(j));
            plotValueVectorsCRC.push_back(strain->CalculateStrainsPlot(j, lmNode, 3));
            plotValueVectorsLNG.push_back(strain->CalculateStrainsPlot(j, lmNode, 4));
        }

        for (int j = 0; j < 3; j++) {

            QString fileName;
            std::vector<std::vector<double>> plotValueVectors;
            if (j == 0) {
                fileName = "SQZ.csv";
                plotValueVectors = plotValueVectorsSQZ;
            } else if (j == 1) {
                fileName = "CRC.csv";
                plotValueVectors = plotValueVectorsCRC;
            } else {
                fileName = "LNG.csv";
                plotValueVectors = plotValueVectorsLNG;
            }//_if
            std::ofstream file;
            file.open(directory.toStdString() + "/" + fileName.toStdString());

            std::vector<double> values;
            for (int s = 0; s < 16; s++) {
                for (int f = 0; f < 10; f++)
                    values.push_back(plotValueVectors[f][s]);
                //Append the curve to the file
                for (size_t z = 0; z < values.size(); z++) {
                    file << values.at(z);
                    if (z == values.size() - 1) file << endl;
                    else file << ",";
                }
                values.clear();
            }//_for
            file.close();
        }
        qDebug() << "CASE" << i << "done!";
    }
}
