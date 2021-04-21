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
 * http://www.cemrg.co.uk/
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
#include <mitkMeshUtil.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkExtractUnstructuredGrid.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <vtkPolyDataWriter.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkPolygon.h>
#include <vtkDelaunay2D.h>
#include <vtkLookupTable.h>
#include <mitkProperties.h>
#include <mitkLookupTable.h>
#include <mitkStringProperty.h>
#include <dirent.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkLabelShapeKeepNObjectsImageFilter.h>
#include <CemrgCommandLine.h>
#include <vtkDecimatePro.h>
#include <QMessageBox>
#include <itkSimpleContourExtractorImageFilter.h>
#include <vtkLinearExtrusionFilter.h>
#include <vtkSphereSource.h>
#include <itkResampleImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkEllipseSpatialObject.h>
#include <itkSpatialObjectToImageFilter.h>
#include <itkImageFileWriter.h>
#include <mitkBooleanOperation.h>
#include <mitkMorphologicalOperations.h>
#include <vtkDistancePolyDataFilter.h>
#include <vtkConnectivityFilter.h>
#include <vtkThreshold.h>
#include <vtkDoubleArray.h>
#include <mitkManualSegmentationToSurfaceFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkCylinderSource.h>
#include <vtkAppendPolyData.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkMath.h>
#include <vtkLineSource.h>
#include <vtkCleanPolyData.h>
#include <vtkClipPolyData.h>
#include <vtkPlane.h>
#include <vtkTriangle.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkCenterOfMass.h>
#include <vtkIterativeClosestPointTransform.h>
#include <vtkGeometryFilter.h>
#include <vtkLandmarkTransform.h>
#include <itkImageDuplicator.h>
#include <vtkFeatureEdges.h>
#include <vtkCellLocator.h>
#include <mitkImageCast.h>
#include <itkScaleTransform.h>
#include <vtkImageResize.h>
#include <vtkImageChangeInformation.h>
#include <itkBinaryCrossStructuringElement.h>
#include <itkGrayscaleErodeImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryMorphologicalOpeningImageFilter.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include "CemrgAtriaClipper.h"
#include "CemrgScar3D.h"
#include "CemrgMeasure.h"


void autoNIIconvert() {
    
    std::string frames;
    ifstream file0("/home/or15/Desktop/process/frames.txt");
    if (file0.is_open())
        getline(file0,frames);
    file0.close();

    for (int i=1; i<=QString::fromStdString(frames).toInt(); i++) {

        unsigned int slices = 0;
        mitk::Image::Pointer img3D = mitk::Image::New();

        std::string line;
        ifstream file1("/home/or15/Desktop/process/" + QString::number(i).toStdString() + "/path.txt");
        if (file1.is_open()) {
            getline(file1,line);
            while (getline(file1,line)) {
                if (line=="path.txt")
                    break;
                slices++;
            }
        }
        file1.close();
        int ctr = 0;
        bool first = true;
        ifstream file2("/home/or15/Desktop/process/" + QString::number(i).toStdString() + "/path.txt");
        if (file2.is_open()) {
            getline(file2,line);
            while (getline(file2,line)) {
                if (line=="path.txt")
                    break;
                mitk::Image::Pointer image = mitk::IOUtil::LoadImage("/home/or15/Desktop/process/" + QString::number(i).toStdString() + "/" + line);
                if (first) {
                    mitk::ImageDescriptor::Pointer dsc = image->GetImageDescriptor();
                    img3D->Initialize(dsc->GetChannelDescriptor(0).GetPixelType(), *image->GetGeometry(), slices, 1);
                    first = false;
                }
                img3D->SetVolume(mitk::ImageReadAccessor(image).GetData(),0,ctr);
                ctr++;
            }
        }
        file2.close();

        mitk::IOUtil::Save(img3D, "/home/or15/Desktop/process/dcm-" + QString::number(i-1).toStdString() + ".nii");
    }//for
}

void calcSphericity() {

    std::unique_ptr<CemrgMeasure> analys = std::unique_ptr<CemrgMeasure>(new CemrgMeasure());
    std::vector<int> list = {184,141,239,167,245,142,139,209,192,189,225,260,182,187,135,191,255,257,149,168,250,198,180,197,146,131,119,213,177,227,157,107,95,89,101,98,91,100,104,99,105};

    ofstream result;
    QString path = "/home/or15/syncdir/Proj/RZ/AutoAtrialSegmentation/CNN/CalcFibrosis/spherResults.txt";
    result.open(path.toStdString(), std::ios_base::app);

    for (int i=0; i<list.size(); i++) {

        QString aPath = "/home/or15/syncdir/Proj/RZ/AutoAtrialSegmentation/CNN/CalcFibrosis/CEMRA5Observer/CNN/Case" + QString::number(list.at(i)) + "/scar.vtk";
        QString mPath = "/home/or15/syncdir/Proj/RZ/AutoAtrialSegmentation/CNN/CalcFibrosis/CEMRA5Observer/GTR/Case" + QString::number(list.at(i)) + "/scar.vtk";

        vtkSmartPointer<vtkPolyData> aMesh = mitk::IOUtil::LoadSurface(aPath.toStdString())->GetVtkPolyData();
        vtkSmartPointer<vtkPolyData> mMesh = mitk::IOUtil::LoadSurface(mPath.toStdString())->GetVtkPolyData();

        double aSpher = analys->GetSphericity(aMesh);
        double mSpher = analys->GetSphericity(mMesh);

        result << aSpher << " " << mSpher << "\n";
    }
    result.close();
}

void convertVTKtoMPS() {

    std::unique_ptr<CemrgMeasure> rr(new CemrgMeasure());
    std::vector <std::tuple<double, double, double>> points;
    QString directory = "/home/or15/Desktop/Marina/Points/Dataset";

    for (int caseNo=17; caseNo<18; caseNo++) {
        if (caseNo==4) continue;
        QString dirc = directory + QString::number(caseNo) + "/PMlower";
        for (int frame=0; frame<10; frame++) {
            points = rr->Deconvert(dirc,frame);
            mitk::PointSet::Pointer set = mitk::PointSet::New();
            for (unsigned int i=0; i<points.size(); i++) {
                mitk::Point3D p;
                p.SetElement(0, std::get<0>(points.at(i)));
                p.SetElement(1, std::get<1>(points.at(i)));
                p.SetElement(2, std::get<2>(points.at(i)));
                set->InsertPoint(i,p);
            }
            mitk::IOUtil::Save(set, (dirc + "/PMlower-" + QString::number(frame) + ".mps").toStdString());
            points.clear();
        }
        dirc = directory + QString::number(caseNo) + "/PMupper";
        for (int frame=0; frame<10; frame++) {
            points = rr->Deconvert(dirc,frame);
            mitk::PointSet::Pointer set = mitk::PointSet::New();
            for (unsigned int i=0; i<points.size(); i++) {
                mitk::Point3D p;
                p.SetElement(0, std::get<0>(points.at(i)));
                p.SetElement(1, std::get<1>(points.at(i)));
                p.SetElement(2, std::get<2>(points.at(i)));
                set->InsertPoint(i,p);
            }
            mitk::IOUtil::Save(set, (dirc + "/PMupper-" + QString::number(frame) + ".mps").toStdString());
            points.clear();
        }
        dirc = directory + QString::number(caseNo) + "/MAA";
        for (int frame=0; frame<10; frame++) {
            points = rr->Deconvert(dirc,frame);
            mitk::PointSet::Pointer set = mitk::PointSet::New();
            for (unsigned int i=0; i<points.size(); i++) {
                mitk::Point3D p;
                p.SetElement(0, std::get<0>(points.at(i)));
                p.SetElement(1, std::get<1>(points.at(i)));
                p.SetElement(2, std::get<2>(points.at(i)));
                set->InsertPoint(i,p);
            }
            mitk::IOUtil::Save(set, (dirc + "/MAA-" + QString::number(frame) + ".mps").toStdString());
            points.clear();
        }
    }
}

void generateAtlas() {

    for (int v=0; v<6; v++) {
        for (int h=0; h<5; h++) {

            //Loop cases
            for (unsigned int i=1; i<=36; i++) {
                //Skip cases
                if (i==4||i==11||i==13||i==18||i==21||i==22||i==23)
                    continue;
                std::string pDirc = "/home/or15/Desktop/Proj/SN/LVMaps/CylOptimised/FlattenedCases";
                std::string path1 = pDirc + "/Case" + QString::number(i).toStdString() + "-flat/Volume" + QString::number(v).toStdString() + "/" + QString::number(h).toStdString() + ".vtk";
                std::string path2 = pDirc + "/transformed-" + QString::number(i).toStdString() + ".vtk";
                std::string comnd = "cp " + path1 + " " + path2;
                system(comnd.c_str());
            }

            std::string mapPath = "/home/or15/Desktop/Proj/SN/LVMaps/CylOptimised/Results/EmptyFlatMap.vtk";
            vtkSmartPointer<vtkPolyData> map = mitk::IOUtil::LoadSurface(mapPath)->GetVtkPolyData();
            vtkSmartPointer<vtkFloatArray> intersectionScore = vtkSmartPointer<vtkFloatArray>::New();
            for (vtkIdType ptID = 0; ptID<map->GetNumberOfPoints(); ptID++)
                intersectionScore->InsertNextTuple1(0);

            //Loop cases
            for (unsigned int i=1; i<=36; i++) {

                //Skip cases
                if (i==4||i==11||i==13||i==18||i==21||i==22||i==23)
                    continue;

                std::string casPath = "/home/or15/Desktop/Proj/SN/LVMaps/CylOptimised/FlattenedCases/transformed-" + QString::number(i).toStdString() + ".vtk";
                vtkSmartPointer<vtkPolyData> cas = mitk::IOUtil::LoadSurface(casPath)->GetVtkPolyData();
                vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
                pointLocator->SetDataSet(cas);
                pointLocator->BuildLocator();

                for (vtkIdType ptID = 0; ptID<map->GetNumberOfPoints(); ptID++) {
                    double* point = map->GetPoint(ptID);
                    vtkIdType ID = pointLocator->FindClosestPoint(point);
                    double score = cas->GetPointData()->GetScalars("unknown")->GetTuple1(ID);
                    double total_score = intersectionScore->GetTuple1(ptID) + score;
                    intersectionScore->InsertTuple1(ptID, total_score);
                }
            }
            for (vtkIdType ptID = 0; ptID<map->GetNumberOfPoints(); ptID++) {
                double total_score = intersectionScore->GetTuple1(ptID);
                intersectionScore->InsertTuple1(ptID, total_score/29.0);
            }

            map->GetPointData()->SetScalars(intersectionScore);
            mitk::Surface::Pointer output = mitk::Surface::New();
            output->SetVtkPolyData(map);
            std::string outPath = "/home/or15/Desktop/Proj/SN/LVMaps/CylOptimised/Results/Map-" + QString::number(v).toStdString() + "-" + QString::number(h+1).toStdString() + ".vtk";
            mitk::IOUtil::Save(output, outPath);

            //Clean up
            std::string clnDirec = "/home/or15/Desktop/Proj/SN/LVMaps/CylOptimised/FlattenedCases";
            std::string clnPaths = clnDirec + "/transformed-*.vtk";
            std::string clnComnd = "rm " + clnPaths;
            system(clnComnd.c_str());
        }
    }
}

void singleCurvesCalculator() {

    int segRatios[3] = {40, 40, 20};

    //Loop cases
    for (unsigned int c=1; c<=36; c++) {

        //Skip cases
        if (c==4||c==11||c==13||c==18||c==21||c==22||c==23)
            continue;

        //Setup the surface
        for (unsigned int v=0; v<6; v++) {
            for (int i=0; i<5; i++) {

                QString directory = "/home/or15/Desktop/Proj/SN/LVMaps/CylOptimised/FlattenedCases/Case" + QString::number(c) + "-flat/Volume" + QString::number(v);
                mitk::DataNode::Pointer lmNode = mitk::IOUtil::LoadDataNode((directory + "/PointSet.mps").toStdString());

                std::unique_ptr<CemrgStrains> strain;
                strain = std::unique_ptr<CemrgStrains>(new CemrgStrains(directory, i));
                strain->ReferenceAHA(lmNode, segRatios);
                mitk::Surface::Pointer surface = strain->FlattenedAHA();
                surface->GetVtkPolyData()->GetCellData()->SetScalars(strain->GetFlatSurfScalars());
                mitk::IOUtil::Save(surface, (directory + "/" + QString::number(i) + ".vtk").toStdString());
            }
        }
        qDebug() << "CASE" << c << "done!";
    }
}

bool cylTest(const double* pt1, const double* pt2, const double* testpt, double length_sq, double radius_sq) {

    double dx, dy, dz;
    double pdx, pdy, pdz;
    double dot, dsq;
    dx = pt2[0] - pt1[0];
    dy = pt2[1] - pt1[1];
    dz = pt2[2] - pt1[2];
    pdx = testpt[0] - pt1[0];
    pdy = testpt[1] - pt1[1];
    pdz = testpt[2] - pt1[2];
    dot = pdx * dx + pdy * dy + pdz * dz;
    if (dot < 0.0 || dot > length_sq) {
        return false;
    } else {
        dsq = (pdx*pdx + pdy*pdy + pdz*pdz) - dot*dot/length_sq;
        if (dsq > radius_sq) return false;
        else return true;
    }
}

void generateLocationMap() {

    //Create a cylinder
    const double radiusAbbott = 5.99/2.0;
    const double heightAbbott = 42.0;
    const double volumeAbbott = 1000;
    const double radiusMedtronic = 6.7/2.0;
    const double heightMedtronic = 25.9;
    const double volumeMedtronic = 800;

    for (unsigned int v=0; v<6; v++) {
        for (int h=1; h<6; h++) {

            //Optimised measures
            std::vector<double> volumes = {800, 720, 648, 583.2, 524.88, 472.39};
            double height = h*5.0;
            double radius = sqrt(volumes.at(v) / (M_PI * height));

            //Loop cases
            for (unsigned int i=1; i<37; i++) {

                //Skip cases
                if (i==4||i==11||i==13||i==18||i==21||i==22||i==23)
                    continue;

                //Prepare heatmap
                QString dirc = "/home/or15/Desktop/Proj/SN/LVMaps/CylOptimised/Case";
                QString path0 = dirc + QString::number(i) + "/transformed-0_decimated.vtk";
                mitk::Surface::Pointer refLV = mitk::IOUtil::LoadSurface(path0.toStdString());
                vtkSmartPointer<vtkIntArray> intersectionScore = vtkSmartPointer<vtkIntArray>::New();
                for (vtkIdType ptID = 0; ptID<refLV->GetVtkPolyData()->GetNumberOfPoints(); ptID++)
                    intersectionScore->InsertNextTuple1(0);

                //Define MV plane
                double centre[3];
                vtkSmartPointer<vtkCenterOfMass> centerOfMassFilter = vtkSmartPointer<vtkCenterOfMass>::New();
                centerOfMassFilter->SetInputData(refLV->GetVtkPolyData());
                centerOfMassFilter->SetUseScalarsAsWeights(false);
                centerOfMassFilter->Update();
                centerOfMassFilter->GetCenter(centre);
                QString path3 = dirc + QString::number(i) + "/PointSet.mps";
                mitk::PointSet::Pointer pts = mitk::IOUtil::LoadPointSet(path3.toStdString());
                double nrm[3] = {0}; double pt1[3] = {0}; double pt2[3] = {0}; double pt3[3] = {0};
                pt1[0] = -pts->GetPoint(1).GetElement(0); pt1[1] = -pts->GetPoint(1).GetElement(1); pt1[2] = pts->GetPoint(1).GetElement(2);
                pt2[0] = -pts->GetPoint(2).GetElement(0); pt2[1] = -pts->GetPoint(2).GetElement(1); pt2[2] = pts->GetPoint(2).GetElement(2);
                pt3[0] = -pts->GetPoint(3).GetElement(0); pt3[1] = -pts->GetPoint(3).GetElement(1); pt3[2] = pts->GetPoint(3).GetElement(2);
                vtkTriangle::ComputeNormal(pt1, pt2, pt3, nrm);
                bool side = false;
                double product = 0;
                double mvPt[3] = {0};
                mvPt[0] = pt1[0] - centre[0];
                mvPt[1] = pt1[1] - centre[1];
                mvPt[2] = pt1[2] - centre[2];
                for (int i = 0; i < 3; i++)
                    product += (mvPt[i])*(nrm[i]);
                if (product>0) side = true;

                //Loop frames
                for (int j=0; j<(i<24?10:20); (i<24?j=j+1:j=j+2)) {

                    //Decimate LV
                    QString path0 = dirc + QString::number(i) + "/transformed-0.vtk";
                    QString pthD0 = dirc + QString::number(i) + "/transformed-0_decimated.vtk";
                    QString mirtk = "/home/or15/Install/MIRTK/MIRTK-build/lib/mirtk-remesh";
                    QString commd = mirtk + " " + path0 + " " + pthD0 + " -edgelength 7 -max-edgelength 10 -min-edgelength 5 -melt-nodes no";
                    if (j==0) std::system(commd.toStdString().c_str());

                    QString pthDN = dirc + QString::number(i) + "/transformed-" + QString::number(j) + "_decimated.vtk";
                    mirtk = "/home/or15/Install/MIRTK/MIRTK-build/lib/mirtk-transform-points";
                    QString doDir = "/home/or15/Desktop/Proj/SN/LVMaps/Aux/Marina/Structures/Points/Dataset" + (i<33 ? QString::number(i) : "EBR" + QString::number(i-32));
                    QString dofil = doDir + (i<33 ? "/C" + QString::number(i) : "/EBR" + QString::number(i-32)) + ".dof";
                    commd = mirtk + " " + pthD0 + " " + pthDN + " -dofin " + dofil + " -ascii -St " + QString::number(j*10);
                    std::system(commd.toStdString().c_str());

                    //Load LV
                    vtkSmartPointer<vtkPolyData> LV = mitk::IOUtil::LoadSurface(pthDN.toStdString())->GetVtkPolyData();

                    //Load Cones
                    QString path2 = dirc + QString::number(i) + "/cone_PM1-" + QString::number(j) + ".vtk";
                    vtkSmartPointer<vtkUnstructuredGridReader> C1G = vtkSmartPointer<vtkUnstructuredGridReader>::New();
                    C1G->SetFileName(path2.toStdString().c_str());
                    C1G->Update();
                    vtkSmartPointer<vtkUnstructuredGrid> C1 = C1G->GetOutput();
                    QString path3 = dirc + QString::number(i) + "/cone_PM2-" + QString::number(j) + ".vtk";
                    vtkSmartPointer<vtkUnstructuredGridReader> C2G = vtkSmartPointer<vtkUnstructuredGridReader>::New();
                    C2G->SetFileName(path3.toStdString().c_str());
                    C2G->Update();
                    vtkSmartPointer<vtkUnstructuredGrid> C2 = C2G->GetOutput();
                    //Load Cylinders
                    QString path4 = dirc + QString::number(i) + "/cylinder_PM1-" + QString::number(j) + ".vtk";
                    vtkSmartPointer<vtkUnstructuredGridReader> P1G = vtkSmartPointer<vtkUnstructuredGridReader>::New();
                    P1G->SetFileName(path4.toStdString().c_str());
                    P1G->Update();
                    vtkSmartPointer<vtkUnstructuredGrid> P1 = P1G->GetOutput();
                    QString path5 = dirc + QString::number(i) + "/cylinder_PM2-" + QString::number(j) + ".vtk";
                    vtkSmartPointer<vtkUnstructuredGridReader> P2G = vtkSmartPointer<vtkUnstructuredGridReader>::New();
                    P2G->SetFileName(path5.toStdString().c_str());
                    P2G->Update();
                    vtkSmartPointer<vtkUnstructuredGrid> P2 = P2G->GetOutput();

                    //Generate normals
                    vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
                    normalGenerator->SetInputData(LV);
                    normalGenerator->ComputePointNormalsOn();
                    normalGenerator->ComputeCellNormalsOff();
                    normalGenerator->AutoOrientNormalsOn();
                    normalGenerator->FlipNormalsOn();
                    normalGenerator->Update();
                    LV = normalGenerator->GetOutput();

                    //Loop the locations
                    for (vtkIdType ptID = 0; ptID<LV->GetNumberOfPoints(); ptID++) {

                        //Ignore checked locations
                        if (fabs(intersectionScore->GetTuple1(ptID)-1)<1E-5) continue;

                        //Retrieve starting location
                        double strPos[3] = {0};
                        strPos[0] = LV->GetPoint(ptID)[0];
                        strPos[1] = LV->GetPoint(ptID)[1];
                        strPos[2] = LV->GetPoint(ptID)[2];

                        //Ignore points above MV
                        double product = 0;
                        double mvPt[3] = {0};
                        mvPt[0] = pt1[0] - strPos[0];
                        mvPt[1] = pt1[1] - strPos[1];
                        mvPt[2] = pt1[2] - strPos[2];
                        for (int i = 0; i < 3; i++)
                            product += (mvPt[i])*(nrm[i]);
                        if (side) {
                            if (product<0) continue;
                        } else {
                            if (product>0) continue;
                        }//_if

                        //Retrieve point normal and end position
                        double normal[3] = {0};
                        normal[0] = LV->GetPointData()->GetArray("Normals")->GetTuple(ptID)[0];
                        normal[1] = LV->GetPointData()->GetArray("Normals")->GetTuple(ptID)[1];
                        normal[2] = LV->GetPointData()->GetArray("Normals")->GetTuple(ptID)[2];
                        double endPos[3] = {0};
                        endPos[0] = strPos[0] + height * normal[0];
                        endPos[1] = strPos[1] + height * normal[1];
                        endPos[2] = strPos[2] + height * normal[2];

                        //Check cone1
                        for (vtkIdType ID = 0; ID<C1->GetNumberOfPoints(); ID++) {
                            //Retrieve test location
                            double tstPoint[3] = {0.0};
                            tstPoint[0] = -C1->GetPoint(ID)[0];
                            tstPoint[1] = -C1->GetPoint(ID)[1];
                            tstPoint[2] =  C1->GetPoint(ID)[2];
                            //Collision in filtered points
                            bool result = cylTest(strPos, endPos, tstPoint, pow(height,2), pow(radius,2));
                            if (result) {
                                intersectionScore->InsertTuple1(ptID, 1);
                                goto nextLocation;
                            }
                        }

                        //Check cone2
                        for (vtkIdType ID = 0; ID<C2->GetNumberOfPoints(); ID++) {
                            //Retrieve test location
                            double tstPoint[3] = {0.0};
                            tstPoint[0] = -C2->GetPoint(ID)[0];
                            tstPoint[1] = -C2->GetPoint(ID)[1];
                            tstPoint[2] =  C2->GetPoint(ID)[2];
                            //Collision in filtered points
                            bool result = cylTest(strPos, endPos, tstPoint, pow(height,2), pow(radius,2));
                            if (result) {
                                intersectionScore->InsertTuple1(ptID, 1);
                                goto nextLocation;
                            }
                        }

                        //Check papillary1
                        for (vtkIdType ID = 0; ID<P1->GetNumberOfPoints(); ID++) {
                            //Retrieve test location
                            double tstPoint[3] = {0.0};
                            tstPoint[0] = -P1->GetPoint(ID)[0];
                            tstPoint[1] = -P1->GetPoint(ID)[1];
                            tstPoint[2] =  P1->GetPoint(ID)[2];
                            //Collision in filtered points
                            bool result = cylTest(strPos, endPos, tstPoint, pow(height,2), pow(radius,2));
                            if (result) {
                                intersectionScore->InsertTuple1(ptID, 1);
                                goto nextLocation;
                            }
                        }

                        //Check papillary2
                        for (vtkIdType ID = 0; ID<P2->GetNumberOfPoints(); ID++) {
                            //Retrieve test location
                            double tstPoint[3] = {0.0};
                            tstPoint[0] = -P2->GetPoint(ID)[0];
                            tstPoint[1] = -P2->GetPoint(ID)[1];
                            tstPoint[2] =  P2->GetPoint(ID)[2];
                            //Collision in filtered points
                            bool result = cylTest(strPos, endPos, tstPoint, pow(height,2), pow(radius,2));
                            if (result) {
                                intersectionScore->InsertTuple1(ptID, 1);
                                goto nextLocation;
                            }
                        }

                        //Check endo of LV                        
                        for (vtkIdType ID = 0; ID<LV->GetNumberOfPoints(); ID++) {
                            //Retrieve test location
                            double tstPoint[3] = {0.0};
                            tstPoint[0] = LV->GetPoint(ID)[0];
                            tstPoint[1] = LV->GetPoint(ID)[1];
                            tstPoint[2] = LV->GetPoint(ID)[2];
                            //Ignore points above MV
                            double product = 0;
                            double mvPt[3] = {0};
                            mvPt[0] = pt1[0] - tstPoint[0];
                            mvPt[1] = pt1[1] - tstPoint[1];
                            mvPt[2] = pt1[2] - tstPoint[2];
                            for (int i = 0; i < 3; i++)
                                product += (mvPt[i])*(nrm[i]);
                            if (side) {
                                if (product<0) continue;
                            } else {
                                if (product>0) continue;
                            }
                            //Ignore too close points
                            double distance = sqrt(vtkMath::Distance2BetweenPoints(strPos, tstPoint));
                            if (distance<(radius*2*1.3)) continue;
                            //Collision in filtered points
                            bool result = cylTest(strPos, endPos, tstPoint, pow(height,2), pow(radius,2));
                            if (result) {
                                intersectionScore->InsertTuple1(ptID, 1);
                                goto nextLocation;
                            }
                        }
nextLocation:
                        continue;
                    }//_loop_locations
                }//_loop_frames

                //Store results
                refLV->GetVtkPolyData()->GetPointData()->AddArray(intersectionScore);
                vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
                plane->SetNormal(nrm);
                plane->SetOrigin(pt1);
                vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
                clipper->SetClipFunction(plane);
                clipper->SetInputData(refLV->GetVtkPolyData());
                if (side) clipper->InsideOutOn();
                else clipper->InsideOutOff();
                clipper->Update();
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
                refLV->SetVtkPolyData(cleaner->GetOutput());
                mitk::IOUtil::Save(refLV, (dirc + QString::number(i) + "/CylMap_" + QString::number(v) + "_" + QString::number(h) + ".vtk").toStdString());
                qDebug() << "Case" << i << v << h;
            }
        }//_heights
    }//_volumes
}

void newCTveinClipperBasedonDistance() {

    //Test cases
    QString direct = "/home/or15/Desktop/process/";
    QString segBPPath = direct + "PVeinsCroppedImage.nii";
    QString segLAPath = direct + "leftAtrialWall900201.nrrd";
    mitk::Image::Pointer segBP = mitk::IOUtil::LoadImage(segBPPath.toStdString());
    mitk::Image::Pointer segLA = mitk::IOUtil::LoadImage(segLAPath.toStdString());

    //Conversion to ITK
    typedef itk::Image<float, 3> ImageTypeFloat;
    typedef itk::Image<unsigned char, 3> ImageTypeChar;
    ImageTypeChar::Pointer segBPinITK = ImageTypeChar::New();
    ImageTypeChar::Pointer segLAinITK = ImageTypeChar::New();
    CastToItkImage(segBP, segBPinITK);
    CastToItkImage(segLA, segLAinITK);

    //Distance map
    typedef itk::SignedMaurerDistanceMapImageFilter<ImageTypeChar, ImageTypeFloat> DistanceFilterType;
    DistanceFilterType::Pointer distanceFilter = DistanceFilterType::New();
    distanceFilter->SetInput(segBPinITK);
    distanceFilter->SetSquaredDistance(false);
    distanceFilter->SetInsideIsPositive(false);
    distanceFilter->SetUseImageSpacing(true);
    try {
        distanceFilter->Update();
    } catch ( ... ) {
        return;
    }

    //Threshold
    typedef itk::ThresholdImageFilter<ImageTypeFloat> ThreshFilterType;
    ThreshFilterType::Pointer threshFilter1 = ThreshFilterType::New();
    threshFilter1->SetInput(distanceFilter->GetOutput());
    threshFilter1->ThresholdOutside(0,5);
    threshFilter1->SetOutsideValue(0);
    threshFilter1->Update();
    ThreshFilterType::Pointer threshFilter2 = ThreshFilterType::New();
    threshFilter2->SetInput(threshFilter1->GetOutput());
    threshFilter2->ThresholdAbove(0);
    threshFilter2->SetOutsideValue(1);
    threshFilter2->Update();

    //Clip
    typedef itk::ImageRegionIteratorWithIndex<ImageTypeChar> ItType1;
    typedef itk::ImageRegionIteratorWithIndex<ImageTypeFloat> ItType2;
    ItType1 it1(segLAinITK, segLAinITK->GetRequestedRegion());
    ItType2 it2(threshFilter2->GetOutput(), threshFilter2->GetOutput()->GetRequestedRegion());
    it1.GoToBegin();
    for (it2.GoToBegin(); !it2.IsAtEnd(); ++it2) {
        if (static_cast<int>(it2.Get()) == 0)
                it1.Set(0);
        ++it1;
    }//_for

    //Output
    QString distPath = direct + "output.nii";
    mitk::IOUtil::Save(mitk::ImportItkImage(segLAinITK), distPath.toStdString());
}

void veinClipperBasedonDistance() {

    //Test cases
    QString direct = "/home/or15/Desktop/Test/";
    std::vector<int> testList = {72,77,25,37,81,96,46,99,39,65,58,12,97,88,70,87,36,21,83,9,103,106,67,64,47,44};
    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());

    for (int i=0; i<testList.size(); i++) {
        //Find veins
        QString index = QString::number(testList.at(i));
        QString segLAPath = direct + "LA_Predictions/thrsh_image_" + index + ".nii.gz";
        QString segVNPath = direct + "VN_Predictions/thrsh_image_" + index + ".nii.gz";
        mitk::Image::Pointer segLA = mitk::IOUtil::LoadImage(segLAPath.toStdString());
        mitk::Image::Pointer segVN = mitk::IOUtil::LoadImage(segVNPath.toStdString());

        bool smooth(true);
        bool applyMedian(true);
        unsigned int medianKernelSize(3);
        float gaussianSD(1.0);
        QString meshLA = direct + "/LA_" + QString::number(i) + ".vtk";
        mitk::ManualSegmentationToSurfaceFilter::Pointer surfaceFilter1 = mitk::ManualSegmentationToSurfaceFilter::New();
        surfaceFilter1->SetInput(segLA);
        surfaceFilter1->SetThreshold(0.5);
        surfaceFilter1->SetUseGaussianImageSmooth(smooth);
        surfaceFilter1->SetSmooth(smooth);
        surfaceFilter1->InterpolationOn();
        surfaceFilter1->SetGaussianStandardDeviation(gaussianSD);
        surfaceFilter1->SetMedianFilter3D(applyMedian);
        surfaceFilter1->SetMedianKernelSize(medianKernelSize, medianKernelSize, medianKernelSize);
        surfaceFilter1->SetDecimate(mitk::ImageToSurfaceFilter::NoDecimation);
        surfaceFilter1->UpdateLargestPossibleRegion();
        mitk::Surface::Pointer temp1 = surfaceFilter1->GetOutput();
        vtkSmartPointer<vtkPolyData> polyData1 = temp1->GetVtkPolyData();
        polyData1->SetVerts(nullptr);
        polyData1->SetLines(nullptr);
        vtkSmartPointer<vtkPolyDataNormals> normalsGen1 = vtkSmartPointer<vtkPolyDataNormals>::New();
        normalsGen1->AutoOrientNormalsOn();
        normalsGen1->FlipNormalsOff();
        normalsGen1->SetInputData(polyData1);
        normalsGen1->Update();
        temp1->SetVtkPolyData(normalsGen1->GetOutput());
        mitk::IOUtil::SaveSurface(temp1, meshLA.toStdString());
        //cmd->ExecuteSurf(direct, segLAPath, 1, 0.5, 0, 0);
        //std::rename((direct + "/segmentation.vtk").toStdString().c_str(), meshLA.toStdString().c_str());
        QString meshVN = direct + "/VN_" + QString::number(i) + ".vtk";
        mitk::ManualSegmentationToSurfaceFilter::Pointer surfaceFilter2 = mitk::ManualSegmentationToSurfaceFilter::New();
        surfaceFilter2->SetInput(segVN);
        surfaceFilter2->SetThreshold(0.5);
        surfaceFilter2->SetUseGaussianImageSmooth(smooth);
        surfaceFilter2->SetSmooth(smooth);
        surfaceFilter2->InterpolationOn();
        surfaceFilter2->SetGaussianStandardDeviation(gaussianSD);
        surfaceFilter2->SetMedianFilter3D(applyMedian);
        surfaceFilter2->SetMedianKernelSize(medianKernelSize, medianKernelSize, medianKernelSize);
        surfaceFilter2->SetDecimate(mitk::ImageToSurfaceFilter::NoDecimation);
        surfaceFilter2->UpdateLargestPossibleRegion();
        mitk::Surface::Pointer temp2 = surfaceFilter2->GetOutput();
        vtkSmartPointer<vtkPolyData> polyData2 = temp2->GetVtkPolyData();
        polyData2->SetVerts(nullptr);
        polyData2->SetLines(nullptr);
        vtkSmartPointer<vtkPolyDataNormals> normalsGen2 = vtkSmartPointer<vtkPolyDataNormals>::New();
        normalsGen2->AutoOrientNormalsOn();
        normalsGen2->FlipNormalsOff();
        normalsGen2->SetInputData(polyData2);
        normalsGen2->Update();
        temp2->SetVtkPolyData(normalsGen2->GetOutput());
        mitk::IOUtil::SaveSurface(temp2, meshVN.toStdString());
        //cmd->ExecuteSurf(direct, segVNPath, 1, 0.5, 0, 0);
        //std::rename((direct + "/segmentation.vtk").toStdString().c_str(), meshVN.toStdString().c_str());

        vtkSmartPointer<vtkPolyData> model1 = mitk::IOUtil::LoadSurface(meshLA.toStdString())->GetVtkPolyData();
        vtkSmartPointer<vtkPolyData> model2 = mitk::IOUtil::LoadSurface(meshVN.toStdString())->GetVtkPolyData();
        vtkSmartPointer<vtkDistancePolyDataFilter> distanceFilter = vtkSmartPointer<vtkDistancePolyDataFilter>::New();
        distanceFilter->SetInputData(0, model1);
        distanceFilter->SetInputData(1, model2);
        distanceFilter->Update();

        double valuesRange[2];
        vtkSmartPointer<vtkDoubleArray> distances = vtkSmartPointer<vtkDoubleArray>::New();
        distances = vtkDoubleArray::SafeDownCast(distanceFilter->GetOutput()->GetPointData()->GetArray("Distance"));
        distances->GetValueRange(valuesRange);
        double sum = 0.0;
        double sumDeviation = 0.0;
        double mean = 0.0; double stdv = 0.0;
        for (unsigned int j=0; j<distances->GetNumberOfTuples(); j++)
            sum += distances->GetTuple1(j);
        mean = sum / distances->GetNumberOfTuples();
        for (unsigned int j=0; j<distances->GetNumberOfTuples(); j++)
            sumDeviation += (distances->GetTuple1(j)-mean) * (distances->GetTuple1(j)-mean);
        stdv = std::sqrt(sumDeviation/distances->GetNumberOfTuples());

        vtkSmartPointer<vtkThreshold> threshold = vtkSmartPointer<vtkThreshold>::New();
        threshold->SetInputConnection(distanceFilter->GetOutputPort());
        threshold->ThresholdByUpper(mean+1.5*stdv);
        threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, vtkDataSetAttributes::SCALARS);
        threshold->Update();
        vtkSmartPointer<vtkConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkConnectivityFilter>::New();
        connectivityFilter->SetInputConnection(threshold->GetOutputPort());
        connectivityFilter->SetExtractionModeToAllRegions();
        connectivityFilter->ColorRegionsOn();
        connectivityFilter->Update();
        vtkSmartPointer<vtkDataSetSurfaceFilter> filter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
        filter->SetInputData(connectivityFilter->GetOutput());
        filter->Update();
        mitk::Surface::Pointer surf = mitk::Surface::New();
        surf->SetVtkPolyData(filter->GetOutput());
        mitk::IOUtil::Save(surf, "/home/or15/Desktop/Test/test_" + QString::number(i).toStdString() + ".vtk");
    }
}

void autoScarCalculationsBasedonShape() {

    QString direct = "/home/cc14/samples_UT_ascii";
    QString lgePath = direct + "/dcm-LGE.nii";
    typedef itk::Image<short, 3> ImageTypeCHAR;
    typedef itk::Image<short, 3> ImageTypeSHRT;

    //Scar projection
    int minStep = -1;
    int maxStep = 3;
    int methodType = 2;
    std::unique_ptr<CemrgScar3D> scar(new CemrgScar3D());
    scar->SetMinStep(minStep);
    scar->SetMaxStep(maxStep);
    scar->SetMethodType(methodType);
    ImageTypeCHAR::Pointer segITK = ImageTypeCHAR::New();
    CastToItkImage(mitk::IOUtil::LoadImage((direct + "/PVeinsCroppedImage.nii").toStdString()), segITK);
    ImageTypeSHRT::Pointer lgeITK = ImageTypeSHRT::New();
    CastToItkImage(mitk::IOUtil::LoadImage(lgePath.toStdString()), lgeITK);
    itk::ResampleImageFilter<ImageTypeCHAR, ImageTypeCHAR>::Pointer resampleFilter;
    resampleFilter = itk::ResampleImageFilter<ImageTypeCHAR, ImageTypeCHAR>::New();
    resampleFilter->SetInput(segITK);
    resampleFilter->SetReferenceImage(lgeITK);
    resampleFilter->SetUseReferenceImage(true);
    resampleFilter->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<ImageTypeCHAR>::New());
    resampleFilter->SetDefaultPixelValue(0);
    resampleFilter->UpdateLargestPossibleRegion();
    segITK = resampleFilter->GetOutput();
    mitk::IOUtil::Save(mitk::ImportItkImage(segITK), (direct + "/PVeinsCroppedImage.nii").toStdString());
    scar->SetScarSegImage(mitk::ImportItkImage(segITK));

    //Thresholding
    int vxls = 3;
    int threshType = 1;
    double value1 = 0.97;
    double value2 = 1.61;
    double value3 = 3.3;
    typedef itk::Image<float, 3> ImageType;
    typedef itk::BinaryBallStructuringElement<ImageTypeCHAR::PixelType, 3> BallType;
    typedef itk::GrayscaleErodeImageFilter<ImageTypeCHAR, ImageType, BallType> ErosionFilterType;
    BallType binaryBall;
    binaryBall.SetRadius(vxls);
    binaryBall.CreateStructuringElement();
    ErosionFilterType::Pointer erosionFilter = ErosionFilterType::New();
    erosionFilter->SetInput(segITK);
    erosionFilter->SetKernel(binaryBall);
    erosionFilter->UpdateLargestPossibleRegion();
    mitk::Image::Pointer roiImage = mitk::Image::New();
    roiImage = mitk::ImportItkImage(erosionFilter->GetOutput())->Clone();
    ImageType::Pointer lgeFloat = ImageType::New();
    CastToItkImage(mitk::IOUtil::LoadImage(lgePath.toStdString()), lgeFloat);
    double mean = 0.0, stdv = 0.0;
    scar->CalculateMeanStd(mitk::ImportItkImage(lgeFloat), roiImage, mean, stdv);
    double thresh1 = mean*value1;
    double thresh2 = mean*value2;
    double thresh3 = mean+value3*stdv;
    QString prodPath = direct + mitk::IOUtil::GetDirectorySeparator();

    for (int i=1; i<=16; i++) {

        //mitk::Surface::Pointer scarShell = scar->Scar3D(direct.toStdString(), mitk::ImportItkImage(lgeITK),i);
        //mitk::IOUtil::Save(scarShell, (direct + "/scar" + QString::number(i) + ".vtk").toStdString());

        double percentage1 = scar->Thresholding(thresh1);
        double percentage2 = scar->Thresholding(thresh2);
        double percentage3 = scar->Thresholding(thresh3);

        ofstream prodFile1;
        prodFile1.open((prodPath + "scar" + QString::number(i) + ".txt").toStdString());
        prodFile1 << threshType << "\n";
        prodFile1 << mean << "\n";
        prodFile1 << stdv << "\n";
        prodFile1 << "\n";
        prodFile1 << value1 << "\n";
        prodFile1 << thresh1 << "\n";
        prodFile1 << percentage1 << "\n";
        prodFile1 << "\n";
        prodFile1 << value2 << "\n";
        prodFile1 << thresh2 << "\n";
        prodFile1 << percentage2 << "\n";
        prodFile1 << "\n";
        prodFile1 << value3 << "\n";
        prodFile1 << thresh3 << "\n";
        prodFile1 << percentage3 << "\n";
        prodFile1.close();
    }

}

void autoScarCalculations() {

    //Test cases CEMRA and LGE (0)
    std::vector<int> testList = {184,141,239,167,245,142,139,209,192,189,225,260,182,187,135,191,255,257,149,168,250,198,180,197,146,131,119,213,177,227,157,107,95,89,101,98,91,100,104,99,105};

    //Image types
    typedef itk::Image<short, 3> ImageTypeSHRT;
    typedef itk::Image<short, 3> ImageTypeCHAR;
    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());

//    for (int i=0; i<testList.size(); i++) {
    for (int i=0; i<100; i++) {

        //Prepare dirctories (1)
        //QString direct = "/home/or15/syncdir/Proj/RZ/AutoAtrialSegmentation/CNN/CalcFibrosis/CEMRA5Observer/CNN/Case" + QString::number(testList.at(i));
        QString direct = "/home/or15/Desktop/Proj/RZ/AutoAtrialSegmentation/CNN/CalcFibrosis/PredictionsChallenge";

        //Get LGE and MRA image (2)
//        QString lgePath = direct + "/dcm-LGE.nii";
//        QString mraPath = direct + "/dcm-MRA.nii";
        QString cnnPath = direct + "/output_" + QString::number(i) + ".nii";

//        //Adjust CNN label to MRA space
//        mitk::Image::Pointer mraIMG = mitk::IOUtil::LoadImage(mraPath.toStdString());
//        mitk::Image::Pointer cnnIMG = mitk::IOUtil::LoadImage(cnnPath.toStdString());
//        double origin[3]; double spacing[3];
//        mraIMG->GetGeometry()->GetOrigin().ToArray(origin);
//        mraIMG->GetGeometry()->GetSpacing().ToArray(spacing);
//        vtkSmartPointer<vtkImageResize> resizeFilter = vtkSmartPointer<vtkImageResize>::New();
//        resizeFilter->SetResizeMethodToOutputDimensions();
//        resizeFilter->SetOutputDimensions(mraIMG->GetDimension(0), mraIMG->GetDimension(1), mraIMG->GetDimension(2));
//        resizeFilter->InterpolateOff();
//        resizeFilter->SetInputData(cnnIMG->GetVtkImageData());
//        resizeFilter->Update();
//        vtkSmartPointer<vtkImageChangeInformation> changeFilter = vtkSmartPointer<vtkImageChangeInformation>::New();
//        changeFilter->SetInputConnection(resizeFilter->GetOutputPort());
//        changeFilter->SetOutputSpacing(spacing);
//        changeFilter->SetOutputOrigin(origin);
//        changeFilter->Update();
//        mitk::Image::Pointer cnnLA = mitk::Image::New();
//        cnnIMG->Initialize(changeFilter->GetOutput());
//        cnnIMG->SetVolume(changeFilter->GetOutput()->GetScalarPointer());

//        //Image Registration (3)
//        mitk::IOUtil::Save(cnnIMG, (direct + "/LA.nii").toStdString());
//        cmd->ExecuteTransformation(direct, "LA.nii", "LA-reg.nii");

        //Clean the segmentation
        typedef itk::ImageRegionIteratorWithIndex<ImageTypeCHAR> ItType;
        ImageTypeCHAR::Pointer orgSegImage = ImageTypeCHAR::New();
        //CastToItkImage(mitk::IOUtil::LoadImage((direct + "/LA-reg.nii").toStdString()), orgSegImage);
        CastToItkImage(mitk::IOUtil::LoadImage(cnnPath.toStdString()), orgSegImage);
        typedef itk::ConnectedComponentImageFilter<ImageTypeCHAR, ImageTypeCHAR> ConnectedComponentImageFilterType;
        ConnectedComponentImageFilterType::Pointer connected1 = ConnectedComponentImageFilterType::New();
        connected1->SetInput(orgSegImage);
        connected1->Update();
        typedef itk::LabelShapeKeepNObjectsImageFilter<ImageTypeCHAR> LabelShapeKeepNObjImgFilterType;
        LabelShapeKeepNObjImgFilterType::Pointer lblShpKpNObjImgFltr1 = LabelShapeKeepNObjImgFilterType::New();
        lblShpKpNObjImgFltr1->SetInput(connected1->GetOutput());
        lblShpKpNObjImgFltr1->SetBackgroundValue(0);
        lblShpKpNObjImgFltr1->SetNumberOfObjects(1);
        lblShpKpNObjImgFltr1->SetAttribute(LabelShapeKeepNObjImgFilterType::LabelObjectType::NUMBER_OF_PIXELS);
        lblShpKpNObjImgFltr1->Update();
        using DuplicatorType = itk::ImageDuplicator<ImageTypeCHAR>;
        DuplicatorType::Pointer duplicator = DuplicatorType::New();
        duplicator->SetInputImage(lblShpKpNObjImgFltr1->GetOutput());
        duplicator->Update();
        ItType itDUP(duplicator->GetOutput(), duplicator->GetOutput()->GetRequestedRegion());
        for (itDUP.GoToBegin(); !itDUP.IsAtEnd(); ++itDUP)
            if ((int)itDUP.Get() != 0)
                itDUP.Set(1);
        //QString segCleanPath = direct + "/prodClean.nii";
        QString segCleanPath = direct + "/prodClean_" + QString::number(i) + ".nii";
        mitk::IOUtil::Save(mitk::ImportItkImage(duplicator->GetOutput()), segCleanPath.toStdString());

//        //Vein clipping mesh
//        QString output1 = cmd->ExecuteSurf(direct, segCleanPath, 1, .5, 0, 10);
//        mitk::Surface::Pointer shell = mitk::IOUtil::LoadSurface(output1.toStdString());
//        vtkSmartPointer<vtkDecimatePro> deci = vtkSmartPointer<vtkDecimatePro>::New();
//        deci->SetInputData(shell->GetVtkPolyData());
//        deci->SetTargetReduction(0.1);
//        deci->PreserveTopologyOn();
//        deci->Update();
//        shell->SetVtkPolyData(deci->GetOutput());
//        vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
//        vtkSmartPointer<vtkPolyData> pd = shell->Clone()->GetVtkPolyData();
//        for (int i=0; i<pd->GetNumberOfPoints(); i++) {
//            double* point = pd->GetPoint(i);
//            point[0] = -point[0];
//            point[1] = -point[1];
//            pd->GetPoints()->SetPoint(i, point);
//        }//_for
//        pointLocator->SetDataSet(pd);
//        pointLocator->BuildLocator();

//        //Separate veins
//        ImageTypeCHAR::Pointer veinsSegImage = ImageTypeCHAR::New();
//        veinsSegImage = lblShpKpNObjImgFltr1->GetOutput();
//        ItType itORG(orgSegImage, orgSegImage->GetRequestedRegion());
//        ItType itVEN(veinsSegImage, veinsSegImage->GetRequestedRegion());
//        itORG.GoToBegin();
//        for (itVEN.GoToBegin(); !itVEN.IsAtEnd(); ++itVEN) {
//            if ((int)itVEN.Get() != 0)
//                itVEN.Set((int)itORG.Get());
//            ++itORG;
//        }
//        for (itVEN.GoToBegin(); !itVEN.IsAtEnd(); ++itVEN)
//            if ((int)itVEN.Get() != 2)
//                itVEN.Set(0);
//        typedef itk::BinaryCrossStructuringElement<ImageTypeCHAR::PixelType, 3> CrossType;
//        typedef itk::BinaryMorphologicalOpeningImageFilter<ImageTypeCHAR, ImageTypeCHAR, CrossType> MorphFilterType;
//        CrossType binaryCross;
//        binaryCross.SetRadius(1.0);
//        binaryCross.CreateStructuringElement();
//        MorphFilterType::Pointer morphFilter = MorphFilterType::New();
//        morphFilter->SetInput(veinsSegImage);
//        morphFilter->SetKernel(binaryCross);
//        morphFilter->SetForegroundValue(2);
//        morphFilter->SetBackgroundValue(0);
//        morphFilter->UpdateLargestPossibleRegion();
//        veinsSegImage = morphFilter->GetOutput();
//        mitk::IOUtil::Save(mitk::ImportItkImage(veinsSegImage), (direct + "/prodVeins.nii").toStdString());
//        ConnectedComponentImageFilterType::Pointer connected2 = ConnectedComponentImageFilterType::New();
//        connected2->SetInput(veinsSegImage);
//        connected2->Update();
//        typedef itk::RelabelComponentImageFilter<ImageTypeCHAR, ImageTypeCHAR> RelabelFilterType;
//        RelabelFilterType::Pointer relabeler = RelabelFilterType::New();
//        relabeler->SetInput(connected2->GetOutput());
//        relabeler->Update();
//        mitk::IOUtil::Save(mitk::ImportItkImage(relabeler->GetOutput()), (direct + "/prodSepratedVeins.nii").toStdString());

//        //Find vein landmark
//        veinsSegImage = relabeler->GetOutput();
//        ItType itLMK(veinsSegImage, veinsSegImage->GetRequestedRegion());
//        vtkSmartPointer<vtkIdList> pickedSeedIds = vtkSmartPointer<vtkIdList>::New();
//        pickedSeedIds->Initialize();
//        std::vector<std::vector<double>> veinsCentre;
//        const int nveins = static_cast<int>(connected2->GetObjectCount());
//        for (int j=0; j<nveins; j++) {
//            int ctrVeinsVoxels = 0;
//            std::vector<double> veinLandmark(3, 0.0);
//            for (itLMK.GoToBegin(); !itLMK.IsAtEnd(); ++itLMK) {
//                if ((int)itLMK.Get() == (j+1)) {
//                    ImageTypeCHAR::PointType point;
//                    veinsSegImage->TransformIndexToPhysicalPoint(itLMK.GetIndex(), point);
//                    veinLandmark[0] += point[0];
//                    veinLandmark[1] += point[1];
//                    veinLandmark[2] += point[2];
//                    ctrVeinsVoxels++;
//                }
//            }//_for
//            veinLandmark[0] /= ctrVeinsVoxels;
//            veinLandmark[1] /= ctrVeinsVoxels;
//            veinLandmark[2] /= ctrVeinsVoxels;
//            veinsCentre.push_back(veinLandmark);
//        }
//        for (int j=0; j<nveins; j++) {
//            double veinLandmark[3];
//            veinLandmark[0] = veinsCentre.at(j)[0];
//            veinLandmark[1] = veinsCentre.at(j)[1];
//            veinLandmark[2] = veinsCentre.at(j)[2];
//            vtkIdType id = pointLocator->FindClosestPoint(veinLandmark);
//            pickedSeedIds->InsertNextId(id);
//        }
//        std::vector<int> pickedSeedLabels;
//        for (int j=0; j<nveins; j++)
//            pickedSeedLabels.push_back(21);

//        //Clip the veins
//        std::unique_ptr<CemrgAtriaClipper> clipper(new CemrgAtriaClipper(direct, shell));
//        clipper->ComputeCtrLines(pickedSeedLabels, pickedSeedIds, false);
//        clipper->ComputeCtrLinesClippers(pickedSeedLabels);
//        clipper->ClipVeinsImage(pickedSeedLabels, mitk::ImportItkImage(duplicator->GetOutput()), false);

//        //Create a mesh from clipped segmentation of veins
//        QString output2 = cmd->ExecuteSurf(direct, (direct + "/PVeinsCroppedImage.nii"), 1, .5, 0, 10);
//        mitk::Surface::Pointer LAShell = mitk::IOUtil::LoadSurface(output2.toStdString());

//        //Clip the mitral valve
//        ImageTypeCHAR::Pointer mvImage = ImageTypeCHAR::New();
//        CastToItkImage(mitk::IOUtil::LoadImage(segCleanPath.toStdString()), mvImage);
//        ItType itMVI1(mvImage, mvImage->GetRequestedRegion());
//        itORG.GoToBegin();
//        for (itMVI1.GoToBegin(); !itMVI1.IsAtEnd(); ++itMVI1) {
//            if ((int)itMVI1.Get() != 0)
//                itMVI1.Set((int)itORG.Get());
//            ++itORG;
//        }
//        for (itMVI1.GoToBegin(); !itMVI1.IsAtEnd(); ++itMVI1)
//            if ((int)itMVI1.Get() != 3)
//                itMVI1.Set(0);
//        typedef itk::ConnectedComponentImageFilter<ImageTypeCHAR, ImageTypeCHAR> ConnectedComponentImageFilterType;
//        ConnectedComponentImageFilterType::Pointer connected3 = ConnectedComponentImageFilterType::New();
//        connected3->SetInput(mvImage);
//        connected3->Update();
//        typedef itk::LabelShapeKeepNObjectsImageFilter<ImageTypeCHAR> LabelShapeKeepNObjImgFilterType;
//        LabelShapeKeepNObjImgFilterType::Pointer lblShpKpNObjImgFltr2 = LabelShapeKeepNObjImgFilterType::New();
//        lblShpKpNObjImgFltr2->SetInput(connected3->GetOutput());
//        lblShpKpNObjImgFltr2->SetBackgroundValue(0);
//        lblShpKpNObjImgFltr2->SetNumberOfObjects(1);
//        lblShpKpNObjImgFltr2->SetAttribute(LabelShapeKeepNObjImgFilterType::LabelObjectType::NUMBER_OF_PIXELS);
//        lblShpKpNObjImgFltr2->Update();
//        mvImage = lblShpKpNObjImgFltr2->GetOutput();
//        mitk::IOUtil::Save(mitk::ImportItkImage(mvImage), (direct + "/prodMVI.nii").toStdString());
//        ItType itMVI2(mvImage, mvImage->GetRequestedRegion());
//        ImageTypeCHAR::Pointer clnSegImage = ImageTypeCHAR::New();
//        CastToItkImage(mitk::IOUtil::LoadImage(segCleanPath.toStdString()), clnSegImage);
//        ItType itCLN(clnSegImage, clnSegImage->GetRequestedRegion());
//        itORG.GoToBegin();
//        itMVI2.GoToBegin();
//        for (itCLN.GoToBegin(); !itCLN.IsAtEnd(); ++itCLN) {
//            if ((int)itCLN.Get() != 0) {
//                if ((int)itORG.Get() == 3 && (int)itMVI2.Get() == 0) {
//                    ++itORG; ++itMVI2;
//                    continue;
//                }
//                itCLN.Set((int)itORG.Get());
//            }//_if
//            ++itORG; ++itMVI2;
//        }//_for
//        vtkSmartPointer<vtkIntArray> cutScores = vtkSmartPointer<vtkIntArray>::New();
//        for (vtkIdType ID = 0; ID<LAShell->GetVtkPolyData()->GetNumberOfPoints(); ID++) {
//            ImageTypeCHAR::PointType pointXYZ;
//            ImageTypeCHAR::IndexType voxelXYZ;
//            pointXYZ[0] = -LAShell->GetVtkPolyData()->GetPoint(ID)[0];
//            pointXYZ[1] = -LAShell->GetVtkPolyData()->GetPoint(ID)[1];
//            pointXYZ[2] =  LAShell->GetVtkPolyData()->GetPoint(ID)[2];
//            clnSegImage->TransformPhysicalPointToIndex(pointXYZ, voxelXYZ);
//            bool bodyLA = true;
//            ImageTypeCHAR::IndexType voxelNHG;
//            for (int a=-1; a<=1 && bodyLA; a++) {
//                for (int b=-1; b<=1 && bodyLA; b++) {
//                    for (int c=-1; c<=1 && bodyLA; c++) {
//                        voxelNHG[0] = voxelXYZ[0] + a;
//                        voxelNHG[1] = voxelXYZ[1] + b;
//                        voxelNHG[2] = voxelXYZ[2] + c;
//                        if (static_cast<int>(clnSegImage->GetPixel(voxelNHG)) == 3) bodyLA = false;
//                    }
//                }
//            }
//            if (!bodyLA && static_cast<int>(clnSegImage->GetPixel(voxelXYZ)) != 1) cutScores->InsertNextTuple1(3);
//            else cutScores->InsertNextTuple1(0);
//        }
//        LAShell->GetVtkPolyData()->GetPointData()->SetScalars(cutScores);
//        vtkSmartPointer<vtkClipPolyData> clip = vtkSmartPointer<vtkClipPolyData>::New();
//        clip->SetInputData(LAShell->GetVtkPolyData());
//        clip->InsideOutOn();
//        clip->Update();
//        vtkSmartPointer<vtkPolyDataConnectivityFilter> lrgRegion = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
//        lrgRegion->SetInputConnection(clip->GetOutputPort());
//        lrgRegion->SetExtractionModeToLargestRegion();
//        lrgRegion->Update();
//        vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
//        cleaner->SetInputConnection(lrgRegion->GetOutputPort());
//        cleaner->Update();
//        LAShell->SetVtkPolyData(cleaner->GetOutput());
//        mitk::IOUtil::Save(LAShell, (direct + "/segmentation.vtk").toStdString());

//        double mean = 0.0, stdv = 0.0;
//        for (int step=1; step<=6; step++) {

//            //Scar projection
//            int minStep = -1;
//            int maxStep = step;
//            int methodType = 2;
//            std::unique_ptr<CemrgScar3D> scar(new CemrgScar3D());
//            scar->SetMinStep(minStep);
//            scar->SetMaxStep(maxStep);
//            scar->SetMethodType(methodType);
//            ImageTypeCHAR::Pointer segITK = ImageTypeCHAR::New();
//            CastToItkImage(mitk::IOUtil::LoadImage((direct + "/PVeinsCroppedImage.nii").toStdString()), segITK);
//            ImageTypeSHRT::Pointer lgeITK = ImageTypeSHRT::New();
//            CastToItkImage(mitk::IOUtil::LoadImage(lgePath.toStdString()), lgeITK);
//            //itk::ResampleImageFilter<ImageTypeCHAR, ImageTypeCHAR>::Pointer resampleFilter;
//            //resampleFilter = itk::ResampleImageFilter<ImageTypeCHAR, ImageTypeCHAR>::New();
//            //resampleFilter->SetInput(segITK);
//            //resampleFilter->SetReferenceImage(lgeITK);
//            //resampleFilter->SetUseReferenceImage(true);
//            //resampleFilter->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<ImageTypeCHAR>::New());
//            //resampleFilter->SetDefaultPixelValue(0);
//            //resampleFilter->UpdateLargestPossibleRegion();
//            //segITK = resampleFilter->GetOutput();
//            //mitk::IOUtil::Save(mitk::ImportItkImage(segITK), (direct + "/PVeinsCroppedImage.nii").toStdString());
//            scar->SetScarSegImage(mitk::ImportItkImage(segITK));
//            mitk::Surface::Pointer scarShell = scar->Scar3D(direct.toStdString(), mitk::ImportItkImage(lgeITK));
//            mitk::IOUtil::Save(scarShell, (direct + "/scar_" + QString::number(step) + ".vtk").toStdString());

//            //Thresholding
//            int vxls = 3;
//            int threshType = 1;
//            double value1 = 0.97;
//            double value2 = 1.61;
//            double value3 = 3.3;

//            if (step==1) {
//                typedef itk::Image<float, 3> ImageType;
//                typedef itk::BinaryBallStructuringElement<ImageTypeCHAR::PixelType, 3> BallType;
//                typedef itk::GrayscaleErodeImageFilter<ImageTypeCHAR, ImageType, BallType> ErosionFilterType;
//                BallType binaryBall;
//                binaryBall.SetRadius(vxls);
//                binaryBall.CreateStructuringElement();
//                ErosionFilterType::Pointer erosionFilter = ErosionFilterType::New();
//                erosionFilter->SetInput(segITK);
//                erosionFilter->SetKernel(binaryBall);
//                erosionFilter->UpdateLargestPossibleRegion();
//                mitk::Image::Pointer roiImage = mitk::Image::New();
//                roiImage = mitk::ImportItkImage(erosionFilter->GetOutput())->Clone();
//                ImageType::Pointer lgeFloat = ImageType::New();
//                CastToItkImage(mitk::IOUtil::LoadImage(lgePath.toStdString()), lgeFloat);
//                scar->CalculateMeanStd(mitk::ImportItkImage(lgeFloat), roiImage, mean, stdv);
//            }

//            double thresh1 = mean*value1;
//            double thresh2 = mean*value2;
//            double thresh3 = mean+value3*stdv;
//            double percentage1 = scar->Thresholding(thresh1);
//            double percentage2 = scar->Thresholding(thresh2);
//            double percentage3 = scar->Thresholding(thresh3);
//            QString prodPath = direct + mitk::IOUtil::GetDirectorySeparator();
//            ofstream prodFile1;
//            prodFile1.open((prodPath + "prodThresholds_" + QString::number(step) + ".txt").toStdString());
//            prodFile1 << threshType << "\n";
//            prodFile1 << mean << "\n";
//            prodFile1 << stdv << "\n";
//            prodFile1 << "\n";
//            prodFile1 << value1 << "\n";
//            prodFile1 << thresh1 << "\n";
//            prodFile1 << percentage1 << "\n";
//            prodFile1 << "\n";
//            prodFile1 << value2 << "\n";
//            prodFile1 << thresh2 << "\n";
//            prodFile1 << percentage2 << "\n";
//            prodFile1 << "\n";
//            prodFile1 << value3 << "\n";
//            prodFile1 << thresh3 << "\n";
//            prodFile1 << percentage3 << "\n";
//            prodFile1.close();

//        }//_for_steps
    }
}

void autoMeshDeformer() {

    QString dirc;
    QString path;
    QString dofi;
    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());

    for (int i=33; i<=33; i++) {
        //dirc = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_5/HighDose5/Test" + QString::number(i);
        //path = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_5/HighDose5/Test" + QString::number(i) + "/Original/LV-0.nii";
        //cmd->ExecuteSurf(dirc, path, 1, 0.5, 0, 10);
        //dirc = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_10/HighDose10/Test" + QString::number(i);
        //path = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_10/HighDose10/Test" + QString::number(i) + "/Original/LV-0.nii";
        //cmd->ExecuteSurf(dirc, path, 1, 0.5, 0, 10);

        dirc = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_5/HighDose5/Test" + QString::number(i);
        path = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_5/HighDose5/Test" + QString::number(i) + "/segmentation.vtk";
        dofi = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_5/HighDose5/Test" + QString::number(i) + "/tsffd.dof";
        cmd->ExecuteApplying(dirc, path, 0, dofi, 20, 1);
        dirc = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_5/LowDose5/Test" + QString::number(i);
        path = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_5/HighDose5/Test" + QString::number(i) + "/segmentation.vtk";
        dofi = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_5/HighDose5/Test" + QString::number(i) + "/rigid.dof";
        cmd->ExecuteApplying(dirc, path, 0, dofi, 1, 1);
        std::rename((dirc + "/transformed-0.vtk").toStdString().c_str() , (dirc + "/segmentation_reg.vtk").toStdString().c_str());
        dirc = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_5/LowDose5/Test" + QString::number(i);
        path = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_5/LowDose5/Test" + QString::number(i) + "/segmentation_reg.vtk";
        dofi = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_5/LowDose5/Test" + QString::number(i) + "/tsffd.dof";
        cmd->ExecuteApplying(dirc, path, 0, dofi, 20, 1);

        //dirc = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_10/HighDose10/Test" + QString::number(i);
        //path = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_10/HighDose10/Test" + QString::number(i) + "/segmentation.vtk";
        //dofi = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_10/HighDose10/Test" + QString::number(i) + "/tsffd.dof";
        //cmd->ExecuteApplying(dirc, path, 0, dofi, 10, 1);
        //dirc = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_10/LowDose10/Test" + QString::number(i);
        //path = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_10/HighDose10/Test" + QString::number(i) + "/segmentation.vtk";
        //dofi = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_10/HighDose10/Test" + QString::number(i) + "/rigid.dof";
        //cmd->ExecuteApplying(dirc, path, 0, dofi, 1, 1);
        //std::rename((dirc + "/transformed-0.vtk").toStdString().c_str() , (dirc + "/segmentation_reg.vtk").toStdString().c_str());
        //dirc = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_10/LowDose10/Test" + QString::number(i);
        //path = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_10/LowDose10/Test" + QString::number(i) + "/segmentation_reg.vtk";
        //dofi = "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_10/LowDose10/Test" + QString::number(i) + "/tsffd.dof";
        //cmd->ExecuteApplying(dirc, path, 0, dofi, 10, 1);
    }
}

void downsampleImages() {

    typedef itk::Image<short,3> ImageType;
    typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleImageFilterType;
    typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> NearestInterpolatorType;

    for (int i=34; i<=34; i++) {
        for (int j=0; j<=19; j++) {
            std::string path =
                    "/home/or15/Desktop/ToDo/Pigs/Tests/Incr_5/LowDose5/Test" +
                    QString::number(i).toStdString() + "/dcm-" +
                    QString::number(j).toStdString() + ".nii";
            mitk::Image::Pointer image = mitk::IOUtil::LoadImage(path);
            ImageType::Pointer itkImage = ImageType::New();
            mitk::CastToItkImage(image, itkImage);
            ResampleImageFilterType::Pointer downsampler = ResampleImageFilterType::New();
            downsampler->SetInput(itkImage);
            NearestInterpolatorType::Pointer interpolator = NearestInterpolatorType::New();
            downsampler->SetInterpolator(interpolator);
            downsampler->SetDefaultPixelValue(0);
            ResampleImageFilterType::SpacingType spacing = itkImage->GetSpacing();
            downsampler->SetOutputSpacing(spacing);
            downsampler->SetOutputOrigin(itkImage->GetOrigin());
            downsampler->SetOutputDirection(itkImage->GetDirection());
            ResampleImageFilterType::SizeType size = itkImage->GetLargestPossibleRegion().GetSize();
            for (int i=0; i<3; ++i)
                size[i] /= 1.21;
            downsampler->SetSize(size);
            downsampler->UpdateLargestPossibleRegion();
            image = mitk::ImportItkImage(downsampler->GetOutput())->Clone();
            mitk::IOUtil::SaveImage(image, path);
        }
    }
}

void convertfromDICOM() {

    for (int j=33; j<=36; j++) {
        for (int i=1; i<=21; i++) {
            std::string intPath = "/home/or15/Desktop/ToDo/Pigs/Data/" + QString::number(j).toStdString() + "/_retro_5%/LowDose/" + QString::number(i).toStdString() + "/";
            DIR* dirp = opendir(intPath.c_str()); struct dirent * dp;
            dp = readdir(dirp); dp = readdir(dirp); dp = readdir(dirp);
            intPath = intPath + dp->d_name; closedir(dirp);
            mitk::Image::Pointer image = mitk::IOUtil::LoadImage(intPath);
            std::string outPath = "/home/or15/Desktop/ToDo/Pigs/Tests/LowDose5/Test" + QString::number(j).toStdString() + "/dcm-" + QString::number(i-1).toStdString() + ".nii";
            mitk::IOUtil::SaveImage(image, outPath);
        }
    }
}

void createMVImagesVTK() {

    constexpr unsigned int Dimension = 3;
    using ImageType = itk::Image<unsigned char, Dimension>;

    std::string direct, line1, line2;
    ifstream file("/home/or15/Desktop/Proj/RZ/AutoAtrialSegmentation/CNN/ProjectFiles/DIRECT_paths.txt");

    int ctr = 0;
    if (file.is_open()) {
        while (getline(file,direct)) {

            //Skip cases
            ctr = ctr + 1;
            if (ctr<243) continue;

            //Adjust folder path
            std::string expr = "/data/or15/AutoAtrialSegmentation/";
            direct = direct.replace(direct.find(expr), expr.length(), "/home/or15/Desktop/Proj/RZ/AutoAtrialSegmentation/CNN/").c_str();
            line1 = direct + "/BP_VN.nii";
            line2 = direct + "/prodMVCLandmarks.mps";

            //Load in the image
            mitk::Image::Pointer segImage = mitk::IOUtil::LoadImage(line1);
            ImageType::Pointer segItkImage = ImageType::New();
            CastToItkImage(segImage, segItkImage);

            //Load in the pointset
            mitk::PointSet::Pointer pointSet = mitk::IOUtil::LoadPointSet(line2);
            double x_c = 0;
            double y_c = 0;
            double z_c = 0;
            for(int i=0; i<pointSet->GetSize(); i++) {
                x_c += pointSet->GetPoint(i).GetElement(0);
                y_c += pointSet->GetPoint(i).GetElement(1);
                z_c += pointSet->GetPoint(i).GetElement(2);
            }//_for
            x_c /= pointSet->GetSize();
            y_c /= pointSet->GetSize();
            z_c /= pointSet->GetSize();
            double distance[pointSet->GetSize()];
            for(int i=0; i<pointSet->GetSize(); i++) {
                double x_d = pointSet->GetPoint(i).GetElement(0) - x_c;
                double y_d = pointSet->GetPoint(i).GetElement(1) - y_c;
                double z_d = pointSet->GetPoint(i).GetElement(2) - z_c;
                distance[i] = sqrt(pow(x_d,2) + pow(y_d,2) + pow(z_d,2));
            }//_for
            double radius = *std::max_element(distance, distance + pointSet->GetSize());

            //Create the sphere
            vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
            sphere->SetCenter(x_c, y_c, z_c);
            sphere->SetRadius(radius);
            sphere->SetThetaResolution(100);
            sphere->Update();

            //Prepare empty image
            vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();
            double spacing[3];
            segImage->GetVtkImageData()->GetSpacing(spacing);
            whiteImage->SetSpacing(spacing);
            int dimensions[3];
            segImage->GetVtkImageData()->GetDimensions(dimensions);
            whiteImage->SetDimensions(dimensions);
            whiteImage->SetExtent(0, dimensions[0] - 1, 0, dimensions[1] - 1, 0, dimensions[2] - 1);
            double origin[3];
            segImage->GetGeometry()->GetOrigin().ToArray(origin);
            whiteImage->SetOrigin(origin);
            whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);
            unsigned char otval = 0;
            unsigned char inval = 255;
            vtkIdType count = whiteImage->GetNumberOfPoints();
            for (vtkIdType i = 0; i < count; ++i)
                whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);

            //Poly to image stencil
            vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
            pol2stenc->SetInputData(sphere->GetOutput());
            pol2stenc->SetOutputOrigin(origin);
            pol2stenc->SetOutputSpacing(spacing);
            pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
            pol2stenc->Update();

            //Set white image and set the background
            vtkSmartPointer<vtkImageStencil> imgstenc = vtkSmartPointer<vtkImageStencil>::New();
            imgstenc->SetInputData(whiteImage);
            imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
            imgstenc->ReverseStencilOff();
            imgstenc->SetBackgroundValue(otval);
            imgstenc->Update();

            //VTK to ITK conversion
            mitk::Image::Pointer sphereImg = mitk::Image::New();
            sphereImg->Initialize(imgstenc->GetOutput());
            sphereImg->SetVolume(imgstenc->GetOutput()->GetScalarPointer());
            ImageType::Pointer sphereItkImage = ImageType::New();
            CastToItkImage(sphereImg, sphereItkImage);
            itk::ResampleImageFilter<ImageType, ImageType>::Pointer resampleFilter;
            resampleFilter = itk::ResampleImageFilter<ImageType, ImageType >::New();
            resampleFilter->SetInput(sphereItkImage);
            resampleFilter->SetReferenceImage(segItkImage);
            resampleFilter->SetUseReferenceImage(true);
            resampleFilter->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<ImageType>::New());
            resampleFilter->SetDefaultPixelValue(0);
            resampleFilter->UpdateLargestPossibleRegion();
            sphereItkImage = resampleFilter->GetOutput();

            //Create the intersection
            mitk::Image::Pointer result;
            mitk::Image::Pointer LA = mitk::ImportItkImage(segItkImage);
            mitk::Image::Pointer MV = mitk::ImportItkImage(sphereItkImage);
            mitk::BooleanOperation booleanOperation(mitk::BooleanOperation::Intersection, LA, MV);
            result = booleanOperation.GetResult();

            //Write to file
            std::size_t found = line1.find_last_of("/");
            std::string path = line1.substr(0,found) + "/" + "MVimg.nii";
            mitk::IOUtil::Save(result, path);
            qDebug() << path.c_str();
        }
    }
}

void calcMorphAuto() {

    std::string line;
    ifstream file("/home/or15/Desktop/Proj/MD/CEMRGResults/files.txt");
    std::unique_ptr<CemrgCommandLine> cmd(new CemrgCommandLine());
    if (file.is_open()) {
        while (getline (file,line)) {
            QString directory = QString::fromStdString(line);
            //Start analysing automatically
            QString path = directory + mitk::IOUtil::GetDirectorySeparator() + "AnalyticBloodpool.nii";
            mitk::Image::Pointer analyticImage = mitk::IOUtil::LoadImage(path.toStdString());
            if (analyticImage) {
                //Loop through labelled image
                typedef itk::Image<short, 3> ImageType;
                typedef itk::ImageRegionIteratorWithIndex<ImageType> ItType;
                ImageType::Pointer analyticItkImage = ImageType::New();
                CastToItkImage(analyticImage, analyticItkImage);
                ItType itLbl(analyticItkImage, analyticItkImage->GetRequestedRegion());
                for (itLbl.GoToBegin(); !itLbl.IsAtEnd(); ++itLbl)
                    if ((int)itLbl.Get() == 19 || (int)itLbl.Get() == 20)
                        itLbl.Set(0);
                //Relabel the components to separate bloodpool and appendage based on the number of voxels
                typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> ConnectedComponentImageFilterType;
                ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
                connected->SetInput(analyticItkImage);
                connected->Update();
                typedef itk::RelabelComponentImageFilter<ImageType, ImageType> RelabelFilterType;
                RelabelFilterType::Pointer relabeler = RelabelFilterType::New();
                relabeler->SetInput(connected->GetOutput());
                relabeler->Update();
                //Keep the selected labels
                const int bpLabel = 1;
                const int apLabel = 2;
                typedef itk::ThresholdImageFilter<ImageType> ThresholdImageFilterType;
                ThresholdImageFilterType::Pointer thresholdFilter1 = ThresholdImageFilterType::New();
                thresholdFilter1->SetInput(analyticItkImage);
                thresholdFilter1->ThresholdOutside(bpLabel, bpLabel);
                thresholdFilter1->SetOutsideValue(0);
                ThresholdImageFilterType::Pointer thresholdFilter2 = ThresholdImageFilterType::New();
                thresholdFilter2->SetInput(relabeler->GetOutput());
                thresholdFilter2->ThresholdOutside(apLabel, apLabel);
                thresholdFilter2->SetOutsideValue(0);
                //Import to MITK image
                mitk::Image::Pointer bp = mitk::ImportItkImage(thresholdFilter1->GetOutput());
                mitk::Image::Pointer ap = mitk::ImportItkImage(thresholdFilter2->GetOutput());
                int iter = 1;
                float th = 0.5;
                int blur = 0;
                int smth = 0;
                float ds = 0.99;
                QString path1 = directory + mitk::IOUtil::GetDirectorySeparator() + "bp.nii";
                mitk::IOUtil::Save(bp, path1.toStdString());
                QString output = cmd->ExecuteSurf(directory, path1, iter, th, blur, smth);
                mitk::Surface::Pointer shell = mitk::IOUtil::LoadSurface(output.toStdString());
                vtkSmartPointer<vtkDecimatePro> deci = vtkSmartPointer<vtkDecimatePro>::New();
                deci->SetInputData(shell->GetVtkPolyData());
                deci->SetTargetReduction(ds);
                deci->PreserveTopologyOn();
                deci->Update();
                shell->SetVtkPolyData(deci->GetOutput());
                mitk::Surface::Pointer surfLA = shell->Clone();
                QString path2 = directory + mitk::IOUtil::GetDirectorySeparator() + "ap.nii";
                mitk::IOUtil::Save(ap, path2.toStdString());
                output = cmd->ExecuteSurf(directory, path2, iter, th, blur, smth);
                shell = mitk::IOUtil::LoadSurface(output.toStdString());
                deci->SetInputData(shell->GetVtkPolyData());
                deci->SetTargetReduction(ds);
                deci->PreserveTopologyOn();
                deci->Update();
                shell->SetVtkPolyData(deci->GetOutput());
                mitk::Surface::Pointer surfAP = shell->Clone();
                //Volume and surface calculations
                std::unique_ptr<CemrgMeasure> morphAnal = std::unique_ptr<CemrgMeasure>(new CemrgMeasure());
                double surfceLA = morphAnal->calcSurfaceMesh(surfLA);
                double volumeLA = morphAnal->calcVolumeMesh(surfLA);
                double surfceAP = morphAnal->calcSurfaceMesh(surfAP);
                double volumeAP = morphAnal->calcVolumeMesh(surfAP);
                double sphereci = morphAnal->GetSphericity(surfLA->GetVtkPolyData());
                //Store in text file
                ofstream morphResult;
                QString morphPath = directory + mitk::IOUtil::GetDirectorySeparator() + "morphResultsSV.txt";
                morphResult.open(morphPath.toStdString(), std::ios_base::app);
                morphResult << "SA" << " " << surfceLA << "\n";
                morphResult << "VA" << " " << volumeLA << "\n";
                morphResult << "SP" << " " << surfceAP << "\n";
                morphResult << "VP" << " " << volumeAP << "\n";
                morphResult << "SC" << " " << sphereci << "\n";
                morphResult.close();
            }
            qDebug() << line.c_str();
        }
        file.close();
    }//_if
}

void create4DImagesfromDICOM() {

    for (int i=2; 2<=36; i++) {
        std::string intPath = "/home/or15/Desktop/Proj/BS/AtrialMotion/" + QString::number(i).toStdString() + "/Mri Heart/";
        DIR* dirp = opendir(intPath.c_str()); struct dirent * dp;
        dp = readdir(dirp); dp = readdir(dirp); dp = readdir(dirp);
        intPath = intPath + dp->d_name; closedir(dirp);
        mitk::Image::Pointer image = mitk::IOUtil::LoadImage(intPath);
        std::string outPath = "/home/or15/Desktop/Proj/BS/AtrialMotion/" + QString::number(i).toStdString() + "/4D.nii";
        mitk::IOUtil::SaveImage(image, outPath);
    }
}

//#include <mitkIDICOMTagsOfInterest.h>
mitk::DataNode::Pointer testDICOMtagOFinterest(QList<mitk::DataNode::Pointer> nodes) {

    mitk::DataNode::Pointer node = nodes.at(0);
    mitk::BaseData::Pointer data = node->GetData();
    mitk::PropertyList* list = data->GetPropertyList();
    mitk::BaseProperty* property = list->GetProperty("files");
    mitk::StringLookupTableProperty::Pointer filesProp = dynamic_cast<mitk::StringLookupTableProperty *>(property);
    if (filesProp.IsNull()) {
        mitkThrow() << "No property with dicom file path.";
        return NULL;
    }
    mitk::StringLookupTable filesLut = filesProp->GetValue();
    const mitk::StringLookupTable::LookupTableType& lookUpTableMap = filesLut.GetLookupTable();

    QStringList seriesToLoad;
    for (auto iterator : lookUpTableMap) {
        seriesToLoad.push_back(QString::fromStdString(iterator.second));
    }

    //    mitk::IDICOMTagsOfInterest* result = nullptr;
    //    std::vector<us::ServiceReference<mitk::IDICOMTagsOfInterest> > toiRegisters = us::GetModuleContext()->GetServiceReferences<mitk::IDICOMTagsOfInterest>();
    //    if (!toiRegisters.empty()) {
    //        if (toiRegisters.size() > 1)
    //            MITK_WARN << "Multiple DICOM tags of interest services found. Using just one.";
    //        result = us::GetModuleContext()->GetService<mitk::IDICOMTagsOfInterest>(toiRegisters.front());
    //    }

    std::vector<mitk::BaseData::Pointer> baseDatas = mitk::IOUtil::Load(seriesToLoad.front().toStdString());
    for (const auto &data : baseDatas) {

        mitk::DataNode::Pointer node = mitk::DataNode::New();
        node->SetData(data);
        return node;
    }

    //    IDICOMTagsOfInterest *toiService = GetDicomTagsOfInterestService();
    //    if (toiService != nullptr) {
    //        toiService->AddTagOfInterest(DICOMSegmentationConstants::SEGMENT_SEQUENCE_PATH());
    //    }
}

//int ctr = -5;
void CreateSyntImage() {

    //    //Prepare empty image
    //    vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();
    //    double spacing[3] = {0.83,0.83,0.83};
    //    whiteImage->SetSpacing(spacing);
    //    int dimensions[3] = {25,25,25};
    //    whiteImage->SetDimensions(dimensions);
    //    int extents[3] = {dimensions[0]-1, dimensions[1]-1, dimensions[2]-1};
    //    whiteImage->SetExtent(0, extents[0], 0, extents[1], 0, extents[2]);
    //    double origin[3] = {0,0,0};
    //    whiteImage->SetOrigin(origin);
    //    whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);
    //    //Colouring
    //    for (vtkIdType i = 0; i < whiteImage->GetNumberOfPoints(); ++i)
    //        whiteImage->GetPointData()->GetScalars()->SetTuple1(i, 0);
    //    //Tracking test
    //    ctr+=5;
    //    for(int x=10; x<15; x++) {
    //        for(int y=10+ctr; y<15+ctr; y++) {
    //            for(int z=10+ctr; z<15+ctr; z++) {
    //                vtkIdType idx = x + dimensions[0] * (y + dimensions[1] * z);
    //                whiteImage->GetPointData()->GetScalars()->SetTuple1(idx, 50);
    //            }
    //        }
    //    }
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
    //    mitk::Image::Pointer image = mitk::Image::New();
    //    image->Initialize(whiteImage);
    //    image->SetVolume(whiteImage->GetScalarPointer());
    //    mitk::IOUtil::Save(image, "/home/or15/Desktop/test.nii");
}

void ColourMesh() {

    for (int m=0; m<10; m++) {
        mitk::Surface::Pointer mesh;
        mesh = mitk::IOUtil::LoadSurface("/home/or15/Downloads/SqzMeshes/Mesh" + QString::number(m).toStdString() + ".vtk");
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
        for (int i=0; i<mesh->GetVtkPolyData()->GetNumberOfCells(); i++) {
            double sqz = mesh->GetVtkPolyData()->GetCellData()->GetScalars()->GetTuple1(i);
            double dcolor[3];
            colorLookupTable->GetColor(sqz, dcolor);
            unsigned char color[3];
            for(unsigned int j = 0; j < 3; j++)
                color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
            colors->InsertNextTupleValue(color);
        }//_for

        mesh->GetVtkPolyData()->GetCellData()->SetScalars(colors);
        mitk::IOUtil::SaveSurface(mesh, "/home/or15/Downloads/SqzMeshes/MeshColoured" + QString::number(m).toStdString() + ".vtk");
    }//_for
}

void binariseImage() {

    typedef itk::Image<double, 3> ImageType;
    mitk::Image::Pointer image = mitk::IOUtil::LoadImage("/home/or15/Downloads/wallThicknessNRRD/leftAtrialWall900101_UPsampled.nrrd");
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

    mitk::IOUtil::SaveImage(mitk::ImportItkImage(imageITK), "/home/or15/Downloads/wallThicknessNRRD/leftAtrialWall900101_UPsampled_mask.nrrd");
}

void MeshError() {

    for (int i=1; i<=18; i++) {
        if (i==4) continue;
        QString pathImage = "/home/or15/Desktop/Proj/RZ/StrainsWork/Analysis/MeshError/C" + QString::number(i) + "/segmentation.nii";
        typedef itk::Image<unsigned char, 3> ImageType;
        typedef itk::ImageRegionIteratorWithIndex<ImageType> ItType;
        ImageType::Pointer segItkImage = ImageType::New();
        CastToItkImage(mitk::IOUtil::LoadImage(pathImage.toStdString()), segItkImage);
        typedef itk::SimpleContourExtractorImageFilter<ImageType, ImageType> SimpleContourExtractorImageFilterType;
        SimpleContourExtractorImageFilterType::Pointer contourFilter = SimpleContourExtractorImageFilterType::New();
        contourFilter->SetInput(segItkImage);
        contourFilter->SetInputBackgroundValue(0);
        contourFilter->SetInputForegroundValue(1);
        contourFilter->Update();
        ImageType::Pointer ctrItkImage = contourFilter->GetOutput();
        QString pathMesh = "/home/or15/Desktop/Proj/RZ/StrainsWork/Analysis/MeshError/C" + QString::number(i) + "/segmentation.vtk";
        mitk::Surface::Pointer surface = mitk::IOUtil::LoadSurface(pathMesh.toStdString());
        vtkSmartPointer<vtkPolyData> pd = surface->GetVtkPolyData();
        for (int i=0; i<pd->GetNumberOfPoints(); i++) {
            double* point = pd->GetPoint(i);
            point[0] = -point[0];
            point[1] = -point[1];
            pd->GetPoints()->SetPoint(i, point);
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
                distances.push_back(sqrt(pow(x_d,2) + pow(y_d,2) + pow(z_d,2)));
            }//_if
        }//_for

        ofstream file;
        QString pathOutput = "/home/or15/Desktop/Proj/RZ/StrainsWork/Analysis/MeshError/C" + QString::number(i) + "/distances.csv";
        file.open(pathOutput.toStdString());
        for (size_t j=0; j<distances.size(); j++) {
            file << distances.at(j);
            if (j == distances.size()-1) file << endl;
            else file << ",";
        }
        qDebug() << "Case" << i;
    }//_i_18
}

void curvesCalculator() {

    for (int i=1; i<=24; i++) {
        if (i==4) continue;

        QString directory = "/home/or15/Desktop/Proj/RZ/StrainsWork/Strains/DEEDS/C" + QString::number(i);
        if (i>18) directory = "/home/or15/Desktop/Proj/RZ/StrainsWork/Strains/DEEDS/N" + QString::number(i-18);
        mitk::DataNode::Pointer lmNode = mitk::IOUtil::LoadDataNode((directory + "/PointSet.mps").toStdString());
        int segRatios[3] = {40, 40, 20};

        std::unique_ptr<CemrgStrains> strain;
        strain = std::unique_ptr<CemrgStrains>(new CemrgStrains(directory, 0));
        strain->ReferenceAHA(lmNode, segRatios);
        std::vector<std::vector<double>> plotValueVectorsSQZ;
        std::vector<std::vector<double>> plotValueVectorsCRC;
        std::vector<std::vector<double>> plotValueVectorsLNG;

        for (int j=0; j<10; j++) {
            plotValueVectorsSQZ.push_back(strain->CalculateSqzPlot(j));
            plotValueVectorsCRC.push_back(strain->CalculateStrainsPlot(j, lmNode, 3));
            plotValueVectorsLNG.push_back(strain->CalculateStrainsPlot(j, lmNode, 4));
        }

        for (int j=0; j<3; j++) {

            QString fileName;
            std::vector<std::vector<double>> plotValueVectors;
            if (j==0) {
                fileName = "SQZ.csv";
                plotValueVectors = plotValueVectorsSQZ;
            } else if (j==1) {
                fileName = "CRC.csv";
                plotValueVectors = plotValueVectorsCRC;
            } else {
                fileName = "LNG.csv";
                plotValueVectors = plotValueVectorsLNG;
            }//_if
            ofstream file;
            file.open(directory.toStdString() + mitk::IOUtil::GetDirectorySeparator() + fileName.toStdString());

            std::vector<double> values;
            for (int s=0; s<16; s++) {
                for (int f=0; f<10; f++)
                    values.push_back(plotValueVectors[f][s]);
                //Append the curve to the file
                for (size_t z=0; z<values.size(); z++) {
                    file << values.at(z);
                    if (z == values.size()-1) file << endl;
                    else file << ",";
                }
                values.clear();
            }//_for
            file.close();
        }
        qDebug() << "CASE" << i << "done!";
    }
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

void flattening() {

    mitk::Surface::Pointer surface = mitk::IOUtil::LoadSurface("/home/or15/Downloads/LGECART30MIN-Scar.vtk");

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
    } catch(itk::ExceptionObject & excp) {
    }
    MeshType::Pointer mesh = filter->GetOutput();

    surface->SetVtkPolyData(MeshUtil<MeshType>::MeshToPolyData(mesh));
    mitk::IOUtil::Save(surface, "/home/or15/Downloads/test.vtk");
}

void RRcalcsAuto() {

    std::vector<double> values;
    std::unique_ptr<CemrgMeasure> rr(new CemrgMeasure());
    std::vector <std::tuple<double, double, double>> points;
    std::string lineD;

    //1
    ifstream dirc("/home/or15/Desktop/Proj/RZ/StrainsWork/TrackTSFFD/Transforms/dirPaths.txt");
    std::vector<int> DS = {1,3,7,11,14,15};

    if (dirc.is_open()) {
        for (int ds=0; ds<DS.size(); ds++) {
            dirc.clear();
            dirc.seekg(0, ios::beg);
            while (getline(dirc,lineD)) {

                //2
                std::string prePath;
                auto const posSL2 = lineD.find_last_of('/');
                auto const posSL1 = lineD.find_last_of('/',posSL2-10);
                prePath = std::string(lineD.substr(0, posSL1+1).c_str()) + "ValidationSet/Dataset" + std::to_string(DS.at(ds)) + lineD.substr(posSL2);

                //3
                for (int i=0; i<1; i++) {
                    for (int j=0; j<1; j++) {

                        //4
                        std::stringstream stream;
                        stream << "/Perm-BEO-SWO";// << std::fixed << std::setprecision(1) << j/10.0;
                        //stream << "/Perm-BE" << std::fixed << std::setprecision(2) << i/100.0 << "-SW" << std::fixed << std::setprecision(2) << j/100.0;
                        std::string permStr = stream.str();

                        for (int frame=0; frame<10; frame++) {
                            QString path = QString::fromStdString(prePath) + QString::fromStdString(permStr);
                            qDebug() << path;
                            points = rr->Deconvert(path, frame);
                            if (QString::fromStdString(prePath).endsWith("PP",Qt::CaseInsensitive))
                                values.push_back(rr->CalcPerimeter(points));
                            else if (QString::fromStdString(prePath).endsWith("MAA",Qt::CaseInsensitive))
                                values.push_back(rr->CalcArea(points));
                            else
                                values.push_back(rr->CalcDistance(points));
                        }
                        ofstream file;
                        QString path = QString::fromStdString(prePath) + QString::fromStdString(permStr) + ".csv";
                        file.open(path.toStdString());
                        for (size_t z=0; z<values.size(); z++) {
                            file << values.at(z);
                            if (z == values.size()-1) file << endl;
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
    ifstream file("/home/or15/Desktop/Proj/Tom/Points/Tools/filPaths.txt");
    ifstream dirc("/home/or15/Desktop/Proj/Tom/Points/Tools/dirPaths.txt");
    if (file.is_open() && dirc.is_open()) {
        while (getline(file,lineF)) {

            getline(dirc,lineD);
            mitk::PointSet::Pointer MIPS = mitk::IOUtil::LoadPointSet(lineF);
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

    std::string inputFilename = "/home/or15/Downloads/transformed-0.vtk";
    std::string otputFilename = "/home/or15/Downloads/transformed-0.vtk";

    mitk::BaseData::Pointer meshData = mitk::IOUtil::Load(inputFilename).at(0);
    mitk::UnstructuredGrid::Pointer mitkVtkGrid = dynamic_cast<mitk::UnstructuredGrid*>(meshData.GetPointer());
    vtkSmartPointer<vtkUnstructuredGrid> vtkGrid = mitkVtkGrid->GetVtkUnstructuredGrid();
    //for (vtkIdType i=0; i<vtkGrid->GetNumberOfPoints(); i++) {
    //    double* point = vtkGrid->GetPoint(i);
    //    point[0] = -point[0];
    //    point[1] = -point[1];
    //    point[2] = -point[2];
    //    vtkGrid->GetPoints()->SetPoint(i, point);
    //}//_for
    //vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    //writer->SetInputData(vtkGrid);
    //writer->SetFileName(otputFilename.c_str());
    //writer->Write();
    vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
    geometryFilter->SetInputData(vtkGrid);
    geometryFilter->Update();
    vtkPolyData* polydata = geometryFilter->GetOutput();
    mitk::Surface::Pointer surf = mitk::Surface::New();
    surf->SetVtkPolyData(polydata);
    mitk::IOUtil::Save(surf, otputFilename);

    /**
    inputFilename = "/home/or15/Desktop/Proj/Sander/ControlCMR/VTKs/transformed-5.vtk";
    otputFilename = "/home/or15/Desktop/Proj/Sander/ControlCMR/Results/testest5.vtk";
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
        for (int i=0; i<output->GetNumberOfPoints(); i++) {
            double* point = output->GetPoint(i);
            point[0] = -point[0];
            point[1] = -point[1];
            output->GetPoints()->SetPoint(i, point);
        }//_for

        // Write out
        vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();
        double bounds[6];
        output->GetBounds(bounds);

        double voxel=1;
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
        whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);

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

        mitk::Image::Pointer manIm = mitk::IOUtil::LoadImage("/home/or15/Desktop/Proj/JB/OrodNormals/n1/Manuals/man-segmentation.nii");
        mitk::Image::Pointer autIm = mitk::IOUtil::LoadImage("/home/or15/Desktop/Proj/JB/OrodNormals/n1/Manuals/aut-segmentation.nii");

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
        for (int i=0; i<11; i++) {
            mitk::BaseData::Pointer meshData = mitk::IOUtil::Load(
                        "/home/or15/Desktop/Proj/JB/TomTec/transformed-"+
                        QString::number(i).toStdString()+".vtk").at(0);
            mitk::UnstructuredGrid::Pointer mitkVtkGrid = dynamic_cast<mitk::UnstructuredGrid*>(meshData.GetPointer());
            vtkSmartPointer<vtkUnstructuredGrid> vtkGrid = mitkVtkGrid->GetVtkUnstructuredGrid();

            vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
            surfaceFilter->SetInputData(vtkGrid);
            surfaceFilter->Update();
            vtkSmartPointer<vtkPolyData> polydata = surfaceFilter->GetOutput();
            mitk::Surface::Pointer converted = mitk::Surface::New();
            converted->SetVtkPolyData(polydata);
            mitk::IOUtil::Save(converted,
                               "/home/or15/Desktop/Proj/JB/TomTec/transformed-"+
                               QString::number(i).toStdString()+".vtk");
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
        std::tuple<double,double,double> set0(s0->GetPoint(0).GetElement(0),s0->GetPoint(0).GetElement(1),s0->GetPoint(0).GetElement(2));
        std::tuple<double,double,double> set1(s1->GetPoint(0).GetElement(0),s1->GetPoint(0).GetElement(1),s1->GetPoint(0).GetElement(2));
        std::tuple<double,double,double> set2(s2->GetPoint(0).GetElement(0),s2->GetPoint(0).GetElement(1),s2->GetPoint(0).GetElement(2));
        std::tuple<double,double,double> set3(s3->GetPoint(0).GetElement(0),s3->GetPoint(0).GetElement(1),s3->GetPoint(0).GetElement(2));
        std::tuple<double,double,double> set4(s4->GetPoint(0).GetElement(0),s4->GetPoint(0).GetElement(1),s4->GetPoint(0).GetElement(2));
        std::tuple<double,double,double> set5(s5->GetPoint(0).GetElement(0),s5->GetPoint(0).GetElement(1),s5->GetPoint(0).GetElement(2));
        std::tuple<double,double,double> set6(s6->GetPoint(0).GetElement(0),s6->GetPoint(0).GetElement(1),s6->GetPoint(0).GetElement(2));
        std::tuple<double,double,double> set7(s7->GetPoint(0).GetElement(0),s7->GetPoint(0).GetElement(1),s7->GetPoint(0).GetElement(2));
        std::tuple<double,double,double> set8(s8->GetPoint(0).GetElement(0),s8->GetPoint(0).GetElement(1),s8->GetPoint(0).GetElement(2));
        std::tuple<double,double,double> set9(s9->GetPoint(0).GetElement(0),s9->GetPoint(0).GetElement(1),s9->GetPoint(0).GetElement(2));
        QString path = "/home/or15/Desktop/Proj/JB/OrodNormals/n1/TrackingPoints";
        CemrgMeasure* measure = new CemrgMeasure();
        std::tuple<double,double,double> tra0 = measure->Deconvert(path,0).at(0);
        std::tuple<double,double,double> tra1 = measure->Deconvert(path,1).at(0);
        std::tuple<double,double,double> tra2 = measure->Deconvert(path,2).at(0);
        std::tuple<double,double,double> tra3 = measure->Deconvert(path,3).at(0);
        std::tuple<double,double,double> tra4 = measure->Deconvert(path,4).at(0);
        std::tuple<double,double,double> tra5 = measure->Deconvert(path,5).at(0);
        std::tuple<double,double,double> tra6 = measure->Deconvert(path,6).at(0);
        std::tuple<double,double,double> tra7 = measure->Deconvert(path,7).at(0);
        std::tuple<double,double,double> tra8 = measure->Deconvert(path,8).at(0);
        std::tuple<double,double,double> tra9 = measure->Deconvert(path,9).at(0);
        //qDebug() << measure->CalcDist3D(set0,tra0);
        //qDebug() << measure->CalcDist3D(set1,tra1);
        //qDebug() << measure->CalcDist3D(set2,tra2);
        //qDebug() << measure->CalcDist3D(set3,tra3);
        //qDebug() << measure->CalcDist3D(set4,tra4);
        //qDebug() << measure->CalcDist3D(set5,tra5);
        //qDebug() << measure->CalcDist3D(set6,tra6);
        //qDebug() << measure->CalcDist3D(set7,tra7);
        //qDebug() << measure->CalcDist3D(set8,tra8);
        //qDebug() << measure->CalcDist3D(set9,tra9);
    }
    /**
      * Manual segmentation test
      **/
    if (false) {
        //sanity check
        vtkSmartPointer<vtkCubeSource> mCT = vtkSmartPointer<vtkCubeSource>::New(); mCT->Update();
        vtkSmartPointer<vtkCubeSource> aCT = vtkSmartPointer<vtkCubeSource>::New(); aCT->Update();
        for (int i=0; i<4; i++) {
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
        mitk::Surface::Pointer mS = mST;//mitk::IOUtil::LoadSurface("/home/or15/Desktop/Proj/JB/OrodNormals/n1/Manuals/man-segmentation.vtk");
        mitk::Surface::Pointer aS = aST;//mitk::IOUtil::LoadSurface("/home/or15/Desktop/Proj/JB/OrodNormals/n1/Manuals/aut-segmentation.vtk");
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
        mitk::IOUtil::Save(rS,"/home/or15/Desktop/Proj/JB/OrodNormals/n1/Manuals/res-segmentation.vtk");
        qDebug() << "Filter Done!";

        //        double areaM = 0;
        //        CemrgStrains* str = new CemrgStrains();
        //        for (vtkIdType cellID = 0; cellID < mPD->GetNumberOfCells(); cellID++)
        //            areaM += str->GetCellArea(mPD, cellID);
        //        double areaR = 0;
        //        for (vtkIdType cellID = 0; cellID < rPD->GetNumberOfCells(); cellID++)
        //            areaR += str->GetCellArea(rPD, cellID);
        //        qDebug() << "Difference in %:" << (areaR*100.0)/areaM;
    }
}

void testMonoContract() {

    std::string path;
    float scales[10] = {1.00,0.95,0.90,0.85,0.80,0.80,0.85,0.90,0.95,1.00};

    //    path = "/home/or15/Desktop/Proj/JB/StrainTests/VTKs/MONOLV/transformed-0.vtk";
    //    mitk::Surface::Pointer surf = mitk::IOUtil::LoadSurface(path);
    //    vtkSmartPointer<vtkPolyData> pd = surf->GetVtkPolyData();
    //    double x_s = 0; double y_s = 0; double z_s = 0;
    //    for (int i=0; i<pd->GetNumberOfPoints(); i++) {
    //        x_s = x_s + pd->GetPoints()->GetPoint(i)[0];
    //        y_s = y_s + pd->GetPoints()->GetPoint(i)[1];
    //        z_s = z_s + pd->GetPoints()->GetPoint(i)[2];
    //    }//_for
    //    x_s /= pd->GetNumberOfPoints();
    //    y_s /= pd->GetNumberOfPoints();
    //    z_s /= pd->GetNumberOfPoints();
    //    double gravity[3] = {x_s, y_s, z_s};
    //    for (int i=0; i<pd->GetNumberOfPoints(); i++) {
    //        double* point = pd->GetPoint(i);
    //        point[0] = point[0] - gravity[0];
    //        point[1] = point[1] - gravity[1];
    //        point[2] = point[2] - gravity[2];
    //        pd->GetPoints()->SetPoint(i, point);
    //    }
    //    surf->SetVtkPolyData(pd);
    //    mitk::IOUtil::Save(surf, path);
    //    path = "/home/or15/Desktop/Proj/JB/StrainTests/VTKs/MONOLV/MONOpoints.mps";
    //    mitk::PointSet::Pointer points = mitk::IOUtil::LoadPointSet(path);
    //    for (int i=0; i<6; i++) {
    //        mitk::Point3D point;
    //        qDebug() << points->GetPoint(i).GetElement(2) << gravity[2];
    //        point.SetElement(0, (-1*points->GetPoint(i).GetElement(0) - gravity[0])*-1);
    //        point.SetElement(1, (-1*points->GetPoint(i).GetElement(1) - gravity[1])*-1);
    //        point.SetElement(2, points->GetPoint(i).GetElement(2) - gravity[2]);
    //        points->InsertPoint(i,point);
    //        qDebug() << points->GetPoint(i).GetElement(2) << gravity[2];
    //    }//for
    //    mitk::IOUtil::Save(points, path);


    for (int j=1; j<10; j++) {

        path = "/home/or15/Desktop/Proj/JB/StrainTests/VTKs/MONOLV/transformed-0.vtk";
        mitk::Surface::Pointer surf = mitk::IOUtil::LoadSurface(path);

        vtkSmartPointer<vtkPolyData> pd = surf->GetVtkPolyData();
        for (int i=0; i<pd->GetNumberOfPoints(); i++) {
            double* point = pd->GetPoint(i);
            point[0] = point[0] * scales[j];
            point[1] = point[1] * scales[j];
            point[2] = point[2] * scales[j];
            pd->GetPoints()->SetPoint(i, point);
        }
        surf->SetVtkPolyData(pd);

        path = "/home/or15/Desktop/Proj/JB/StrainTests/VTKs/MONOLV/transformed-" + std::to_string(j) + ".vtk";
        mitk::IOUtil::Save(surf, path);
    }
}

/**
 * Aux Tests for Strains
 * */
void auxTestsStrain() {

    /**
     * Test Strains
     **/
    if (false) {

        //        std::vector<mitk::Point3D> lm = ConvertMPS(lmNode);
        //        mitk::Surface::Pointer surf = ReadVTKMesh(meshNo); ZeroVTKMesh(lm.at(0), surf);
        //        mitk::Point3D centre = Circlefit3d(ZeroPoint(lm.at(0),lm.at(1)), ZeroPoint(lm.at(0),lm.at(2)), ZeroPoint(lm.at(0),lm.at(3)));
        //        mitk::Matrix<double,3,3> rotationMat = CalcRotationMatrix(centre, ZeroPoint(lm.at(0),lm.at(5))); RotateVTKMesh(rotationMat, surf);
        //        mitk::IOUtil::Save(surf, "/home/or15/Documents/MATLAB/" + std::to_string(meshNo) + ".vtk");

        //        std::vector<double> lamdas;
        //        if (meshNo!=0) {
        //            double val;
        //            std::ifstream infile("/home/or15/Documents/MATLAB/" + std::to_string(meshNo) + "-max.csv");
        //            while (infile >> val) {
        //                lamdas.push_back(val);
        //            }
        //        } else {
        //            for (vtkIdType cellID = 0; cellID < surf->GetVtkPolyData()->GetNumberOfCells(); cellID++)
        //                lamdas.push_back(0);
        //        }

        //        std::vector<double> strainRCL(16,0);
        //        for (vtkIdType cellID = 0; cellID < surf->GetVtkPolyData()->GetNumberOfCells(); cellID++) {
        //            //Ignore non AHA segments
        //            if (refCellLabels[cellID] == 0)
        //                continue;
        //            //Prepare plot values
        //            strainRCL.at(refCellLabels[cellID]-1) += lamdas.at(cellID);
        //        }//_for

        //        for (int i=0; i<16; i++)
        //            strainRCL.at(i) /= std::count(refCellLabels.begin(), refCellLabels.end(), i+1);
        //        return strainRCL;
    }

    /**
      TEST
      **/
    //qDebug() << "RCTR IS " << RCTR.GetElement(0) << RCTR.GetElement(1) << RCTR.GetElement(2);

    /**
     * Strain Test coordinates
     **/
    //    controller++;
    //    if(controller%100==0 && false) {
    //        //qDebug() << "Radial" << radiAxis.GetElement(0) << radiAxis.GetElement(1) << radiAxis.GetElement(2);
    //        //qDebug() << "Cross2" << Cross(circAxis,longAxis).GetElement(0) << Cross(circAxis,longAxis).GetElement(1) << Cross(circAxis,longAxis).GetElement(2);
    //        mitk::Surface::Pointer axisX = mitk::Surface::New();
    //        mitk::Surface::Pointer axisY = mitk::Surface::New();
    //        mitk::Surface::Pointer axisZ = mitk::Surface::New();
    //        vtkSmartPointer<vtkLineSource> line1 = vtkSmartPointer<vtkLineSource>::New();
    //        vtkSmartPointer<vtkLineSource> line2 = vtkSmartPointer<vtkLineSource>::New();
    //        vtkSmartPointer<vtkLineSource> line3 = vtkSmartPointer<vtkLineSource>::New();
    //        mitk::Point3D centre;
    //        double x = pt1[0] + pt2[0] + pt3[0]; centre.SetElement(0, x/3);
    //        double y = pt1[1] + pt2[1] + pt3[1]; centre.SetElement(1, y/3);
    //        double z = pt1[2] + pt2[2] + pt3[2]; centre.SetElement(2, z/3);

    //        line1->SetPoint1(centre.GetElement(0),centre.GetElement(1),centre.GetElement(2));
    //        line1->SetPoint2(centre.GetElement(0)+radiAxis.GetElement(0),
    //                         centre.GetElement(1)+radiAxis.GetElement(1),
    //                         centre.GetElement(2)+radiAxis.GetElement(2));
    //        line1->Update();
    //        axisX->SetVtkPolyData(line1->GetOutput());
    //        line2->SetPoint1(centre.GetElement(0),centre.GetElement(1),centre.GetElement(2));
    //        line2->SetPoint2(centre.GetElement(0)+longAxis.GetElement(0),
    //                         centre.GetElement(1)+longAxis.GetElement(1),
    //                         centre.GetElement(2)+longAxis.GetElement(2));
    //        line2->Update();
    //        axisY->SetVtkPolyData(line2->GetOutput());
    //        line3->SetPoint1(centre.GetElement(0),centre.GetElement(1),centre.GetElement(2));
    //        line3->SetPoint2(centre.GetElement(0)+circAxis.GetElement(0),
    //                         centre.GetElement(1)+circAxis.GetElement(1),
    //                         centre.GetElement(2)+circAxis.GetElement(2));
    //        line3->Update();
    //        axisZ->SetVtkPolyData(line3->GetOutput());

    //        temporaries.push_back(axisX);
    //        temporaries.push_back(axisY);
    //        temporaries.push_back(axisZ);}

    //    /**
    //     * STRAIN TEST
    //     **/
    //    if (false) {
    //        for (unsigned int i=0; i<strain->temporaries.size(); i++) {
    //            mitk::DataNode::Pointer aNode = mitk::DataNode::New();
    //            aNode->SetData(strain->temporaries.at(i));
    //            this->GetDataStorage()->Add(aNode, node);
    //        }
    //    }
}
