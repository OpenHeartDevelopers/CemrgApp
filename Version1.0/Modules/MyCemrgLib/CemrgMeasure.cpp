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
 * Measurements Tools for MITK
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

//Vtk
#include <vtkPolyData.h>
#include <vtkMath.h>
#include <vtkIdList.h>
#include <vtkMassProperties.h>

//Qmitk
#include <mitkIOUtil.h>

//Qt
#include "CemrgMeasure.h"


void CemrgMeasure::Convert(QString dir, mitk::DataNode::Pointer node) {

    mitk::BaseData::Pointer data = node->GetData();
    mitk::PointSet::Pointer set = dynamic_cast<mitk::PointSet*>(data.GetPointer());

    int items = 0;
    for (mitk::PointSet::PointsIterator it = set->Begin(); it != set->End(); ++it) {
        items++;
    }

    ofstream file;
    double x, y, z;
    file.open(dir.toStdString() + mitk::IOUtil::GetDirectorySeparator() + "input.vtk");

    //Header
    file << "# vtk DataFile Version 3.0" << endl;
    file << "vtk output" << endl;
    file << "ASCII" << endl;
    file << "DATASET POLYDATA" << endl;
    file << "POINTS " << items << " float" << endl;

    //Body
    items = 0;
    for (mitk::PointSet::PointsIterator it = set->Begin(); it != set->End(); ++it) {
        x = it.Value().GetElement(0) * -1;
        y = it.Value().GetElement(1) * -1;
        z = it.Value().GetElement(2);
        items+=3;

        file << x << " " << y << " " << z << " ";
        if (items%9==0)
            file << endl;
    }//for
    file.close();
}

std::vector <std::tuple<double, double, double>> CemrgMeasure::Deconvert(QString dir, int noFile) {

    double x, y, z;
    std::string line;
    unsigned int items;
    std::vector<std::string> tokens;
    std::vector <std::tuple<double, double, double>> points;
    ifstream file(dir.toStdString() + mitk::IOUtil::GetDirectorySeparator() + "transformed-" + std::to_string(noFile) + ".vtk");

    if (file.is_open()) {

        //Skip header
        for (int i=0; i<5; i++)
            getline(file,line);

        //Read main
        while (getline(file,line)) {
            items = 0;
            tokens = Split(line,tokens);
            while (items<tokens.size()) {
                x = std::stod(tokens.at(items + 0)) * -1;
                y = std::stod(tokens.at(items + 1)) * -1;
                z = std::stod(tokens.at(items + 2));
                items+=3;
                points.push_back(std::tuple<double, double, double>(x,y,z));
            }
        }
        file.close();
    }//_if

    return points;
}

double CemrgMeasure::CalcDistance(std::vector <std::tuple<double, double, double>>& points) {

    double dist = 0;

    if (points.size() != 2)
        dist = -1;
    else
        dist = CalcDist3D(points.at(0), points.at(1));

    return dist;
}

double CemrgMeasure::CalcPerimeter(std::vector <std::tuple<double, double, double>>& points) {

    double peri = 0;

    if (points.size() < 3) {
        peri = -1;
    } else {
        for (unsigned int i=0; i<points.size(); i++) {
            if (i < points.size() - 1)
                peri = peri + CalcDist3D(points.at(i), points.at(i+1));
            else
                peri = peri + 0; //CalcDist3D(points.at(i), points.at(0)); closed chain perimeter
        }//_for
    }

    return peri;
}

double CemrgMeasure::CalcArea(std::vector<std::tuple<double, double, double>>& points) {

    double area = 0;
    std::tuple<double, double, double> mean = CalcMean(points);

    if (points.size() < 3) {
        area = -1;
    } else {
        for(unsigned int i=0; i<points.size(); i++) {
            if (i < points.size() - 1)
                area = area + Heron(points.at(i), points.at(i+1), mean);
            else
                area = area + Heron(points.at(i), points.at(0), mean);
        }//_for
    }

    return area;
}

mitk::Point3D CemrgMeasure::FindCentre(mitk::PointSet::Pointer pointset) {

    std::vector <std::tuple<double, double, double>> points;
    for (int i=0; i<pointset->GetSize(); i++) {

        std::tuple<double, double, double> point;
        std::get<0>(point) = pointset->GetPoint(i).GetElement(0);
        std::get<1>(point) = pointset->GetPoint(i).GetElement(1);
        std::get<2>(point) = pointset->GetPoint(i).GetElement(2);
        points.push_back(point);
    }//_for

    mitk::Point3D centrePoint;
    std::tuple<double, double, double> centre = CalcMean(points);
    centrePoint[0] = std::get<0>(centre);
    centrePoint[1] = std::get<1>(centre);
    centrePoint[2] = std::get<2>(centre);

    return centrePoint;
}

double CemrgMeasure::GetSphericity(vtkPolyData* LAC_poly) {

    double LACA;
    double LAC_mc[3];
    double AR, sigma, Sphericity;

    //containers for centre of mass, area, etc.
    double** TiMC;
    double* TiA;

    TiMC = new double*[LAC_poly->GetNumberOfCells()];
    for (int i=0; i<LAC_poly->GetNumberOfCells(); i++) {
        TiMC[i] = new double[3];
    }

    TiA = new double[LAC_poly->GetNumberOfCells()];

    GetCentreOfMassOfEachT(LAC_poly, TiMC);
    GetArea(LAC_poly, TiA, LACA);
    GetCentreOfMass(TiMC, TiA, LAC_poly->GetNumberOfCells(), LACA, LAC_mc);
    GetAverageRadius(TiMC, TiA, LAC_poly->GetNumberOfCells(), LAC_mc,  LACA, AR);
    LASphericity(TiMC, TiA, LAC_poly->GetNumberOfCells(), LAC_mc, LACA, AR, sigma, Sphericity);

    return Sphericity;
}

double CemrgMeasure::calcVolumeMesh(mitk::Surface::Pointer surface) {

    vtkSmartPointer<vtkMassProperties> mass = vtkSmartPointer<vtkMassProperties>::New();
    mass->SetInputData(surface->GetVtkPolyData());
    mass->Update();
    return mass->GetVolume();
}

double CemrgMeasure::calcSurfaceMesh(mitk::Surface::Pointer surface) {

    vtkSmartPointer<vtkMassProperties> mass = vtkSmartPointer<vtkMassProperties>::New();
    mass->SetInputData(surface->GetVtkPolyData());
    mass->Update();
    return mass->GetSurfaceArea();
}

/********************************************
 *        Private Members Defintions        *
 ********************************************/

std::tuple<double, double, double> CemrgMeasure::CalcMean(std::vector <std::tuple<double, double, double>>& points) {

    double x_s = 0;
    double y_s = 0;
    double z_s = 0;

    for(unsigned int i=0; i<points.size(); i++) {
        x_s = x_s + std::get<0>(points.at(i));
        y_s = y_s + std::get<1>(points.at(i));
        z_s = z_s + std::get<2>(points.at(i));
    }//_for

    x_s /= points.size();
    y_s /= points.size();
    z_s /= points.size();
    //Mean point
    std::tuple<double, double, double> mean(x_s, y_s, z_s);

    return mean;
}

double CemrgMeasure::CalcDist3D(
        std::tuple<double, double, double>& pointA,
        std::tuple<double, double, double>& pointB) {

    double x_d = std::get<0>(pointA) - std::get<0>(pointB);
    double y_d = std::get<1>(pointA) - std::get<1>(pointB);
    double z_d = std::get<2>(pointA) - std::get<2>(pointB);

    //Distance between two points
    double distance = sqrt(pow(x_d,2) + pow(y_d,2) + pow(z_d,2));

    return distance;
}

double CemrgMeasure::Heron(
        std::tuple<double, double, double>& pointA,
        std::tuple<double, double, double>& pointB,
        std::tuple<double, double, double>& centre) {

    //3D distances
    double ab = CalcDist3D(pointA, pointB);
    double ac = CalcDist3D(pointA, centre);
    double bc = CalcDist3D(pointB, centre);

    //Calculate area of triangle
    double speri = (ab+ac+bc) / 2;
    double tArea = sqrt(speri*(speri-ab)*(speri-ac)*(speri-bc));

    return tArea;
}

std::vector<std::string>& CemrgMeasure::Split(const std::string& str, std::vector<std::string>& elements) {

    elements.clear();
    std::string item;
    std::stringstream ss(str);
    while (getline(ss, item, ' '))
        elements.push_back(item);
    return elements;
}

void CemrgMeasure::GetArea(vtkPolyData* polys, double* TiA, double& LACA) {

    double p1[3], p2[3], p3[3];
    double p1_p2[3];		// p1 minus p2
    double p2_p3[3];		// p2 minus p3
    double crossProduct[3];
    LACA = 0;
    vtkSmartPointer<vtkIdList> cell_points = vtkSmartPointer<vtkIdList>::New();

    for (int i=0; i<polys->GetNumberOfCells(); i++) {

        //vtkIdType neighbor_point;
        polys->GetCellPoints(i, cell_points);

        // Get all three vertices
        vtkIdType neighbor_point_id = cell_points->GetId(0);
        polys->GetPoint( neighbor_point_id, p1);

        neighbor_point_id = cell_points->GetId(1);
        polys->GetPoint( neighbor_point_id, p2);

        neighbor_point_id = cell_points->GetId(2);
        polys->GetPoint( neighbor_point_id, p3);

        vtkMath::Subtract(p1, p2, p1_p2);
        vtkMath::Subtract(p2, p3, p2_p3);

        vtkMath::Cross(p1_p2, p2_p3, crossProduct);
        TiA[i] = 0.5*vtkMath::Norm(crossProduct);
        LACA += TiA[i];
    }
    cout << "Total area LACA = " << LACA << endl;
}

void CemrgMeasure::GetCentreOfMass(double** TiMC, double* TiA, int TiMC_size, double LACA, double* LAC_mc) {

    LAC_mc[0] = 0; LAC_mc[1] = 0; LAC_mc[2] = 0;
    double contrib;

    for (int i=0; i<TiMC_size; i++) {

        contrib = TiA[i]/LACA;
        LAC_mc[0]+=(contrib)*TiMC[i][0];
        LAC_mc[1]+=(contrib)*TiMC[i][1];
        LAC_mc[2]+=(contrib)*TiMC[i][2];
    }
}

void CemrgMeasure::GetCentreOfMassOfEachT(vtkPolyData* polys, double** TiMC) {

    int num_points = 0;
    vtkIdType num_cell_points;
    double cX, cY, cZ;
    double cP[3];
    vtkSmartPointer<vtkIdList> cell_points = vtkSmartPointer<vtkIdList>::New();

    for (int i=0; i<polys->GetNumberOfCells(); i++) {

        cX=0; cY=0; cZ=0; num_points=0;
        vtkIdType neighbor_point;
        polys->GetCellPoints(i, cell_points);
        num_cell_points = cell_points->GetNumberOfIds();

        for (neighbor_point = 0; neighbor_point < num_cell_points; ++neighbor_point) {

            // Get the neighbour point id
            vtkIdType neighbor_point_id = cell_points->GetId(neighbor_point);
            // Get the neighbour point position
            polys->GetPoint( neighbor_point_id, cP );
            cX+=cP[0]; cY+=cP[1]; cZ+=cP[2];
            num_points++;
        }
        cX/=num_points;
        cY/=num_points;
        cZ/=num_points;
        TiMC[i][0] = cX;
        TiMC[i][1] = cY;
        TiMC[i][2] = cZ;
    }
}

void CemrgMeasure::GetAverageRadius(double** TiMC, double* TiA, double TiMC_size, double* LAC_mc, double LACA, double& AR) {

    double contrib;
    double a[3], b[3], a_minus_b[3];
    AR = 0;

    for (int i=0; i<TiMC_size; i++) {

        a[0] = LAC_mc[0]; a[1] = LAC_mc[1]; a[2] = LAC_mc[2];
        b[0] = TiMC[i][0]; b[1] = TiMC[i][1]; b[2] = TiMC[i][2];
        vtkMath::Subtract(a,b,a_minus_b);
        contrib = TiA[i]/LACA;
        AR += contrib*vtkMath::Norm(a_minus_b);
    }
}

void CemrgMeasure::LASphericity(
        double** TiMC, double* TiA, double TiMC_size, double* LAC_mc, double LACA, double AR, double& sigma, double& Sphericity) {

    double contrib;
    double a[3], b[3], a_minus_b[3];
    double phi=0, eps=0, CVS=0;

    for (int i=0; i<TiMC_size; i++) {

        contrib = TiA[i]/LACA;
        a[0] = LAC_mc[0]; a[1] = LAC_mc[1]; a[2] = LAC_mc[2];
        b[0] = TiMC[i][0]; b[1] = TiMC[i][1]; b[2] = TiMC[i][2];
        vtkMath::Subtract(a,b,a_minus_b);
        phi = (vtkMath::Norm(a_minus_b) - AR)*(vtkMath::Norm(a_minus_b) - AR);
        eps += contrib*phi;
    }

    sigma = sqrt(eps);
    CVS = sigma/AR;
    Sphericity = 100*(1-CVS);
}
