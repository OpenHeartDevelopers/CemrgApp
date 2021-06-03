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
 * Measurements Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef CemrgMeasure_h
#define CemrgMeasure_h

#include <MitkCemrgAppModuleExports.h>
#include <mitkPoint.h>
#include <mitkPointSet.h>
#include <mitkDataNode.h>
#include <mitkSurface.h>
#include <vtkPolyData.h>
#include <QDebug>

class MITKCEMRGAPPMODULE_EXPORT CemrgMeasure {

public:
    typedef std::tuple<double, double, double> Point;
    typedef std::vector<Point> Points;

    //Point to Point Tools
    void Convert(QString dir, mitk::DataNode::Pointer);
    Points Deconvert(QString dir, int noFile);
    double CalcDistance(Points& points);
    double CalcPerimeter(Points& points);
    double CalcArea(Points& points);
    mitk::Point3D FindCentre(mitk::PointSet::Pointer pointset);

    //Sphericity Tools
    double GetSphericity(vtkPolyData* poly);

    //Mesh Mass Tools
    double calcVolumeMesh(mitk::Surface::Pointer surface);
    double calcSurfaceMesh(mitk::Surface::Pointer surface);

private:

    //Point to Point Tools
    Point CalcMean(Points& points);
    double CalcDist3D(Point& pointA, Point& pointB);
    double Heron(Point& pointA, Point& pointB, Point& centre);
    std::vector<std::string>& Split(const std::string& str, std::vector<std::string>& elements);

    //Sphericity Tools
    void GetArea(vtkPolyData* polys, double* TiA, double& LACA);
    void GetCentreOfMass(double** TiMC, double* TiA, int TiMC_size, double LACA, double* LAC_mc);
    void GetCentreOfMassOfEachT(vtkPolyData* polys, double** TiMC);
    void GetAverageRadius(double** TiMC, double* TiA, double TiMC_size, double* LAC_mc, double LACA, double& AR);
    void LASphericity(double** TiMC, double* TiA, double TiMC_size, double* LAC_mc, double LACA, double AR, double& sigma, double& Sphericity);
};

#endif // CemrgMeasure_h
