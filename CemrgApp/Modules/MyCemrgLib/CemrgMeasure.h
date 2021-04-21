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

#ifndef CemrgMeasure_h
#define CemrgMeasure_h

#include <MyCemrgLibExports.h>
#include <QDebug>


class MyCemrgLib_EXPORT CemrgMeasure {

public:

    //Point to Point Tools
    void Convert(QString dir, mitk::DataNode::Pointer);
    std::vector <std::tuple<double, double, double>> Deconvert(QString dir, int noFile);
    double CalcDistance(std::vector <std::tuple<double, double, double>>& points);
    double CalcPerimeter(std::vector <std::tuple<double, double, double>>& points);
    double CalcArea(std::vector <std::tuple<double, double, double>>& points);
    mitk::Point3D FindCentre(mitk::PointSet::Pointer pointset);

    //Sphericity Tools
    double GetSphericity(vtkPolyData* poly);

    //Mesh Mass Tools
    double calcVolumeMesh(mitk::Surface::Pointer surface);
    double calcSurfaceMesh(mitk::Surface::Pointer surface);

private:

    //Point to Point Tools
    std::tuple<double, double, double> CalcMean(std::vector <std::tuple<double, double, double>>& points);
    double CalcDist3D(std::tuple<double, double, double>& pointA, std::tuple<double, double, double>& pointB);
    double Heron(std::tuple<double, double, double>& pointA, std::tuple<double, double, double>& pointB, std::tuple<double, double, double>& centre);
    std::vector<std::string>& Split(const std::string& str, std::vector<std::string>& elements);

    //Sphericity Tools
    void GetArea(vtkPolyData* polys, double* TiA, double& LACA);
    void GetCentreOfMass(double** TiMC, double* TiA, int TiMC_size, double LACA, double* LAC_mc);
    void GetCentreOfMassOfEachT(vtkPolyData* polys, double** TiMC);
    void GetAverageRadius(double** TiMC, double* TiA, double TiMC_size, double* LAC_mc, double LACA, double& AR);
    void LASphericity(double** TiMC, double* TiA, double TiMC_size, double* LAC_mc, double LACA, double AR, double& sigma, double& Sphericity);
};

#endif // CemrgMeasure_h
