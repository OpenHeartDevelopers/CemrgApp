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
 * Power Transmitter Calculations Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * angela.lee@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef CemrgPower_h
#define CemrgPower_h

// Qmitk
#include <mitkImage.h>
#include <mitkPointSet.h>
#include <MitkCemrgAppModuleExports.h>

// VTK
#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include <vtkMatrix3x3.h>

class MITKCEMRGAPPMODULE_EXPORT CemrgPower {

public:

    CemrgPower();
    CemrgPower(QString dir, int ribSpacing);

    mitk::Surface::Pointer MapPowerTransmitterToLandmarks(mitk::DataNode::Pointer lmNode);
    mitk::Surface::Pointer CalculateAcousticIntensity(mitk::Surface::Pointer endoMesh);
    mitk::Surface::Pointer ReferenceAHA(mitk::PointSet::Pointer lmNode, mitk::Surface::Pointer refSurface);

private:

    QString projectDirectory;
    int currentRibSpacing;

    double sinc(const double x);
    void normalise(double v[]);
    void crossProduct(double a[], double b[], double product[]);
    double dotProduct(double a[], double b[]);
    std::vector<mitk::Point3D> ConvertMPS(mitk::DataNode::Pointer node);
    void fcn_RotationToUnity(const double v[], vtkSmartPointer<vtkMatrix3x3>& RotationMatrix);
    void fcn_RotationFromTwoVectors(double a[], double b[], vtkSmartPointer<vtkMatrix3x3>& RotationMatrix);

    // Copied from Strain for AHA mapping
    mitk::Matrix<double, 3, 3> CalcRotationMatrix(mitk::Point3D point1, mitk::Point3D point2);
    mitk::Point3D Circlefit3d(mitk::Point3D point1, mitk::Point3D point2, mitk::Point3D point3);
    mitk::Point3D RotatePoint(mitk::Matrix<double, 3, 3> rotationMatrix, mitk::Point3D point);
    void RotateVTKMesh(mitk::Matrix<double, 3, 3> rotationMatrix, mitk::Surface::Pointer surface);
    mitk::Point3D ZeroPoint(mitk::Point3D apex, mitk::Point3D point);
    void ZeroVTKMesh(mitk::Point3D apex, mitk::Surface::Pointer surface);
};

#endif // CemrgPower_h
