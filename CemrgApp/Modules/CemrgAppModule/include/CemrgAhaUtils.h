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
 * Power Transmitter Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef CemrgAhaUtils_h
#define CemrgAhaUtils_h

// Qmitk
#include <mitkImage.h>
#include <mitkPointSet.h>
#include <MitkCemrgAppModuleExports.h>

// VTK
//#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
//#include <vtkMatrix3x3.h>

// #include <MyCemrgLibExports.h>


class MITKCEMRGAPPMODULE_EXPORT CemrgAhaUtils {

public:

    CemrgAhaUtils();
    mitk::Surface::Pointer ReferenceAHA(mitk::DataNode::Pointer lmNode, mitk::Surface::Pointer refSurface);

private:

    // Copied from Strain for AHA mapping
    mitk::Matrix<double, 3, 3> CalcRotationMatrix(mitk::Point3D point1, mitk::Point3D point2);
    mitk::Point3D Circlefit3d(mitk::Point3D point1, mitk::Point3D point2, mitk::Point3D point3);
    mitk::Point3D RotatePoint(mitk::Matrix<double, 3, 3> rotationMatrix, mitk::Point3D point);
    void RotateVTKMesh(mitk::Matrix<double, 3, 3> rotationMatrix, mitk::Surface::Pointer surface);
    mitk::Point3D ZeroPoint(mitk::Point3D apex, mitk::Point3D point);
    void ZeroVTKMesh(mitk::Point3D apex, mitk::Surface::Pointer surface);
};

#endif // CemrgAhaUtils_h
