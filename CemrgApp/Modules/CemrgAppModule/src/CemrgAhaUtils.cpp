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
 * Power transmitter Calculations Tools for MITK
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
 * angela.lee@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// Qmitk
#include <mitkIOUtil.h>
#include <mitkProperties.h>

// Qt
#include <QMessageBox>
#include <QDebug>

// C++ Standard
#include <vector>

// VTK
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>


#include "CemrgAhaUtilsUtils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

CemrgAhaUtils::CemrgAhaUtils() {
}

mitk::Surface::Pointer CemrgAhaUtils::ReferenceAHA(mitk::DataNode::Pointer lmNode, mitk::Surface::Pointer refSurface) {

    //Read the mesh data
    vtkSmartPointer<vtkPolyData> pd = refSurface->GetVtkPolyData();

    //Prepare landmarks
    std::vector<mitk::Point3D> LandMarks = ConvertMPS(lmNode);
    // Siemens 4 LM = [apex, basecenter, RV1, RV2]; Siemens 7 LM = [apex, baseMV1, baseMV2, baseMV3, RV1, RV2, apex (just to make it 7 and different from manual LM)]
    // Manual 6 LM = [apex, MV1, MV2, MV3, RV1, RV2]
    if (LandMarks.size() != 4 && LandMarks.size() != 6 && LandMarks.size() != 7) {
        QMessageBox::warning(NULL, "Attention", "Error creating path:\n" + aPath);


        return refSurface;
    }
    mitk::Point3D centre, RIV1, RIV2, APEX, MIV1, MIV2, MIV3;
    if (LandMarks.size() >= 6) {
        APEX = LandMarks.at(0);
        MIV1 = LandMarks.at(1);
        MIV2 = LandMarks.at(2);
        MIV3 = LandMarks.at(3);
        RIV1 = LandMarks.at(4);
        RIV2 = LandMarks.at(5);
        //Calcaulte a circle through the mitral valve points
        centre = Circlefit3d(ZeroPoint(APEX, MIV1), ZeroPoint(APEX, MIV2), ZeroPoint(APEX, MIV3));
        // Zero all points relative to apex
    } else if (LandMarks.size() == 4) {
        APEX = LandMarks.at(0);
        RIV1 = LandMarks.at(2);
        RIV2 = LandMarks.at(3);
        centre = ZeroPoint(APEX, LandMarks.at(1));
    }

    //Zero all points relative to apex
    RIV1 = ZeroPoint(APEX, RIV1);
    RIV2 = ZeroPoint(APEX, RIV2);
    ZeroVTKMesh(APEX, refSurface);
    APEX = ZeroPoint(APEX, APEX);

    //Calculate a circle through the mitral valve points
    mitk::Point3D RCTR;

    //Define Rotation matrix
    mitk::Matrix<double, 3, 3> rotationMat = CalcRotationMatrix(centre, RIV2);

    //Rotate mesh to new frame
    RotateVTKMesh(rotationMat, refSurface);
    //Angle RV cusp 2
    double RVangle1 = atan2(RIV1.GetElement(1), RIV1.GetElement(0));
    double RVangle2 = atan2(RIV2.GetElement(1), RIV2.GetElement(0));
    double appendAngle;
    //Assuming RV angle1 < RV angle 2
    if (RVangle1 < RVangle2) {
        RVangle1 = atan2(RIV2.GetElement(1), RIV2.GetElement(0));
        RVangle2 = atan2(RIV1.GetElement(1), RIV1.GetElement(0));
    }
    MITK_INFO << ("RVangles: (1) " + QString::number(RVangle1) + ", (2) " + QString::number(RVangle2)).toStdString();

    if ((LandMarks.size() == 6)) { // only do this for manual segmentation, with 6 points
        appendAngle = -RVangle1;
        //appendAngle -= freeA - ( M_PI / 3 );
    } else {
        if (RVangle1 > 0) {
            appendAngle = -((M_PI - RVangle1) / 2 + RVangle1) + M_PI / 3;
        } else {
            appendAngle = -RVangle1 / 2 + M_PI / 3;
        }
    }
    qDebug() << "appendAngle " << appendAngle;
    //Point angles
    std::vector<double> pAngles;
    for (int i = 0; i < pd->GetNumberOfPoints(); i++) {
        double pAngle;
        double* pt = pd->GetPoint(i);
        pAngle = atan2(pt[1], pt[0]) + appendAngle;
        pAngle = pAngle * (pAngle > 0 ? 1 : 0) + (2 * M_PI + pAngle) * (pAngle < 0 ? 1 : 0);
        pAngles.push_back(pAngle);
    }
    //Setup flattened AHA mesh
    mitk::Surface::Pointer flatSurface = refSurface->Clone();
    vtkSmartPointer<vtkPolyData> poly = flatSurface->GetVtkPolyData();
    for (int i = 0; i < poly->GetNumberOfPoints(); i++) {
        double* point = poly->GetPoint(i);
        double  radii = point[2];
        double  theta = pAngles.at(i);
        point[0] = radii * cos(theta);
        point[1] = radii * sin(theta);
        point[2] = 0;
        poly->GetPoints()->SetPoint(i, point);
    }//_for
    poly->BuildLinks();

    //Return to the original position and rotation
    APEX = ZeroPoint(LandMarks.at(0), APEX);
    RotateVTKMesh(rotationMat.GetInverse().as_matrix(), refSurface);
    ZeroVTKMesh(APEX, refSurface);
    return refSurface;
}

/**************************************************************************************************
 *************** PRIVATE FUNCTIONS ****************************************************************
 **************************************************************************************************/


std::vector<mitk::Point3D> CemrgAhaUtils::ConvertMPS(mitk::DataNode::Pointer node) {

    std::vector<mitk::Point3D> points;
    mitk::BaseData::Pointer data = node->GetData();
    mitk::PointSet::Pointer set = dynamic_cast<mitk::PointSet*>(data.GetPointer());

    if (set.IsNull())
        return points;
    for (mitk::PointSet::PointsIterator it = set->Begin(); it != set->End(); ++it) {
        mitk::Point3D point;
        point.SetElement(0, it.Value().GetElement(0));
        point.SetElement(1, it.Value().GetElement(1));
        point.SetElement(2, it.Value().GetElement(2));
        points.push_back(point);
    }//for

    return points;
}

/*************************************************************************************************
 * Copied from CemrgStrain
 *************************************************************************************************/
mitk::Matrix<double, 3, 3> CemrgAhaUtils::CalcRotationMatrix(mitk::Point3D point1, mitk::Point3D point2) {

    //X Axis
    mitk::Matrix<double, 1, 3> vec;
    vec[0][0] = point1.GetElement(0);
    vec[0][1] = point1.GetElement(1);
    vec[0][2] = point1.GetElement(2);
    double theta_x = atan(vec[0][1] / vec[0][2]);

    mitk::Matrix<double, 3, 3> R_x;
    R_x[0][0] = 1;
    R_x[0][1] = 0;
    R_x[0][2] = 0;
    R_x[1][0] = 0;
    R_x[1][1] = cos(theta_x);
    R_x[1][2] = -sin(theta_x);
    R_x[2][0] = 0;
    R_x[2][1] = sin(theta_x);
    R_x[2][2] = cos(theta_x);

    //Y Axis
    mitk::Matrix<double, 3, 1> vecX;
    vecX = R_x * vec.GetTranspose();
    double theta_y = atan(-vecX[0][0] / vecX[2][0]);

    mitk::Matrix<double, 3, 3> R_y;
    R_y[0][0] = cos(theta_y);
    R_y[0][1] = 0;
    R_y[0][2] = sin(theta_y);
    R_y[1][0] = 0;
    R_y[1][1] = 1;
    R_y[1][2] = 0;
    R_y[2][0] = -sin(theta_y);
    R_y[2][1] = 0;
    R_y[2][2] = cos(theta_y);

    //Z Axis
    point2 = RotatePoint(R_y * R_x, point2);
    double theta_z = atan(-point2.GetElement(1) / point2.GetElement(0));

    mitk::Matrix<double, 3, 3> R_z;
    R_z[0][0] = cos(theta_z);
    R_z[0][1] = -sin(theta_z);
    R_z[0][2] = 0;
    R_z[1][0] = sin(theta_z);
    R_z[1][1] = cos(theta_z);
    R_z[1][2] = 0;
    R_z[2][0] = 0;
    R_z[2][1] = 0;
    R_z[2][2] = 1;

    //Rotation Matrix
    mitk::Matrix<double, 3, 3> R;
    R = R_z * R_y * R_x;
    return R;
}

mitk::Point3D CemrgAhaUtils::Circlefit3d(mitk::Point3D point1, mitk::Point3D point2, mitk::Point3D point3) {

    //v1, v2 describe the vectors from p1 to p2 and p3, resp.
    mitk::Point3D v1;
    mitk::Point3D v2;
    for (int i = 0; i < 3; i++) {
        v1.SetElement(i, point2.GetElement(i) - point1.GetElement(i));
        v2.SetElement(i, point3.GetElement(i) - point1.GetElement(i));
    }

    //l1, l2 describe the lengths of those vectors
    double l1 = sqrt(pow(v1.GetElement(0), 2) + pow(v1.GetElement(1), 2) + pow(v1.GetElement(2), 2));
    double l2 = sqrt(pow(v2.GetElement(0), 2) + pow(v2.GetElement(1), 2) + pow(v2.GetElement(2), 2));

    //v1n, v2n describe the normalized vectors v1 and v2
    mitk::Point3D v1n = v1;
    mitk::Point3D v2n = v2;
    for (int i = 0; i < 3; i++) {
        v1n.SetElement(i, v1n.GetElement(i) / l1);
        v2n.SetElement(i, v2n.GetElement(i) / l2);
    }

    //nv describes the normal vector on the plane of the circle
    mitk::Point3D nv;
    nv.SetElement(0, v1n.GetElement(1) * v2n.GetElement(2) - v1n.GetElement(2) * v2n.GetElement(1));
    nv.SetElement(1, v1n.GetElement(2) * v2n.GetElement(0) - v1n.GetElement(0) * v2n.GetElement(2));
    nv.SetElement(2, v1n.GetElement(0) * v2n.GetElement(1) - v1n.GetElement(1) * v2n.GetElement(0));

    //v2nb: orthogonalization of v2n against v1n
    double dotp = v2n.GetElement(0) * v1n.GetElement(0) +
        v2n.GetElement(1) * v1n.GetElement(1) +
        v2n.GetElement(2) * v1n.GetElement(2);

    mitk::Point3D v2nb = v2n;
    for (int i = 0; i < 3; i++) {
        v2nb.SetElement(i, v2nb.GetElement(i) - dotp * v1n.GetElement(i));
    }

    //Normalize v2nb
    double l2nb = sqrt(pow(v2nb.GetElement(0), 2) + pow(v2nb.GetElement(1), 2) + pow(v2nb.GetElement(2), 2));
    for (int i = 0; i < 3; i++) {
        v2nb.SetElement(i, v2nb.GetElement(i) / l2nb);
    }

    //Calculate 2d coordinates of points in each plane
    mitk::Point2D p3_2d;
    for (int i = 0; i < 3; i++) {
        p3_2d.SetElement(0, p3_2d.GetElement(0) + v2.GetElement(i) * v1n.GetElement(i));
        p3_2d.SetElement(1, p3_2d.GetElement(1) + v2.GetElement(i) * v2nb.GetElement(i));
    }

    //Calculate the fitting circle
    double a = l1;
    double b = p3_2d.GetElement(0);
    double c = p3_2d.GetElement(1);
    double t = .5 * (a - b) / c;
    double scale1 = b / 2 + c * t;
    double scale2 = c / 2 - b * t;

    //centre
    mitk::Point3D centre;
    for (int i = 0; i < 3; i++) {
        double val = point1.GetElement(i) + (scale1 * v1n.GetElement(i)) + (scale2 * v2nb.GetElement(i));
        centre.SetElement(i, val);
    }
    // qDebug() << centre.GetElement(0) << centre.GetElement(1) << centre.GetElement(2);

    /*radius
    double radius = sqrt(pow(centre.GetElement(0) - point1.GetElement(0),2) +
                         pow(centre.GetElement(1) - point1.GetElement(1),2) +
                         pow(centre.GetElement(2) - point1.GetElement(2),2));*/
    return centre;
}


mitk::Point3D CemrgAhaUtils::RotatePoint(mitk::Matrix<double, 3, 3> rotationMatrix, mitk::Point3D point) {

    mitk::Matrix<double, 1, 3> vec;
    mitk::Matrix<double, 3, 1> ans;

    vec[0][0] = point.GetElement(0);
    vec[0][1] = point.GetElement(1);
    vec[0][2] = point.GetElement(2);
    ans = rotationMatrix * vec.GetTranspose();
    point.SetElement(0, ans[0][0]);
    point.SetElement(1, ans[1][0]);
    point.SetElement(2, ans[2][0]);

    return point;
}

void CemrgAhaUtils::RotateVTKMesh(mitk::Matrix<double, 3, 3> rotationMatrix, mitk::Surface::Pointer surface) {

    //Retrieve the data
    vtkSmartPointer<vtkPolyData> pd = surface->GetVtkPolyData();

    //Rotate mesh into new coordiantes
    for (int i = 0; i < pd->GetNumberOfPoints(); i++) {
        mitk::Point3D point;
        double* pt = pd->GetPoint(i);
        point.SetElement(0, pt[0]);
        point.SetElement(1, pt[1]);
        point.SetElement(2, pt[2]);
        point = RotatePoint(rotationMatrix, point);
        pt[0] = point.GetElement(0);
        pt[1] = point.GetElement(1);
        pt[2] = point.GetElement(2);
        pd->GetPoints()->SetPoint(i, pt);
    }
}

mitk::Point3D CemrgAhaUtils::ZeroPoint(mitk::Point3D apex, mitk::Point3D point) {

    //Zero relative to the apex
    point.SetElement(0, point.GetElement(0) - apex.GetElement(0));
    point.SetElement(1, point.GetElement(1) - apex.GetElement(1));
    point.SetElement(2, point.GetElement(2) - apex.GetElement(2));
    return point;
}

void CemrgAhaUtils::ZeroVTKMesh(mitk::Point3D apex, mitk::Surface::Pointer surface) {

    //Retrieve the data
    vtkSmartPointer<vtkPolyData> pd = surface->GetVtkPolyData();
    for (int i = 0; i < pd->GetNumberOfPoints(); i++) {
        double* point = pd->GetPoint(i);
        point[0] = point[0] - apex.GetElement(0);
        point[1] = point[1] - apex.GetElement(1);
        point[2] = point[2] - apex.GetElement(2);
        pd->GetPoints()->SetPoint(i, point);
    }
}
