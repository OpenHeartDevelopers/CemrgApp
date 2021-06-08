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
 * CEMRGAPPMODULE TESTS
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// C++ Standard
#include <vector>
#include <tuple>
#include <string>
#include <array>
#include <utility>
#include <memory>

// Qt
#include <QtTest/QtTest>
#include <QFileInfo>

// MITK
#include <mitkIOUtil.h>
#include <mitkTestingConfig.h>

// VTK
#include <vtkPolyDataNormals.h>

// CemrgApp
#include <CemrgMeasure.h>

using namespace std;

class TestCemrgMeasure : public QObject {

    Q_OBJECT

private:
    unique_ptr<CemrgMeasure> cemrgMeasure { new CemrgMeasure() };

    // These data are under MITK-Data in the Build directory
    array<pair<const char*, mitk::Surface::Pointer>, 5> surfaceData { {
        {MITK_DATA_DIR "/IGT-Data/Book.stl", nullptr},
        {MITK_DATA_DIR "/IGT-Data/ClaronTool.stl", nullptr},
        {MITK_DATA_DIR "/IGT-Data/EMTool.stl", nullptr},
        {MITK_DATA_DIR "/ball.stl", nullptr},
        {MITK_DATA_DIR "/binary.stl", nullptr}
    } };

private slots:
    void initTestCase();
    void cleanupTestCase();

    void CalcDistance_data();
    void CalcDistance();

    void CalcPerimeter_data();
    void CalcPerimeter();

    void CalcArea_data();
    void CalcArea();

    void FindCentre_data();
    void FindCentre();

    void GetSphericity_data();
    void GetSphericity();

    void calcVolumeMesh_data();
    void calcVolumeMesh();

    void calcSurfaceMesh_data();
    void calcSurfaceMesh();

    void Convert_data();
    void Convert();

    void Deconvert_data();
    void Deconvert();
};

Q_DECLARE_METATYPE(CemrgMeasure::Points)
Q_DECLARE_METATYPE(mitk::PointSet::Pointer)
Q_DECLARE_METATYPE(mitk::Point3D)
Q_DECLARE_METATYPE(vtkPolyData*)
Q_DECLARE_METATYPE(mitk::Surface::Pointer)
Q_DECLARE_METATYPE(mitk::DataNode::Pointer)
Q_DECLARE_METATYPE(string)
