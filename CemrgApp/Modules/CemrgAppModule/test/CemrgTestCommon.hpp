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
#include <mitkDataNode.h>
#include <mitkImageCast.h>
#include <mitkImagePixelReadAccessor.h>

// VTK
#include <vtkPolyDataNormals.h>
#include <vtkLineSource.h>

using namespace std;

struct CemrgTestData {
    static constexpr array<const char*, 5> surfacePaths {
        "Data/Surface/Book.stl",
        "Data/Surface/ClaronTool.stl",
        "Data/Surface/EMTool.stl",
        "Data/Surface/ball.stl",
        "Data/Surface/binary.stl"
    };

    static constexpr const char *strainPath = "Data/Strain";
    static constexpr size_t strainDataSize = 2;

    static constexpr const char *cmdLinePath = "Data/CommandLine";
};
