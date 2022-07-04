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

// CemrgApp
#include "CemrgTestCommon.hpp"
#include <CemrgCommandLine.h>

using namespace std;

class TestCemrgCommandLine: public QObject {

    Q_OBJECT

private:
    unique_ptr<CemrgCommandLine> cemrgCommandLine { new CemrgCommandLine() };

    const QString dataPath = QFINDTESTDATA(CemrgTestData::cmdLinePath);

private slots:
    void initTestCase();
    void cleanupTestCase();

    void ExecuteSurf_data();
    void ExecuteSurf();

    void ExecuteRegistration_data();
    void ExecuteRegistration();

    void ExecuteTransformation_data();
    void ExecuteTransformation();

    void ExecuteSimpleTranslation_data();
    void ExecuteSimpleTranslation();

    void ExecuteTracking_data();
    void ExecuteTracking();

    void ExecuteApplying_data();
    void ExecuteApplying();

    void ExecuteCreateCGALMesh_data();
    void ExecuteCreateCGALMesh();
};
