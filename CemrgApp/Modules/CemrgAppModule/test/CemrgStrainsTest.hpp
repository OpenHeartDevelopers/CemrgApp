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
#include <CemrgStrains.h>

using namespace std;

class TestCemrgStrains: public QObject {

    Q_OBJECT

private:
    unique_ptr<CemrgStrains> cemrgStrains { new CemrgStrains(QFINDTESTDATA(CemrgTestData::strainPath), 0) };

    // Used for preparation of multiple tests
    mitk::DataNode::Pointer ReferenceAHA(const array<int, 3>& segRatios = { 40, 40, 20 }, bool pacingSite = false);

private slots:
    void initTestCase();
    void cleanupTestCase();

    void CalculateGlobalSqzPlot_data();
    void CalculateGlobalSqzPlot();

    void CalculateSqzPlot_data();
    void CalculateSqzPlot();

    void CalculateStrainsPlot_data();
    void CalculateStrainsPlot();

    void CalculateSDI_data();
    void CalculateSDI();

    void ReferenceGuideLines_data();
    void ReferenceGuideLines();
};

Q_DECLARE_METATYPE(vector<double>)
Q_DECLARE_METATYPE(mitk::DataNode::Pointer)
Q_DECLARE_METATYPE(vector<vector<double>>)
Q_DECLARE_METATYPE(vector<mitk::Surface::Pointer>)
