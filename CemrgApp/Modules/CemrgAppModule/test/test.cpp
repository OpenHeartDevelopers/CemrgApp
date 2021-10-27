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

#include "test.hpp"

bool TestSuite::myCondition() {
    return true;
}

void TestSuite::initTestCase() {
    qDebug("Called before everything else.");
}

void TestSuite::myFirstTest() {
    QVERIFY(true);  // check that a condition is satisfied
    QCOMPARE(1, 1); // compare two values
}

void TestSuite::mySecondTest() {
    QVERIFY(myCondition());
    QVERIFY(1 != 2);
}

void TestSuite::cleanupTestCase() {
    qDebug("Called after myFirstTest and mySecondTest.");
}

// TODO: Function name should be the same as the filename
// Therefore, MITK test system can call it
int test(int argc, char *argv[]) {
    QCoreApplication app(argc, argv);
    app.setAttribute(Qt::AA_Use96Dpi, true);
    // TODO: Don't forget to change the class name here
    TestSuite tc;
    QTEST_SET_MAIN_SOURCE_PATH
    return QTest::qExec(&tc, argc, argv);
}
