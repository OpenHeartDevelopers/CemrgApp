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

#include "CemrgMeasureTest.hpp"

typedef CemrgMeasure::Point Point;
typedef CemrgMeasure::Points Points;

void TestCemrgMeasure::initTestCase() {
    // Create surface data
    for (size_t i = 0; i < surfaceData.size(); i++) {
        surfaceData[i].first = QFINDTESTDATA(CemrgTestData::surfacePaths[i]);
        surfaceData[i].second = mitk::IOUtil::Load<mitk::Surface>(surfaceData[i].first.toStdString());
    }
}

void TestCemrgMeasure::cleanupTestCase() {

}

void TestCemrgMeasure::CalcDistance_data() {
    QTest::addColumn<Points>("points");
    QTest::addColumn<double>("result");

    QTest::newRow("1-point test") << Points{{0.0, 0.0, 0.0}} << -1.0;
    QTest::newRow("3-point test") << Points{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {2.0, 2.0, 2.0}} << -1.0;
    QTest::newRow("Test 1") << Points{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}} << 1.7320508075688772935274463415059;
    QTest::newRow("Test 2") << Points{{7.0, 4.0, 3.0}, {17.0, 6.0, 2.0}} << 10.246950765959598383221038680521;
    QTest::newRow("Test 3") << Points{{-4.0, 3.0, -2.0}, {23.0, -2.0, 5.0}} << 28.337254630609507934884031143657;
}

void TestCemrgMeasure::CalcDistance() {
    QFETCH(Points, points);
    QFETCH(double, result);

    QCOMPARE(cemrgMeasure->CalcDistance(points), result);
}

void TestCemrgMeasure::CalcPerimeter_data() {
    QTest::addColumn<Points>("points");
    QTest::addColumn<double>("result");

    QTest::newRow("1-point test") << Points{{0.0, 0.0, 0.0}} << -1.0;
    QTest::newRow("2-point test") << Points{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}} << -1.0;
    QTest::newRow("3-point test") << Points{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}} << 3.4641016151377545870548926830117;
    QTest::newRow("4-point test") << Points{{-1.0, 5.0, 7.0}, {2.0, 3.0, 4.0}, {21.0, -15.0, -20.0}, {-1.0, 5.0, 7.0}} << 80.363148824999240938914956776348;
    QTest::newRow("5-point test") << Points{{-13.0, 8.0, -7.0}, {1.0, 1.0, -1.0}, {-21.0, 15.0, 20.0}, {12.0, 51.0, 72.0}, {-13.0, 8.0, -7.0}} << 214.93578434989057588768898364186;
}

void TestCemrgMeasure::CalcPerimeter() {
    QFETCH(Points, points);
    QFETCH(double, result);

    QCOMPARE(cemrgMeasure->CalcPerimeter(points), result);
}

void TestCemrgMeasure::CalcArea_data() {
    QTest::addColumn<Points>("points");
    QTest::addColumn<double>("result");

    QTest::newRow("1-point test") << Points{{0.0, 0.0, 0.0}} << -1.0;
    QTest::newRow("2-point test") << Points{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}} << -1.0;
    QTest::newRow("3-point test") << Points{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0}} << 0.0000000160123397681644855524243;
    QTest::newRow("4-point test") << Points{{-1.0, 5.0, 7.0}, {2.0, 3.0, 4.0}, {21.0, -15.0, -20.0}, {-1.0, 5.0, 7.0}} << 11.3688170009044213770721398759633;
    QTest::newRow("5-point test") << Points{{-13.0, 8.0, -7.0}, {1.0, 1.0, -1.0}, {-21.0, 15.0, 20.0}, {12.0, 51.0, 72.0}, {-13.0, 8.0, -7.0}} << 1016.3855905904007386197918094694614;
}

void TestCemrgMeasure::CalcArea() {
    QFETCH(Points, points);
    QFETCH(double, result);

    QCOMPARE(cemrgMeasure->CalcArea(points), result);
}

void TestCemrgMeasure::FindCentre_data() {
    QTest::addColumn<mitk::PointSet::Pointer>("pointSet");
    QTest::addColumn<mitk::Point3D>("result");

    // Point and Result values
    const array<tuple<Point, Point>, 8> findCentreData { {
        { {0, 0, 0}, {0, 0, 0} },
        { {1, 1, 1}, {0.5, 0.5, 0.5} },
        { {-1, 5, 7}, {0, 2, 2.66666666666666651864} },
        { {2, 3, 4}, {0.5, 2.25, 3} },
        { {21, -15, -20}, {4.6, -1.2, -1.6} },
        { {-13, 8, -7}, {1.66666666666666674068, 0.33333333333333331483, -2.5} },
        { {12, 51, 72}, {3.14285714285714279370, 7.57142857142857117481, 8.14285714285714234961} },
        { {1, 1, -1}, {2.875, 6.75, 7} } 
    } };

    mitk::PointSet::Pointer pointSet = mitk::PointSet::New();
    mitk::Point3D point;
    Point pt, result;
    for (size_t i = 0; i < findCentreData.size(); i++) {
        tie(pt, result) = findCentreData[i];
        pointSet = pointSet->Clone();
        tie(point[0], point[1], point[2]) = pt;
        pointSet->InsertPoint(point);
        tie(point[0], point[1], point[2]) = result;
        QTest::newRow((to_string(i + 1) + "-point test").c_str()) << pointSet << point;
    }
}

void TestCemrgMeasure::FindCentre() {
    QFETCH(mitk::PointSet::Pointer, pointSet);
    QFETCH(mitk::Point3D, result);

    QCOMPARE(cemrgMeasure->FindCentre(pointSet), result);
}

void TestCemrgMeasure::GetSphericity_data() {
    QTest::addColumn<vtkPolyData*>("polyData");
    QTest::addColumn<double>("result");

    const array<double, CemrgTestData::surfacePaths.size()> sphericityData {
        65.06623281457561347452,
        20.59264742863982178278,
        30.61636022943845603095,
        99.64016354530460262140,
        66.26402857707076066163
    };

    for (size_t i = 0; i < sphericityData.size(); i++)
        QTest::newRow(("Test " + QFileInfo(surfaceData[i].first).fileName().toStdString()).c_str()) << surfaceData[i].second->GetVtkPolyData() << sphericityData[i];
}

void TestCemrgMeasure::GetSphericity() {
    QFETCH(vtkPolyData*, polyData);
    QFETCH(double, result);

    QCOMPARE(cemrgMeasure->GetSphericity(polyData), result);
}

void TestCemrgMeasure::calcVolumeMesh_data() {
    QTest::addColumn<mitk::Surface::Pointer>("surface");
    QTest::addColumn<double>("result");

    const array<double, CemrgTestData::surfacePaths.size()> volumeMeshData {
        1253966.56670495495200157166,
        6649.51557345613582583610,
        1258.14200772278718432062,
        13756.31573527904947695788,
        108985.12497850439103785902
    };

    for (size_t i = 0; i < volumeMeshData.size(); i++)
        QTest::newRow(("Test " + QFileInfo(surfaceData[i].first).fileName().toStdString()).c_str()) << surfaceData[i].second << volumeMeshData[i];
}

void TestCemrgMeasure::calcVolumeMesh() {
    QFETCH(mitk::Surface::Pointer, surface);
    QFETCH(double, result);

    QCOMPARE(cemrgMeasure->calcVolumeMesh(surface), result);
}

void TestCemrgMeasure::calcSurfaceMesh_data() {
    QTest::addColumn<mitk::Surface::Pointer>("surface");
    QTest::addColumn<double>("result");

    const array<double, CemrgTestData::surfacePaths.size()> surfaceMeshData {
        87203.96462370984954759479,
        5998.53801315712826180970,
        1746.24837584279430302558,
        2784.97791560359246432199,
        14338.23953928404807811603
    };

    for (size_t i = 0; i < surfaceMeshData.size(); i++)
        QTest::newRow(("Test " + QFileInfo(surfaceData[i].first).fileName().toStdString()).c_str()) << surfaceData[i].second << surfaceMeshData[i];
}

void TestCemrgMeasure::calcSurfaceMesh() {
    QFETCH(mitk::Surface::Pointer, surface);
    QFETCH(double, result);

    QCOMPARE(cemrgMeasure->calcSurfaceMesh(surface), result);
}


/*****************************************************************************************************/
/***************************************  Conversion Functions ***************************************/
/*****************************************************************************************************/
constexpr size_t conversionDataSize = 5;

const array<const string, conversionDataSize> conversionFileData {
    "# vtk DataFile Version 4.2\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS 12 float\n2.84457 7.69388 30.6308 2.89832 7.45175 30.697 2.99997 6.97862 30.814 \n3.08689 6.58059 30.9175 3.12502 6.39859 30.959 3.23261 5.9521 31.1119 \n3.26323 5.83054 31.1584 3.27528 5.77604 31.1731 3.34243 5.48657 31.2628 \n3.43975 5.07282 31.3958 3.54421 4.63671 31.5429 3.61901 4.32759 31.6499 \n",
    "# vtk DataFile Version 4.2\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS 5 float\n37.7279 5.9931 23.0381 37.9221 6.38497 22.7362 38.1279 6.98015 22.4036\n38.1583 7.11666 22.3511 38.3233 7.97394 22.0579",
    "# vtk DataFile Version 4.2\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS 3 float\n19.7627 4.36692 23.3719 19.7856 3.95139 23.4737 19.7912 3.83002 23.5085",
    "# vtk DataFile Version 4.2\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS 2 float\n37.7279 5.9931 23.0381 37.9221 6.38497 22.7362 \n",
    "# vtk DataFile Version 4.2\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS 4 float\n2.84457 7.69388 30.6308 2.89832 7.45175 30.697 2.99997 6.97862 30.814 \n3.08689 6.58059 30.9175\n"
};

const array<const Points, conversionDataSize> conversionPointData { {
    { {-2.84457, -7.69388, 30.6308}, {-2.89832, -7.45175, 30.697}, {-2.99997, -6.97862, 30.814}, {-3.08689, -6.58059, 30.9175}, {-3.12502, -6.39859, 30.959}, {-3.23261, -5.9521, 31.1119}, {-3.26323, -5.83054, 31.1584}, {-3.27528, -5.77604, 31.1731}, {-3.34243, -5.48657, 31.2628}, {-3.43975, -5.07282, 31.3958}, {-3.54421, -4.63671, 31.5429}, {-3.61901, -4.32759, 31.6499} },
    { {-37.7279, -5.9931, 23.0381}, {-37.9221, -6.38497, 22.7362}, {-38.1279, -6.98015, 22.4036}, {-38.1583, -7.11666, 22.3511}, {-38.3233, -7.97394, 22.0579} },
    { {-19.7627, -4.36692, 23.3719}, {-19.7856, -3.95139, 23.4737}, {-19.7912, -3.83002, 23.5085} },
    { {-37.7279, -5.9931, 23.0381}, {-37.9221, -6.38497, 22.7362} },
    { {-2.84457, -7.69388, 30.6308}, {-2.89832, -7.45175, 30.697}, {-2.99997, -6.97862, 30.814}, {-3.08689, -6.58059, 30.9175} }
} };

void TestCemrgMeasure::Convert_data() {
    QTest::addColumn<QString>("dir");
    QTest::addColumn<mitk::DataNode::Pointer>("dataNode");
    QTest::addColumn<string>("result");

    const QString dir = "./";

    // Create DataNodes
    array<mitk::DataNode::Pointer, conversionDataSize> dataNode;
    for (size_t i = 0; i < conversionDataSize; i++) {
        mitk::PointSet::Pointer pointSet = mitk::PointSet::New();
        for (size_t j = 0; j < conversionPointData[i].size(); j++) {
            mitk::Point3D point;
            tie(point[0], point[1], point[2]) = conversionPointData[i][j];
            pointSet->InsertPoint(point);
        }
        dataNode[i] = mitk::DataNode::New();
        dataNode[i]->SetData(pointSet);
    }

    for (size_t i = 0; i < conversionDataSize; i++)
        QTest::newRow(("Test " + to_string(i + 1)).c_str()) << dir << dataNode[i] << conversionFileData[i];
}

void TestCemrgMeasure::Convert() {
    QFETCH(QString, dir);
    QFETCH(mitk::DataNode::Pointer, dataNode);
    QFETCH(string, result);

    cemrgMeasure->Convert(dir, dataNode);

    // Compare the created file with result
    ifstream convertedFile;
    convertedFile.open(dir.toStdString() + "input.vtk");
    QVERIFY(convertedFile.is_open());
    istringstream resultFile(result);

    // Skip first lines
    string convertedLine, resultLine;
    getline(convertedFile, convertedLine); getline(resultFile, resultLine);
    QVERIFY(convertedFile.good());

    for (int i = 2; convertedFile.good() && resultFile.good(); i++) {
        getline(convertedFile, convertedLine);
        getline(resultFile, resultLine);
        QVERIFY2(QString(convertedLine.c_str()).trimmed() == QString(resultLine.c_str()).trimmed(), ("Line " + to_string(i) + " doesn't match!").c_str());
    }
}

void TestCemrgMeasure::Deconvert_data() {
    QTest::addColumn<QString>("dir");
    QTest::addColumn<int>("fileNo");
    QTest::addColumn<Points>("result");

    const QString dir = "./";

    // Write test data into files
    for (size_t i = 0; i < conversionDataSize; i++) {
        ofstream file;
        file.open(dir.toStdString() + "/transformed-" + to_string(i) + ".vtk");
        file << conversionFileData[i];
    }

    for (size_t i = 0; i < conversionDataSize; i++)
        QTest::newRow(("Test " + to_string(i + 1)).c_str()) << dir << (int)i << conversionPointData[i];
}

void TestCemrgMeasure::Deconvert() {
    QFETCH(QString, dir);
    QFETCH(int, fileNo);
    QFETCH(Points, result);

    QCOMPARE(cemrgMeasure->Deconvert(dir, fileNo), result);
}

int CemrgMeasureTest(int argc, char *argv[]) {
    QCoreApplication app(argc, argv);
    app.setAttribute(Qt::AA_Use96Dpi, true);
    TestCemrgMeasure tc;
    QTEST_SET_MAIN_SOURCE_PATH
    return QTest::qExec(&tc, argc, argv);
}
