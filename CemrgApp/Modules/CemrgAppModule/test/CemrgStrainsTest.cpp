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

#include "CemrgStrainsTest.hpp"

static bool FuzzyCompare(const double& lhs, const double& rhs) {
    constexpr double almostZero = 1e-6;
    if (abs(lhs) < almostZero && abs(rhs) < almostZero)
        return true;
    
    return qFuzzyCompare(lhs, rhs);
}

mitk::DataNode::Pointer TestCemrgStrains::ReferenceAHA(const array<int, 3>& segRatios, bool pacingSite) {
    mitk::PointSet::Pointer pointSet = mitk::IOUtil::Load<mitk::PointSet>((QFINDTESTDATA(CemrgTestData::strainPath) + "/PointSet.mps").toStdString());
    mitk::DataNode::Pointer lmNode = mitk::DataNode::New();
    lmNode->SetData(pointSet);
    cemrgStrains->ReferenceAHA(lmNode, (int*)segRatios.data(), pacingSite);
    return lmNode;
}

void TestCemrgStrains::initTestCase() {

}

void TestCemrgStrains::cleanupTestCase() {

}

void TestCemrgStrains::CalculateGlobalSqzPlot_data() {
    QTest::addColumn<int>("meshNo");
    QTest::addColumn<double>("result");

    const array<double, CemrgTestData::strainDataSize> globalSqzPlotData {
        0,
        -0.08891526089125015297
    };

    for (size_t i = 0; i < globalSqzPlotData.size(); i++)
        QTest::newRow(("Test " + to_string(i + 1)).c_str()) << (int)i << globalSqzPlotData[i];
}

void TestCemrgStrains::CalculateGlobalSqzPlot() {
    QFETCH(int, meshNo);
    QFETCH(double, result);

    QCOMPARE(cemrgStrains->CalculateGlobalSqzPlot(meshNo), result);
}

void TestCemrgStrains::CalculateSqzPlot_data() {
    QTest::addColumn<int>("meshNo");
    QTest::addColumn<vector<double>>("result");

    const array<vector<double>, CemrgTestData::strainDataSize> sqzPlotData { {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {-0.24161484861895057841, -0.24203718408602067913, -0.24189488048215221361, -0.24269039612066256595, -0.23629833084106455221, -0.23554401643925146348, -0.23956905910100148582, -0.24574671778317599968, -0.24545802466742852599, -0.24045512800893234506, -0.23522178212071434555, -0.23486487750507906158, -0.24153056562196148493, -0.23671392457799511622, -0.24116868663725185562, -0.23596466822878597869}
    } };

    // Preparation for tests
    ReferenceAHA();

    for (size_t i = 0; i < sqzPlotData.size(); i++)
        QTest::newRow(("Test " + to_string(i + 1)).c_str()) << (int)i << sqzPlotData[i];
}

void TestCemrgStrains::CalculateSqzPlot() {
    QFETCH(int, meshNo);
    QFETCH(vector<double>, result);

    QVERIFY(equal(begin(result), end(result), begin(cemrgStrains->CalculateSqzPlot(meshNo)), FuzzyCompare));
}

void TestCemrgStrains::CalculateStrainsPlot_data() {
    QTest::addColumn<int>("meshNo");
    QTest::addColumn<mitk::DataNode::Pointer>("lmNode");
    QTest::addColumn<int>("flag");
    QTest::addColumn<vector<double>>("result");

    const array<vector<double>, CemrgTestData::strainDataSize> strainsPlotData { {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {-0.45874886139987286482, -0.38998985665445839999, -0.39095134301528899901, -0.45061430609189162544, -0.48968384803220571522, -0.49392625117662009027, -0.45850249287033251200, -0.35280369788968640732, -0.35498114762612859030, -0.44960627315573070684, -0.49382093507974794688, -0.49665076176258193819, -0.31543513441107573492, -0.43614002435635201849, -0.32295055120045534913, -0.46667719896202092267}
    } };

    // Preparation for tests
    mitk::DataNode::Pointer lmNode = ReferenceAHA();

    for (size_t i = 0; i < strainsPlotData.size(); i++)
        QTest::newRow(("Test " + to_string(i + 1)).c_str()) << (int)i << lmNode << (int)i << strainsPlotData[i];
}

void TestCemrgStrains::CalculateStrainsPlot() {
    QFETCH(int, meshNo);
    QFETCH(mitk::DataNode::Pointer, lmNode);
    QFETCH(int, flag);
    QFETCH(vector<double>, result);

    QVERIFY(equal(begin(result), end(result), begin(cemrgStrains->CalculateStrainsPlot(meshNo, lmNode, flag)), FuzzyCompare));
}

void TestCemrgStrains::CalculateSDI_data() {
    QTest::addColumn<vector<vector<double>>>("valueVectors");
    QTest::addColumn<int>("cycleLengths");
    QTest::addColumn<int>("noFrames");
    QTest::addColumn<double>("result");

    const array<double, 5> sdiData {
        0,
        25,
        16.65,
        12.5,
        10
    };
    
    // Preparation for tests
    ReferenceAHA();

    vector<vector<double>> valueVectors;
    for (size_t i = 0; i < sdiData.size(); i++) {
        valueVectors.push_back(cemrgStrains->CalculateSqzPlot(i % CemrgTestData::strainDataSize));
        QTest::newRow(("Test " + to_string(i + 1)).c_str()) << valueVectors << 1000 << (int)i + 1 << sdiData[i];
    }
}

void TestCemrgStrains::CalculateSDI() {
    QFETCH(vector<vector<double>>, valueVectors);
    QFETCH(int, cycleLengths);
    QFETCH(int, noFrames);
    QFETCH(double, result);

    QCOMPARE(cemrgStrains->CalculateSDI(valueVectors, cycleLengths, noFrames), result);
}

void TestCemrgStrains::ReferenceGuideLines_data() {
    QTest::addColumn<mitk::DataNode::Pointer>("lmNode");
    QTest::addColumn<vector<mitk::Surface::Pointer>>("result");

    const array<vector<tuple<double, double, double>>, 5> pointSetData { {
        {
            {85.4995566059, 85.6483483021, 91.9406230158},
            {23.1505777288, 132.670464778, 141.375157912},
            {44.8106937268, 158.193900374, 135.289413785},
            {20.4255620648, 144.034012592, 115.474941806}
        },
        {
            {85.4995566059, 85.6483483021, 91.9406230158},
            {23.1505777288, 132.670464778, 141.375157912},
            {44.8106937268, 158.193900374, 135.289413785},
            {20.4255620648, 144.034012592, 115.474941806},
            {37.8195032758, 119.896131031, 91.3032690028},
            {49.2203653858, 104.544360923, 140.767393939}
        },
        {
            {-44.8106937268, -158.193900374, -135.289413785},
            {20.4255620648, 144.034012592, 115.474941806},
            {-37.8195032758, -119.896131031, -91.3032690028},
            {49.2203653858, 104.544360923, 140.767393939}
        },
        {
            {-85.4995566059, 85.6483483021, -91.9406230158},
            {-23.1505777288, -132.670464778, -141.375157912},
            {37.8195032758, 119.896131031, 91.3032690028},
            {49.2203653858, -104.544360923, 140.767393939}
        },
        {
            {20.4255620648, 144.034012592, 115.474941806},
            {37.8195032758, 119.896131031, 91.3032690028},
            {49.2203653858, 104.544360923, 140.767393939},
            {85.4995566059, 85.6483483021, 91.9406230158},
            {23.1505777288, 132.670464778, 141.375157912},
            {44.8106937268, 158.193900374, 135.289413785}
        }
    } };

    const array<array<double[2][3], 3>, 5> referenceGuideLinesData { {
        { { { {85.49955749511718750000, 85.64834594726562500000, 91.94062042236328125000}, {23.15057754516601562500, 132.67047119140625000000, 141.37515258789062500000} }, { {44.81069183349609375000, 158.19389343261718750000, 135.28941345214843750000}, {26.44359207153320312500, 130.18695068359375000000, 138.76423645019531250000} }, { {20.42556190490722656250, 144.03401184082031250000, 115.47494506835937500000}, {27.35565567016601562500, 129.49909973144531250000, 138.04109191894531250000} } } },
        { { { {85.49955749511718750000, 85.64834594726562500000, 91.94062042236328125000}, {30.88863945007324218750, 146.29669189453125000000, 130.95314025878906250000} }, { {37.81950378417968750000, 119.89613342285156250000, 91.30326843261718750000}, {54.42469787597656250000, 120.15864562988281250000, 114.13963317871093750000} }, { {49.22036361694335937500, 104.54435729980468750000, 140.76739501953125000000}, {51.91492080688476562500, 122.94588470458984375000, 115.93254852294921875000} } } },
        { { { {-44.81069183349609375000, -158.19389343261718750000, -135.28941345214843750000}, {20.42556190490722656250, 144.03401184082031250000, 115.47494506835937500000} }, { {-37.81950378417968750000, -119.89613342285156250000, -91.30326843261718750000}, {-35.31798934936523437500, -114.21589660644531250000, -98.80001068115234375000} }, { {49.22036361694335937500, 104.54435729980468750000, 140.76739501953125000000}, {18.89675903320312500000, 136.95133972167968750000, 109.59831237792968750000} } } },
        { { { {-85.49955749511718750000, 85.64834594726562500000, -91.94062042236328125000}, {-23.15057754516601562500, -132.67047119140625000000, -141.37515258789062500000} }, { {37.81950378417968750000, 119.89613342285156250000, 91.30326843261718750000}, {-95.71512603759765625000, 121.41880798339843750000, -83.84101867675781250000} }, { {49.22036361694335937500, -104.54435729980468750000, 140.76739501953125000000}, {-41.13645935058593750000, -69.69178771972656250000, -127.11472320556640625000} } } },
        { { { {20.42556190490722656250, 144.03401184082031250000, 115.47494506835937500000}, {60.32823562622070312500, 101.23760223388671875000, 108.75263977050781250000} }, { {23.15057754516601562500, 132.67047119140625000000, 141.37515258789062500000}, {25.26762199401855468750, 138.84080505371093750000, 114.65921020507812500000} }, { {44.81069183349609375000, 158.19389343261718750000, 135.28941345214843750000}, {23.11538887023925781250, 141.14912414550781250000, 115.02178955078125000000} } } }
    } };

    for (size_t i = 0; i < referenceGuideLinesData.size(); i++) {
        // Prepare the input
        mitk::PointSet::Pointer pointSet = mitk::PointSet::New();
        mitk::Point3D point;
        for (auto& pt : pointSetData[i]) {
            tie(point[0], point[1], point[2]) = pt;
            pointSet->InsertPoint(point);
        }
        mitk::DataNode::Pointer lmNode = mitk::DataNode::New();
        lmNode->SetData(pointSet);

        // Prepare the output
        vector<mitk::Surface::Pointer> surfaces;
        for (size_t j = 0; j < referenceGuideLinesData[i].size(); j++) {
            vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
            lineSource->SetPoint1(referenceGuideLinesData[i][j][0][0], referenceGuideLinesData[i][j][0][1], referenceGuideLinesData[i][j][0][2]);
            lineSource->SetPoint2(referenceGuideLinesData[i][j][1][0], referenceGuideLinesData[i][j][1][1], referenceGuideLinesData[i][j][1][2]);
            lineSource->Update();
            mitk::Surface::Pointer surface = mitk::Surface::New();
            surface->SetVtkPolyData(lineSource->GetOutput());
            surfaces.push_back(surface);
        }

        QTest::newRow(("Test " + to_string(i + 1)).c_str()) << lmNode << surfaces;
    }
}

void TestCemrgStrains::ReferenceGuideLines() {
    QFETCH(mitk::DataNode::Pointer, lmNode);
    QFETCH(vector<mitk::Surface::Pointer>, result);

    // for (auto& surface : cemrgStrains->ReferenceGuideLines(lmNode)) {
    //     for (size_t i = 0; i < surface->GetSizeOfPolyDataSeries(); i++) {
    //         for (int j = 0; j < surface->GetVtkPolyData(i)->GetPoints()->GetNumberOfPoints(); j++) {
    //             double point[3];
    //             surface->GetVtkPolyData(i)->GetPoints()->GetPoint(j, point);
    //             qDebug() << fixed << qSetRealNumberPrecision(20) << point[0] << point[1] << point[2];
    //         }
    //     }
    // }
    
    auto referenceGuideLines = cemrgStrains->ReferenceGuideLines(lmNode);
    for (size_t i = 0; i < result.size(); i++)
        QVERIFY(mitk::Equal(*referenceGuideLines[i], *result[i], mitk::eps, true));
}

int CemrgStrainsTest(int argc, char *argv[]) {
    QCoreApplication app(argc, argv);
    app.setAttribute(Qt::AA_Use96Dpi, true);
    TestCemrgStrains tc;
    QTEST_SET_MAIN_SOURCE_PATH
    return QTest::qExec(&tc, argc, argv);
}
