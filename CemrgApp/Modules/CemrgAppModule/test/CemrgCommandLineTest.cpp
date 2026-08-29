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

#include "CemrgCommandLineTest.hpp"

static bool EqualFiles(const QString& filename1, const QString& filename2) {
    // Open files at the end
    ifstream file1(filename1.toStdString(), ifstream::ate | ifstream::binary);
    ifstream file2(filename2.toStdString(), ifstream::ate | ifstream::binary);

    // Check if files are opened
    if (!file1.is_open() || !file2.is_open())
        return false;

    // Different file size
    if (file1.tellg() != file2.tellg())
        return false;

    // Rewind
    file1.seekg(0);
    file2.seekg(0);

    return equal(istreambuf_iterator<char>(file1), istreambuf_iterator<char>(), istreambuf_iterator<char>(file2));
}

static bool CopyDirectory(const QString& src, const QString& dst) {
    // Check if the source directory exists
    QDir dir(src);
    if (!dir.exists())
        return false;

    // Create the destination directory
    dir.mkpath(dst);

    // Copy files
    for (auto& fileName : dir.entryList(QDir::Files)) {
        const QString destFilePath = dst + "/" + fileName;
        if (!QFileInfo(destFilePath).exists()) {
            if (!QFile::copy(src + "/" + fileName, destFilePath))
                return false;
        }
    }

    return true;
}
/*
static QString PrepareSegmentationForCGALMesh(QString dir, QString segmentationFileName) {
    QFileInfo segFileInfo(dir + "/" + segmentationFileName);
    QString preparedSegmentationPath = dir + "/" + segFileInfo.baseName() + ".inr";
    mitk::Image::Pointer image = mitk::IOUtil::Load<mitk::Image>(segFileInfo.absoluteFilePath().toStdString());
    if (!image) {
        qDebug() << "Segmentation file couldn't be loaded: " << segmentationFileName;
        return {};
    }
    try {
        // Convert image to right type
        itk::Image<uint8_t, 3>::Pointer itkImage = itk::Image<uint8_t, 3>::New();
        mitk::CastToItkImage(image, itkImage);
        mitk::CastToMitkImage(itkImage, image);
        // Access image volume
        mitk::ImagePixelReadAccessor<uint8_t, 3> readAccess(image);
        uint8_t *pv = (uint8_t*)readAccess.GetData();
        // Prepare header of inr file
        char header[256];
        const int bitlength = 8;
        const char *btype = "unsigned fixed";
        mitk::Vector3D spacing = image->GetGeometry()->GetSpacing();
        int n = sprintf(header, "#INRIMAGE-4#{\nXDIM=%d\nYDIM=%d\nZDIM=%d\nVDIM=1\nTYPE=%s\nPIXSIZE=%d bits\nCPU=decm\nVX=%6.4f\nVY=%6.4f\nVZ=%6.4f\n", image->GetDimension(0), image->GetDimension(1), image->GetDimension(2), btype, bitlength, spacing.GetElement(0), spacing.GetElement(1), spacing.GetElement(2));
        for (int i = n; i < 252; i++)
            header[i] = '\n';
        header[252] = '#';
        header[253] = '#';
        header[254] = '}';
        header[255] = '\n';
        // Write to binary file
        ofstream segmentationFile(preparedSegmentationPath.toStdString(), ios::out | ios::binary);
        segmentationFile.write(header, sizeof(header));
        segmentationFile.write((char*)pv, image->GetDimension(0) * image->GetDimension(1) * image->GetDimension(2));
        segmentationFile.close();
    } catch (mitk::Exception& ex) {
        qDebug() << ex.GetDescription();
    }
    return QFileInfo(preparedSegmentationPath).fileName();
}
*/
void TestCemrgCommandLine::initTestCase() {
    // Copy MIRTK libraries
    QVERIFY(CopyDirectory("/Externals/MLib", QCoreApplication::applicationDirPath() + "/MLib"));
    QVERIFY(CopyDirectory("/Externals/M3DLib", QCoreApplication::applicationDirPath() + "/M3DLib"));

    // Check the test data directory
    QVERIFY2(QDir(dataPath).exists(), "The test data directory doesn't exist!");
}

void TestCemrgCommandLine::cleanupTestCase() {

}

void TestCemrgCommandLine::ExecuteSurf_data() {
    QTest::addColumn<QString>("segPath");
    QTest::addColumn<float>("threshold");
    QTest::addColumn<int>("blur");
    QTest::addColumn<int>("smoothness");
    QTest::addColumn<QString>("result");

    const array<tuple<QString, QString, int, float, int, int, QString>, 2> surfData { {
        {"sphere_initial.nii", 0.5, 0, 10, "/surf_expected_1.vtk"},
        {"sphere_shifted.nii", 0.5, 0, 10, "/surf_expected_2.vtk"},
    } };

    for (size_t i = 0; i < surfData.size(); i++)
        QTest::newRow(("Test " + to_string(i + 1)).c_str()) << get<0>(surfData[i]) << get<1>(surfData[i]) << get<2>(surfData[i]) << get<3>(surfData[i]) << get<4>(surfData[i]) << get<5>(surfData[i]) << get<6>(surfData[i]);
}

void TestCemrgCommandLine::ExecuteSurf() {
    QFETCH(QString, segPath);
    QFETCH(float, threshold);
    QFETCH(int, blur);
    QFETCH(int, smoothness);
    QFETCH(QString, result);

    QString surfOutput = cemrgCommandLine->ExecuteSurf(dataPath, segPath, threshold, blur, smoothness);
    QVERIFY2(EqualFiles(surfOutput, dataPath + result), "The function output is different from the expected output!");
}

void TestCemrgCommandLine::ExecuteRegistration_data() {
    QTest::addColumn<QString>("fixedFileName");
    QTest::addColumn<QString>("movingFileName");
    QTest::addColumn<QString>("transformFileName");
    QTest::addColumn<QString>("modelName");
    QTest::addColumn<QString>("result");

    const array<QString[5], 1> registrationData { {
        {"sphere_initial.nii", "sphere_shifted.nii", "reg_output.dof", "Rigid", "/reg_expected.dof"}
    } };

    for (size_t i = 0; i < registrationData.size(); i++)
        QTest::newRow(("Test " + to_string(i + 1)).c_str()) << registrationData[i][0] << registrationData[i][1] << registrationData[i][2] << registrationData[i][3] << registrationData[i][4];
}

void TestCemrgCommandLine::ExecuteRegistration() {
    QFETCH(QString, fixedFileName);
    QFETCH(QString, movingFileName);
    QFETCH(QString, transformFileName);
    QFETCH(QString, modelName);
    QFETCH(QString, result);

    cemrgCommandLine->ExecuteRegistration(dataPath, fixedFileName, movingFileName, transformFileName, modelName);
    QVERIFY2(EqualFiles(dataPath + "/" + transformFileName, dataPath + result), "The function output is different from the expected output!");
}

void TestCemrgCommandLine::ExecuteTransformation_data() {
    QTest::addColumn<QString>("imageFileName");
    QTest::addColumn<QString>("outputFileName");
    QTest::addColumn<QString>("transformFileName");
    QTest::addColumn<QString>("result");

    const array<QString[4], 1> transformationData { {
        {"sphere_initial.nii", "transformation_output.nii", "reg_expected.dof", "/transformation_expected.nii"}
    } };

    for (size_t i = 0; i < transformationData.size(); i++)
        QTest::newRow(("Test " + to_string(i + 1)).c_str()) << transformationData[i][0] << transformationData[i][1] << transformationData[i][2] << transformationData[i][3];
}

void TestCemrgCommandLine::ExecuteTransformation() {
    QFETCH(QString, imageFileName);
    QFETCH(QString, outputFileName);
    QFETCH(QString, transformFileName);
    QFETCH(QString, result);

    cemrgCommandLine->ExecuteTransformation(dataPath, imageFileName, outputFileName, transformFileName);
    QVERIFY2(EqualFiles(dataPath + "/" + outputFileName, dataPath + result), "The function output is different from the expected output!");
}

void TestCemrgCommandLine::ExecuteSimpleTranslation_data() {
    QTest::addColumn<QString>("fixedFileName");
    QTest::addColumn<QString>("movingFileName");
    QTest::addColumn<QString>("transformFileName");
    QTest::addColumn<bool>("transformPoints");
    QTest::addColumn<QString>("result");

    const array<tuple<QString, QString, QString, bool, QString>, 1> translationData { {
        {"surf_expected_1.vtk", "surf_expected_2.vtk", "translation_output.dof", true, "/translation_expected.dof"}
    } };

    for (size_t i = 0; i < translationData.size(); i++)
        QTest::newRow(("Test " + to_string(i + 1)).c_str()) << get<0>(translationData[i]) << get<1>(translationData[i]) << get<2>(translationData[i]) << get<3>(translationData[i]) << get<4>(translationData[i]);
}

void TestCemrgCommandLine::ExecuteSimpleTranslation() {
    QFETCH(QString, fixedFileName);
    QFETCH(QString, movingFileName);
    QFETCH(QString, transformFileName);
    QFETCH(bool, transformPoints);
    QFETCH(QString, result);

    cemrgCommandLine->ExecuteSimpleTranslation(dataPath, fixedFileName, movingFileName, transformFileName, transformPoints);
    QVERIFY2(EqualFiles(dataPath + "/" + transformFileName, dataPath + result), "The function output is different from the expected output!");
}

void TestCemrgCommandLine::ExecuteTracking_data() {
    QTest::addColumn<QString>("fixedFileName");
    QTest::addColumn<QString>("movingFileName");
    QTest::addColumn<QString>("transformFileName");
    QTest::addColumn<QString>("result");

    const array<QString[4], 1> trackingData { {
        {"/sphere_initial.nii", "/sphere_shifted.nii", "tracking_output.dof", "/tracking_expected.dof"}
    } };

    for (size_t i = 0; i < trackingData.size(); i++)
        QTest::newRow(("Test " + to_string(i + 1)).c_str()) << trackingData[i][0] << trackingData[i][1] << trackingData[i][2] << trackingData[i][3];
}

void TestCemrgCommandLine::ExecuteTracking() {
    QFETCH(QString, fixedFileName);
    QFETCH(QString, movingFileName);
    QFETCH(QString, transformFileName);
    QFETCH(QString, result);

    QFile::remove(dataPath + "/dcm-0.nii");
    QFile::remove(dataPath + "/dcm-1.nii");
    QFile::copy(dataPath + fixedFileName, dataPath + "/dcm-0.nii");
    QFile::copy(dataPath + movingFileName, dataPath + "/dcm-1.nii");
    cemrgCommandLine->ExecuteTracking(dataPath, "imgTimes.lst", QString(), transformFileName);
    QVERIFY2(EqualFiles(dataPath + "/" + transformFileName, dataPath + result), "The function output is different from the expected output!");
}

void TestCemrgCommandLine::ExecuteApplying_data() {
    QTest::addColumn<QString>("imageFileName");
    QTest::addColumn<double>("iniTime");
    QTest::addColumn<QString>("transformFileName");
    QTest::addColumn<int>("noFrames");
    QTest::addColumn<int>("smoothness");
    QTest::addColumn<QString>("result");

    const array<tuple<QString, double, QString, int, int, QString>, 1> applyingData { {
        {"surf_expected_1.vtk", 0, "tracking_expected.dof", 2, 1, "/applying_expected-"}
    } };

    for (size_t i = 0; i < applyingData.size(); i++)
        QTest::newRow(("Test " + to_string(i + 1)).c_str()) << get<0>(applyingData[i]) << get<1>(applyingData[i]) << get<2>(applyingData[i]) << get<3>(applyingData[i]) << get<4>(applyingData[i]) << get<5>(applyingData[i]);
}

void TestCemrgCommandLine::ExecuteApplying() {
    QFETCH(QString, imageFileName);
    QFETCH(double, iniTime);
    QFETCH(QString, transformFileName);
    QFETCH(int, noFrames);
    QFETCH(int, smoothness);
    QFETCH(QString, result);

    cemrgCommandLine->ExecuteApplying(dataPath, imageFileName, iniTime, transformFileName, noFrames, smoothness);
    for (int i = 0; i < noFrames * smoothness; i++)
        QVERIFY2(EqualFiles(dataPath + "/transformed-" + QString::number(i) + ".vtk", dataPath + result + QString::number(i) + ".vtk"), "The function output is different from the expected output!");
}

void TestCemrgCommandLine::ExecuteCreateCGALMesh_data() {
    QTest::addColumn<QString>("imageFileName");
    QTest::addColumn<QString>("paramsFileName");
    QTest::addColumn<QString>("outputName");

    const array<QString[3], 1> cgalMeshData { {
        {"sphere.inr", "/sphere.par", "sphere"}
    } };

    for (size_t i = 0; i < cgalMeshData.size(); i++)
        QTest::newRow(("Test " + to_string(i + 1)).c_str()) << cgalMeshData[i][0] << cgalMeshData[i][1] << cgalMeshData[i][2];
}

void TestCemrgCommandLine::ExecuteCreateCGALMesh() {
    QFETCH(QString, imageFileName);
    QFETCH(QString, paramsFileName);
    QFETCH(QString, outputName);

    // It creates a different file each time; therefore, we only check whether it exists
    QString cgalMeshOutput = cemrgCommandLine->ExecuteCreateCGALMesh(dataPath, outputName, dataPath + paramsFileName, imageFileName);
    QVERIFY(QFileInfo(cgalMeshOutput).exists());
}

int CemrgCommandLineTest(int argc, char *argv[]) {
    QApplication app(argc, argv);
    app.setAttribute(Qt::AA_Use96Dpi, true);
    QTEST_DISABLE_KEYPAD_NAVIGATION
    TestCemrgCommandLine tc;
    QTEST_SET_MAIN_SOURCE_PATH
    return QTest::qExec(&tc, argc, argv);
}
