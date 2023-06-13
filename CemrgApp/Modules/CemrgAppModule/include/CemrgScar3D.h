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
 * Scar Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef CemrgScar3D_h
#define CemrgScar3D_h

// Qmitk
#include <mitkImage.h>
#include <mitkSurface.h>
#include <mitkPointSet.h>
#include <vtkFloatArray.h>
#include <MitkCemrgAppModuleExports.h>
#include <QString>
#include <QDir>

class MITKCEMRGAPPMODULE_EXPORT CemrgScar3D {

public:
    enum DebugImageType { BACKGROUND=0, CORRIDOR, PROJECTED };

    CemrgScar3D();
    mitk::Surface::Pointer Scar3D(std::string directory, mitk::Image::Pointer lgeImage, std::string segname = "segmentation.vtk");

    mitk::Surface::Pointer ClipMesh3D(mitk::Surface::Pointer surface, mitk::PointSet::Pointer landmarks);
    bool CalculateMeanStd(mitk::Image::Pointer lgeImage, mitk::Image::Pointer roiImage, double& mean, double& stdv);
    double Thresholding(double thresh);
    void SaveScarDebugImage(QString name, QString dir);
    void SaveNormalisedScalars(double divisor, mitk::Surface::Pointer surface, QString name);
    void PrintThresholdingResults(QString dir, std::vector<double> values_vector, int threshType, double mean, double stdv, bool printGuide = true);
    void PrintSingleThresholdingResult(QString dir, double value, int threshType, double mean, double stdv);

    inline double GetMinScalar() const { return minScalar; };
    inline double GetMaxScalar() const { return maxScalar; };
    inline void SetMinStep(int value) { minStep = value; };
    inline void SetMaxStep(int value) { maxStep = value; };
    inline void SetMethodType(int value) { methodType = value; };
    void SetScarSegImage(const mitk::Image::Pointer image);

    inline void SetVoxelBasedProjection(bool value){voxelBasedProjection=value;};
    inline void SetRoiLegacyNormals(bool value){roiLegacyNormals=value;};
    inline void SetRoiRadiusOption(bool value){roiRadiusOption=value;};

private:

    int methodType;
    int minStep, maxStep;
    bool voxelBasedProjection, roiLegacyNormals, roiRadiusOption;
    double minScalar, maxScalar;
    vtkSmartPointer<vtkFloatArray> scalars;

    typedef itk::Image<short, 3> itkImageType;
    itkImageType::Pointer scarSegImage;
    itk::Image<short, 3>::Pointer scarDebugCorridor;

    double GetIntensityAlongNormal(
        itkImageType::Pointer scarImage, itkImageType::Pointer visitedImage,
        double n_x, double n_y, double n_z, double centre_x, double centre_y, double centre_z);
    double GetStatisticalMeasure(
        std::vector<mitk::Point3D> pointsOnAndAroundNormal,
        itkImageType::Pointer scarImage, itkImageType::Pointer visitedImage, int measure);
    void ItkDeepCopy(itkImageType::Pointer input, itkImageType::Pointer output);
};

#endif // CemrgScar3D_h
