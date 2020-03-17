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
 * Atrial Scar Tools for MITK
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
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
#include <MyCemrgLibExports.h>


class MyCemrgLib_EXPORT CemrgScar3D {
    
public:
    
    CemrgScar3D();

    mitk::Surface::Pointer ClipMesh3D(mitk::Surface::Pointer surface, mitk::PointSet::Pointer landmarks);
    mitk::Surface::Pointer Scar3D(std::string directory, mitk::Image::Pointer lgeImage);
    bool CalculateMeanStd(mitk::Image::Pointer lgeImage, mitk::Image::Pointer roiImage, double& mean, double& stdv);
    double Thresholding(double thresh);
    
    double GetMinScalar() const;
    double GetMaxScalar() const;
    void SetMinStep(int value);
    void SetMaxStep(int value);
    void SetMethodType(int value);
    void SetScarSegImage(const mitk::Image::Pointer image);    
    
private:
    
    int methodType;
    int minStep, maxStep;
    double minScalar, maxScalar;
    vtkSmartPointer<vtkFloatArray> scalars;
    typedef itk::Image<short,3> itkImageType;
    itkImageType::Pointer scarSegImage;
    
    double GetIntensityAlongNormal(
            itkImageType::Pointer scarImage, itkImageType::Pointer visitedImage,
            double n_x, double n_y, double n_z, double centre_x, double centre_y, double centre_z);
    double GetStatisticalMeasure(
            std::vector<mitk::Point3D> pointsOnAndAroundNormal,
            itkImageType::Pointer scarImage, itkImageType::Pointer visitedImage, int measure);
    void ItkDeepCopy(itkImageType::Pointer input, itkImageType::Pointer output);
    mitk::Surface::Pointer ReadVTKMesh(std::string meshPath);
};

#endif // CemrgScar3D_h
