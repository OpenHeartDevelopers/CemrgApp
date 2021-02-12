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
 * Atrial Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef CemrgAtrialTools_h
#define CemrgAtrialTools_h

#include <mitkSurface.h>
#include <mitkIOUtil.h>
#include <mitkImage.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkRegularPolygonSource.h>
#include <vmtk/vtkvmtkPolyDataCenterlines.h>
#include <vmtk/vtkvmtkPolyDataCenterlineSections.h>
#include <QString>

#include <itkPoint.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkThresholdImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryMorphologicalOpeningImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkImageRegionIterator.h>
#include <itkAddImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkLabelShapeKeepNObjectsImageFilter.h>

#include "CemrgAtriaClipper.h"

// The following header file is generated by CMake and thus it's located in
// the build directory. It provides an export macro for classes and functions
// that you want to be part of the public interface of your module.
#include <MitkCemrgAppModuleExports.h>

typedef itk::Image<uint16_t,3> ImageType;
typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholdType;
typedef itk::BinaryBallStructuringElement<ImageType::PixelType, 3> StrElType;
typedef itk::BinaryMorphologicalOpeningImageFilter<ImageType, ImageType, StrElType> ImFilterType;
typedef itk::ImageRegionIterator<ImageType> IteratorType;
typedef itk::AddImageFilter<ImageType, ImageType, ImageType> AddFilterType;

typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> ConnectedComponentImageFilterType;
typedef itk::LabelShapeKeepNObjectsImageFilter<ImageType> LabelShapeKeepNObjImgFilterType;

class MITKCEMRGAPPMODULE_EXPORT CemrgAtrialTools {

public:
    CemrgAtrialTools();
    void SetDefaultSegmentationTags();
    ImageType::Pointer LoadImage(QString imagePath);

    inline void SetDebugMode(bool s2d){debugSteps=s2d;};
    inline void SetDebugModeOn(){SetDebugMode(true);};
    inline void SetDebugModeOff(){SetDebugMode(false);};
    inline bool Debugging(){return debugSteps;};

    inline void SetWorkingDirectory(QString wd){directory =wd;};
    inline void SetTagSegmentationName(QString tsn){tagSegName = tsn;};

    // Label setters
    inline void SetMV(int _mv){mv=_mv;};
    inline void SetBODY(int b){abody=b;};
    inline void SetLAAP(int ap){laap=ap;};
    inline void SetLIPV(int li){lipv=li;};
    inline void SetLSPV(int ls){lspv=ls;};
    inline void SetRSPV(int rs){rspv=rs;};
    inline void SetRIPV(int ri){lspv=ri;};

    void SetNaiveSegmentationTags();

    ImageType::Pointer CleanAutomaticSegmentation(QString dir, QString segName="LA-cemrgnet.nii");
    ImageType::Pointer AssignAutomaticLabels(ImageType::Pointer im, QString dir, QString outName="tag-segmentation.nii");
    void GetSurfaceWithTags(ImageType::Pointer im, QString dir, QString outName, double th=0.5, double bl=0.8, double smth=3, double ds=0.5);
    void ClipMitralValveAuto(QString dir, QString mvName, QString outName);

    // helper functions
    ImageType::Pointer ExtractLabel(QString tag, ImageType::Pointer im, uint16_t label, uint16_t filterRadius=1.0);
    ImageType::Pointer AddImage(ImageType::Pointer im1, ImageType::Pointer im2);
    ThresholdType::Pointer ThresholdImage(ImageType::Pointer input, uint16_t thresholdVal);
    ImFilterType::Pointer ImOpen(ImageType::Pointer input, uint16_t radius);
    mitk::Surface::Pointer FlipSegmentation(mitk::Surface::Pointer surf);
    void SaveImageToDisk(ImageType::Pointer im, QString dir, QString imName);

private:

    mitk::Surface::Pointer surface;
    ImageType::Pointer atriumSegmentation;
    ImageType::Pointer body;
    ImageType::Pointer veins;
    ImageType::Pointer mitralvalve;
    std::unique_ptr<CemrgAtriaClipper> clipper;

    QString directory, tagSegName;

    int abody, mv, laap, lipv, lspv, rspv, ripv;
    int autolaap, autolipv, autolspv, autorspv, autoripv;
    std::vector<int> detectedLabels;
    bool segmentationSet, debugSteps;
    //Constant Vein Labels
    //const int LEFTSUPERIORPV  = 11;
    //const int LEFTMIDDLEPV    = 12;
    //const int LEFTINFERIORPV  = 13;
    //const int LEFTCOMMONPV    = 14;
    //const int RIGHTSUPERIORPV = 15;
    //const int RIGHTMIDDLEPV   = 16;
    //const int RIGHTINFERIORPV = 17;
    //const int RIGHTCOMMONPV   = 18;
    //const int APPENDAGECUT    = 19;
    //const int APPENDAGEUNCUT  = 20;
    //const int DEFAULTVALUE    = 21;
};

#endif // CemrgAtrialTools_h
