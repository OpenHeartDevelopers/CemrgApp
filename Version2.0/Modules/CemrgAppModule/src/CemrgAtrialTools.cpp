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

//Qmitk
#include <mitkProgressBar.h>
#include <mitkIOUtil.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkMorphologicalOperations.h>

#include "CemrgAtrialTools.h"

//VTK
#include <vtkPointLocator.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkClipPolyData.h>
#include <vtkImplicitBoolean.h>
#include <vtkCenterOfMass.h>
#include <vtkPlane.h>
#include <vtkSphere.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkLinearExtrusionFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolygon.h>
#include <vtkCellArray.h>

//ITK
#include <itkSubtractImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkLabelShapeKeepNObjectsImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkImageDuplicator.h>
#include <itkResampleImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

//Qt
#include <QtDebug>
#include <QString>
#include <QFileInfo>
#include <QProcess>
#include <QMessageBox>
#include <numeric>

//CemrgAppModule
#include "CemrgScar3D.h"
#include "CemrgCommonUtils.h"
#include "CemrgMeasure.h"


CemrgAtrialTools::CemrgAtrialTools() {
    this->debugSteps = true;
    this->atriumSegmentation = ImageType::New();
    this->tagSegName = "tag-segmentation.nii";

    SetDefaultSegmentationTags();
}

void CemrgAtrialTools::SetDefaultSegmentationTags(){
    abody = 1;
    mv = 10;
    laap = -1;
    lspv = -1;
    lipv = -1;
    rspv = -1;
    ripv = -1;

    autolaap=-1;
    autolspv=-1;
    autolipv=-1;
    autorspv=-1;
    autoripv=-1;

}


ImageType::Pointer CemrgAtrialTools::LoadImage(QString imagePath){
    mitk::Image::Pointer mitkImage = mitk::IOUtil::Load<mitk::Image>(imagePath.toStdString());
    ImageType::Pointer itkImage = ImageType::New();
    mitk::CastToItkImage(mitkImage, itkImage);

    return itkImage;
}

ImageType::Pointer CemrgAtrialTools::CleanAutomaticSegmentation(QString dir, QString segName){
    QString inputPath = dir + mitk::IOUtil::GetDirectorySeparator() + segName;
    ImageType::Pointer orgSegImage = LoadImage(inputPath);

    MITK_INFO << "Extracting clean segmentation.";
    ConnectedComponentImageFilterType::Pointer conn1 = ConnectedComponentImageFilterType::New();
    conn1->SetInput(orgSegImage);
    conn1->Update();

    LabelShapeKeepNObjImgFilterType::Pointer keepObjs = LabelShapeKeepNObjImgFilterType::New();
    keepObjs->SetInput(conn1->GetOutput());
    keepObjs->SetBackgroundValue(0);
    keepObjs->SetNumberOfObjects(1);
    keepObjs->SetAttribute(LabelShapeKeepNObjImgFilterType::LabelObjectType::NUMBER_OF_PIXELS);
    keepObjs->Update();

    ImageType::Pointer atriumCoarse = keepObjs->GetOutput();

    SaveImageToDisk(atriumCoarse, dir, "1_CleanSegmentation.nii");

    MITK_INFO << "Relabel clean segmentation's veins";
    IteratorType segImgIter(atriumCoarse, atriumCoarse->GetLargestPossibleRegion());
    IteratorType ogImgIter(orgSegImage, orgSegImage->GetLargestPossibleRegion());
    while(!segImgIter.IsAtEnd()){
        if(segImgIter.Get() > 0){
            if(ogImgIter.Get() == 3){ // mitral valve
                segImgIter.Set(mv);
            } else{
                segImgIter.Set((int)ogImgIter.Get());
            }
        }
        ++segImgIter;
        ++ogImgIter;
    }

    MITK_INFO << "Extracting and cleaning atrium body";
    uint16_t bodythresh = 1;
    ImageType::Pointer body = ExtractLabel("LA-body", atriumCoarse, bodythresh, 1.0);

    MITK_INFO << "Extracting and cleaning pulmonary veins";
    uint16_t veinsthresh = 2;
    ImageType::Pointer veins = ExtractLabel("Veins", atriumCoarse, veinsthresh, 2.0);

    MITK_INFO << "Extracting and cleaning Mitral Valve";
    uint16_t mvthresh = mv;
    ImageType::Pointer mitralvalve = ExtractLabel("MitralValve", atriumCoarse, mvthresh, 1.0);

    ConnectedComponentImageFilterType::Pointer conn2 = ConnectedComponentImageFilterType::New();
    conn2->SetInput(mitralvalve);
    conn2->Update();

    LabelShapeKeepNObjImgFilterType::Pointer keepObjs2 = LabelShapeKeepNObjImgFilterType::New();
    keepObjs2->SetInput(conn2->GetOutput());
    keepObjs2->SetBackgroundValue(0);
    keepObjs2->SetNumberOfObjects(1);
    keepObjs2->SetAttribute(LabelShapeKeepNObjImgFilterType::LabelObjectType::NUMBER_OF_PIXELS);
    keepObjs2->Update();

    IteratorType bodyIter(body, body->GetLargestPossibleRegion());
    IteratorType veinsIter(veins, veins->GetLargestPossibleRegion());
    IteratorType mvIter(keepObjs2->GetOutput(), keepObjs2->GetOutput()->GetLargestPossibleRegion());
    while (!bodyIter.IsAtEnd()) {
        if(veinsIter.Get() > 0){
            bodyIter.Set(2);
        }
        if (mvIter.Get() > 0) {
            bodyIter.Set(abody);
        }
        ++bodyIter;
        ++veinsIter;
        ++mvIter;
    }

    atriumSegmentation = body;
    SaveImageToDisk(atriumSegmentation, dir, "2_SegmentationWithVeins.nii");
    SaveImageToDisk(keepObjs2->GetOutput(), dir, "prodMVI.nii");

    return atriumSegmentation;
}

ImageType::Pointer CemrgAtrialTools::AssignAutomaticLabels(ImageType::Pointer im, QString dir, QString outName){

    MITK_INFO << "Extracting and cleaning pulmonary veins";

    uint16_t veinsthresh = 2;
    ImageType::Pointer veins = ExtractLabel("Veins", im, veinsthresh, 2.0);
    SaveImageToDisk(veins, dir, "3_ThresholdedVeins.nii");

    MITK_INFO << "Label each detected vein";
    ConnectedComponentImageFilterType::Pointer conn2 = ConnectedComponentImageFilterType::New();
    conn2->SetInput(veins);
    conn2->Update();

    SaveImageToDisk(conn2->GetOutput(), dir, "4_ConnectedComponentsVeins.nii");

    IteratorType imIter(im, im->GetLargestPossibleRegion());
    IteratorType veinsIter(conn2->GetOutput(), conn2->GetOutput()->GetLargestPossibleRegion());

    int idx=0;
    while(!veinsIter.IsAtEnd()){
        if(veinsIter.Get() > 0){
            detectedLabels.push_back((int)veinsIter.Get()+10);
            imIter.Set(detectedLabels[idx]);
            idx++;
        }
        ++imIter;
        ++veinsIter;
    }
    MITK_INFO << "Relabelling finished";

    SetNaiveSegmentationTags();
    SetTagSegmentationName(outName);
    SaveImageToDisk(im, dir, tagSegName);

    segmentationSet = true;
    atriumSegmentation = im;

    return im;
}

void CemrgAtrialTools::GetSurfaceWithTags(ImageType::Pointer im, QString dir, QString outName, double th, double bl, double smth, double ds){
    MITK_INFO << "Extracting Surface";
    mitk::Image::Pointer veinsRelabeledImg = mitk::Image::New();
    mitk::CastToMitkImage(im, veinsRelabeledImg);
    QString path = dir + mitk::IOUtil::GetDirectorySeparator();

    mitk::Surface::Pointer segSurface = CemrgCommonUtils::ExtractSurfaceFromSegmentation(veinsRelabeledImg, th, bl, smth, ds);
    CemrgCommonUtils::FlipXYPlane(segSurface, dir, "segmentation.vtk");

    std::unique_ptr<CemrgScar3D> scar(new CemrgScar3D());
    scar->SetMinStep(-1);
    scar->SetMaxStep(3);
    scar->SetMethodType(2);
    scar->SetScarSegImage(veinsRelabeledImg);

    MITK_INFO << "Projection of image labels onto surface";
    surface = scar->Scar3D(dir.toStdString(), veinsRelabeledImg);

    QString outputPath = dir + mitk::IOUtil::GetDirectorySeparator() + outName;
    mitk::IOUtil::Save(surface, outputPath.toStdString());
    MITK_INFO << ("Saved output shell to " + outputPath).toStdString();

}

void CemrgAtrialTools::ClipMitralValveAuto(QString dir, QString mvName, QString outName){
    MITK_INFO << "[ClipMitralValveAuto]";
    // Make vtk of prodMVI
    MITK_INFO << "Loading Mitral Valve Image (mvi)";
    QString mviPath = dir + mitk::IOUtil::GetDirectorySeparator() + mvName;
    mitk::Image::Pointer mvi = mitk::IOUtil::Load<mitk::Image>(mviPath.toStdString());
    mitk::MorphologicalOperations::Dilate(mvi, 2, mitk::MorphologicalOperations::StructuralElementType::Ball);

    MITK_INFO << "Extract surface from mvi";
    mitk::Surface::Pointer clipperSurf = CemrgCommonUtils::ExtractSurfaceFromSegmentation(mvi, 0.5, 0, 10);

    if(debugSteps){
        QFileInfo fi(mviPath);
        QString mvSurfPath = dir + mitk::IOUtil::GetDirectorySeparator() + fi.baseName() + ".vtk";
        MITK_INFO << ("[debug] Saving surface to " + mvSurfPath).toStdString();
        mitk::IOUtil::Save(clipperSurf, mvSurfPath.toStdString());
    }

    // Implement code from command line tool
    MITK_INFO << "Define implicit function object (mitral valve)";
    vtkSmartPointer<vtkImplicitPolyDataDistance> implicitFn = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    implicitFn->SetInput(clipperSurf->GetVtkPolyData());
    vtkMTimeType mtime = implicitFn->GetMTime();
    std::cout << "MTime: " << mtime<< std::endl ;

    MITK_INFO << "Define clipper object";
    vtkSmartPointer<vtkClipPolyData> mvclipper = vtkSmartPointer<vtkClipPolyData>::New();
    mvclipper->SetClipFunction(implicitFn);
    mvclipper->SetInputData(surface->GetVtkPolyData());
    mvclipper->InsideOutOn();
    // mvclipper->InsideOutOff();
    mvclipper->Update();

    MITK_INFO << "Extract and clean surface mesh.";
    vtkSmartPointer<vtkDataSetSurfaceFilter> surfer = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfer->SetInputData(mvclipper->GetOutput());
    surfer->Update();

    MITK_INFO << "Cleaning...";
    vtkSmartPointer<vtkCleanPolyData> clean = vtkSmartPointer<vtkCleanPolyData>::New();
    clean->SetInputConnection(surfer->GetOutputPort());
    clean->Update();

    MITK_INFO << "Largest region...";
    vtkSmartPointer<vtkPolyDataConnectivityFilter> lrgRegion = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    lrgRegion->SetInputConnection(clean->GetOutputPort());
    lrgRegion->SetExtractionModeToLargestRegion();
    lrgRegion->Update();
    clean = vtkSmartPointer<vtkCleanPolyData>::New();
    clean->SetInputConnection(lrgRegion->GetOutputPort());
    clean->Update();

    surface->SetVtkPolyData(clean->GetOutput());

    if(!outName.isEmpty()){
        QString outPath = dir + mitk::IOUtil::GetDirectorySeparator() + outName;
        MITK_INFO << ("Saving to file: " + outPath).toStdString();
        mitk::IOUtil::Save(surface, outPath.toStdString());
    }
}

//helper functions
ImageType::Pointer CemrgAtrialTools::ExtractLabel(QString tag, ImageType::Pointer im, uint16_t label, uint16_t filterRadius){
    MITK_INFO << ("Thresholding " + tag + " from clean segmentation").toStdString();
    ThresholdType::Pointer thresVeins = ThresholdImage(im, label);

    if(filterRadius > 0){
        // morphological opening step
        MITK_INFO << "Morphological opening - remove speckle noise";
        ImFilterType::Pointer imopen = ImOpen(thresVeins->GetOutput(), filterRadius);
        return imopen->GetOutput();

    } else{
        MITK_INFO << "Morphological opening - ignored";
        return thresVeins->GetOutput();
    }

}

ImageType::Pointer CemrgAtrialTools::AddImage(ImageType::Pointer im1, ImageType::Pointer im2){
    AddFilterType::Pointer sum = AddFilterType::New();
    sum->SetInput1(im1);
    sum->SetInput2(im2);
    sum->Update();

    return sum->GetOutput();
}

ThresholdType::Pointer CemrgAtrialTools::ThresholdImage(ImageType::Pointer input, uint16_t thresholdVal){
    ThresholdType::Pointer thresholdOutput = ThresholdType::New();
    thresholdOutput->SetInput(input);
    thresholdOutput->SetLowerThreshold(thresholdVal);
    thresholdOutput->SetUpperThreshold(thresholdVal);
    thresholdOutput->SetInsideValue(1);
    thresholdOutput->SetOutsideValue(0);
    thresholdOutput->Update();

    return thresholdOutput;
}

ImFilterType::Pointer CemrgAtrialTools::ImOpen(ImageType::Pointer input, uint16_t radius){
    StrElType structuringElement;
    structuringElement.SetRadius(static_cast<unsigned long>(radius));
    structuringElement.CreateStructuringElement();

    ImFilterType::Pointer imOpenFilter = ImFilterType::New();
    imOpenFilter->SetInput(input);
    imOpenFilter->SetKernel(structuringElement);
    imOpenFilter->Update();

    return imOpenFilter;
}

void CemrgAtrialTools::SaveImageToDisk(ImageType::Pointer im, QString dir, QString imName){
    if(debugSteps || tagSegName.compare(imName)==0){
        QString outputPath = dir + mitk::IOUtil::GetDirectorySeparator() + imName;
        MITK_INFO << ("Saving image " + imName + " to: " + dir).toStdString();
        mitk::Image::Pointer outputImg = mitk::Image::New();
        mitk::CastToMitkImage(im, outputImg);
        mitk::IOUtil::Save(outputImg, outputPath.toStdString());
    } else {
        MITK_INFO << "Saving to disk disabled.";
    }
}

// helper functions
void CemrgAtrialTools::SetNaiveSegmentationTags(){
    int numDetectedLabels = detectedLabels.size();
    if(numDetectedLabels>=5){
        autolaap=detectedLabels[0];
        autolipv=detectedLabels[1];
        autolspv=detectedLabels[2];
        autorspv=detectedLabels[3];
        autoripv=detectedLabels[4];
    }
}
