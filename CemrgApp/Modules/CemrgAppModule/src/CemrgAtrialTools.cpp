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
#include <vtkConnectivityFilter.h>
#include <vtkImageResize.h>
#include <vtkImageChangeInformation.h>
#include <vtkGeometryFilter.h>
#include <vtkMassProperties.h>
#include <vtkCellLocator.h>

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
#include "CemrgScarAdvanced.h"
#include "CemrgCommonUtils.h"
#include "CemrgMeasure.h"


CemrgAtrialTools::CemrgAtrialTools() {
    this->debugSteps = true;
    this->debugMessages = false;
    this->surfLoaded = false;
    this->atriumSegmentation = ImageType::New();
    this->tagSegName = "Labelled.nii";

    SetDefaultSegmentationTags();
}

void CemrgAtrialTools::SetDefaultSegmentationTags(){
    abody = 1;
    mv = 10;
    lspv = 11;
    lipv = 13;
    rspv = 15;
    ripv = 17;
    laap = 19;

    naivelaap=-1;
    naivelspv=-1;
    naivelipv=-1;
    naiverspv=-1;
    naiveripv=-1;

}

ShortImageType::Pointer LoadShortImage(QString imagePath){
    mitk::Image::Pointer mitkImage = mitk::IOUtil::Load<mitk::Image>(imagePath.toStdString());
    ShortImageType::Pointer itkImage = ShortImageType::New();
    mitk::CastToItkImage(mitkImage, itkImage);

    return itkImage;
}

ImageType::Pointer CemrgAtrialTools::CemrgLoadImage(QString imagePath, bool binarise){
    mitk::Image::Pointer mitkImage = mitk::IOUtil::Load<mitk::Image>(imagePath.toStdString());
    ImageType::Pointer itkImage = ImageType::New();
    mitk::CastToItkImage(mitkImage, itkImage);

    if(binarise){
        using IteratorType = itk::ImageRegionIteratorWithIndex<ImageType>;
        IteratorType imIter(itkImage, itkImage->GetLargestPossibleRegion());

        imIter.GoToBegin();
        while(!imIter.IsAtEnd()){
            if(imIter.Get()>0){
                imIter.Set(1);
            } else{
                imIter.Set(0);
            }
            ++imIter;
        }
    }

    return itkImage;
}

void CemrgAtrialTools::AdjustSegmentationLabelToImage(QString segImPath, QString imPath, QString outImPath){
    mitk::Image::Pointer im = mitk::IOUtil::Load<mitk::Image>(imPath.toStdString());
    mitk::Image::Pointer segIm = mitk::IOUtil::Load<mitk::Image>(segImPath.toStdString());
    double origin[3]; double spacing[3];
    im->GetGeometry()->GetOrigin().ToArray(origin);
    im->GetGeometry()->GetSpacing().ToArray(spacing);

    vtkSmartPointer<vtkImageResize> resizeFilter = vtkSmartPointer<vtkImageResize>::New();
    resizeFilter->SetResizeMethodToOutputDimensions();
    resizeFilter->SetOutputDimensions(im->GetDimension(0), im->GetDimension(1), im->GetDimension(2));
    resizeFilter->InterpolateOff();
    resizeFilter->SetInputData(segIm->GetVtkImageData());
    resizeFilter->Update();

    vtkSmartPointer<vtkImageChangeInformation> changeFilter = vtkSmartPointer<vtkImageChangeInformation>::New();
    changeFilter->SetInputConnection(resizeFilter->GetOutputPort());
    changeFilter->SetOutputSpacing(spacing);
    changeFilter->SetOutputOrigin(origin);
    changeFilter->Update();

    segIm->Initialize(changeFilter->GetOutput());
    segIm->SetVolume(changeFilter->GetOutput()->GetScalarPointer());

    QString outputPath = (outImPath.isEmpty()) ? segImPath : outImPath;
    mitk::IOUtil::Save(segIm, outputPath.toStdString());
}

void CemrgAtrialTools::ResampleSegmentationLabelToImage(QString segImPath, QString imPath, QString outImPath){
    MITK_INFO(debugMessages) << ("[ResampleSegmentationLabelToImage] Loading " + segImPath).toStdString();
    ImageType::Pointer segItk = CemrgLoadImage(segImPath);
    MITK_INFO(debugMessages) << ("[ResampleSegmentationLabelToImage] Loading " + imPath).toStdString();
    ImageType::Pointer im = CemrgLoadImage(imPath);

    itk::ResampleImageFilter<ImageType, ImageType>::Pointer resampleFilter;
    resampleFilter = itk::ResampleImageFilter<ImageType, ImageType>::New();
    resampleFilter->SetInput(segItk);
    resampleFilter->SetReferenceImage(im);
    resampleFilter->SetUseReferenceImage(true);
    resampleFilter->SetInterpolator(itk::NearestNeighborInterpolateImageFunction<ImageType>::New());
    resampleFilter->SetDefaultPixelValue(0);
    resampleFilter->UpdateLargestPossibleRegion();
    MITK_INFO(debugMessages) << "[ResampleSegmentationLabelToImage] finished resampling";

    QString outputPath = (outImPath.isEmpty()) ? segImPath : outImPath;
    QFileInfo fo(outputPath);

    SaveImageToDisk(resampleFilter->GetOutput(), fo.absolutePath(), fo.baseName());
}

ImageType::Pointer CemrgAtrialTools::RemoveNoiseFromAutomaticSegmentation(QString dir, QString segName){
    QString inputPath = dir + "/" + segName;
    ImageType::Pointer orgSegImage = CemrgLoadImage(inputPath);

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

    ImageType::Pointer outIm = keepObjs->GetOutput();
    IteratorType imIter(outIm, outIm->GetLargestPossibleRegion());

    imIter.GoToBegin();
    while(!imIter.IsAtEnd()){
        float value = (imIter.Get() > 0) ? 1 : 0;
        imIter.Set(value);

        ++imIter;
    }

    return outIm;
}

void CemrgAtrialTools::QuickBinarise(ImageType::Pointer imToBin){
    // binarise to 1, 0
    IteratorType imIter(imToBin, imToBin->GetLargestPossibleRegion());

    imIter.GoToBegin();
    while(!imIter.IsAtEnd()){
        float value = (imIter.Get() > 0) ? 1 : 0;
        imIter.Set(value);

        ++imIter;
    }
}

ImageType::Pointer CemrgAtrialTools::CleanAutomaticSegmentation(QString dir, QString segName, QString cleanName){
    QString inputPath = dir + "/" + segName;
    ImageType::Pointer orgSegImage = CemrgLoadImage(inputPath);
    ImageType::Pointer atriumCoarse = RemoveNoiseFromAutomaticSegmentation(dir, segName);

    if(!cleanName.isEmpty()){
        SaveImageToDisk(atriumCoarse, dir, cleanName);
    }

    MITK_INFO << "[CleanAutomaticSegmentation] Relabel clean segmentation's veins";
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

    if(debugSteps){
        SaveImageToDisk(atriumCoarse, dir, "1_RemoveNoiseWithVeins.nii");
    }

    MITK_INFO << "[CleanAutomaticSegmentation] Extracting and cleaning atrium body";
    uint16_t bodythresh = 1;
    ImageType::Pointer body = ExtractLabel("LA-body", atriumCoarse, bodythresh, 1.0);

    MITK_INFO << "[CleanAutomaticSegmentation] Extracting and cleaning pulmonary veins";
    uint16_t veinsthresh = 2;
    ImageType::Pointer veins = ExtractLabel("Veins", atriumCoarse, veinsthresh, 3.0, 5);

    MITK_INFO << "[CleanAutomaticSegmentation] Extracting and cleaning Mitral Valve";
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
    MITK_INFO << "[CleanAutomaticSegmentation] Saving Mitral Valve";
    if(debugSteps){
        SaveImageToDisk(atriumSegmentation, dir, "2_SegmentationWithVeins.nii");
    }
    SaveImageToDisk(keepObjs2->GetOutput(), dir, "prodMVI.nii");

    return atriumSegmentation;
}

ImageType::Pointer CemrgAtrialTools::AssignAutomaticLabels(ImageType::Pointer im, QString dir, QString outName, bool relabel){

    MITK_INFO << "Extracting and cleaning pulmonary veins";

    uint16_t veinsthresh = 2;
    ImageType::Pointer veins = ExtractLabel("Veins", im, veinsthresh, 3.0, 5);

    MITK_INFO << "Label each detected vein";
    ConnectedComponentImageFilterType::Pointer conn2 = ConnectedComponentImageFilterType::New();
    conn2->SetInput(veins);
    conn2->Update();

    if(debugSteps){
        SaveImageToDisk(veins, dir, "3_ThresholdedVeins.nii");
        SaveImageToDisk(conn2->GetOutput(), dir, "4_ConnectedComponentsVeins.nii");
    }

    IteratorType imIter(im, im->GetLargestPossibleRegion());
    IteratorType veinsIter(conn2->GetOutput(), conn2->GetOutput()->GetLargestPossibleRegion());

    int idx=0;
    int relabelValue = relabel ? 30 : 1;
    while(!veinsIter.IsAtEnd()){
        if(veinsIter.Get() > 0){
            // create large even numbers 32, 34, 36, ... for detected labels that don't overlap with default labels (1,11,13, ...,19)
            int value = ((int)veinsIter.Get())*2+relabelValue;
            detectedLabels.push_back(value);
            imIter.Set(detectedLabels[idx]);
            idx++;
        }
        ++imIter;
        ++veinsIter;
    }
    MITK_INFO << "Relabelling finished";

    SetTagSegmentationName(outName);
    SaveImageToDisk(im, dir, tagSegName);

    segmentationSet = true;
    atriumSegmentation = im;

    return im;
}

ImageType::Pointer CemrgAtrialTools::AssignOstiaLabelsToVeins(ImageType::Pointer im, QString dir, QString outName){
    ImageType::Pointer veins = ThresholdImage(im, 2, 1000); // im >= 2
    ConnectedComponentImageFilterType::Pointer connectedVeins = ConnectedComponentImageFilterType::New();
    connectedVeins->SetInput(veins);
    connectedVeins->Update();

    RelabelFilterType::Pointer labelVeins = RelabelFilterType::New();
    labelVeins->SetInput(connectedVeins->GetOutput());
    labelVeins->Update();

    unsigned long numObjects = labelVeins->GetNumberOfObjects();

    IteratorType liter(labelVeins->GetOutput(), labelVeins->GetOutput()->GetLargestPossibleRegion());
    IteratorType pviter(im, im->GetLargestPossibleRegion());

    liter.GoToBegin();
    pviter.GoToBegin();
    std::vector<int> veinTips(numObjects, 0);
    std::vector<int> veinOstia(numObjects, 0);
    while(!liter.IsAtEnd()){
        int l = (int) liter.Get();
        int pv = (int) pviter.Get();
        if(l>0){
            if(IsLabel(pv)>0){
                veinOstia.at(l-1) = pv;
            } else{
                veinTips.at(l-1) = pv;
                pviter.Set(veinOstia.at(l-1));
            }
        }
        ++liter;
        ++pviter;
    }

    SetTagSegmentationName(outName);
    SaveImageToDisk(im, dir, tagSegName);

    segmentationSet = true;
    atriumSegmentation = im;

    return im;
}

mitk::Image::Pointer CemrgAtrialTools::SurfSegmentation(ImageType::Pointer im, QString dir, QString outName, double th, double bl, double smth, double ds){
    MITK_INFO << "Extracting Surface";
    mitk::Image::Pointer labelSegIm = mitk::Image::New();
    mitk::CastToMitkImage(im, labelSegIm);

    mitk::Surface::Pointer segSurface = CemrgCommonUtils::ExtractSurfaceFromSegmentation(labelSegIm, th, bl, smth, ds);
    CemrgCommonUtils::FlipXYPlane(segSurface, dir, outName);

    return labelSegIm;
}

void CemrgAtrialTools::ProjectTagsOnSurface(ImageType::Pointer im, QString dir, QString outName, double th, double bl, double smth, double ds, bool createSurface){

    mitk::Image::Pointer labelSegIm = mitk::Image::New();
    if(createSurface){
        MITK_INFO << "Extracting Surface";
        labelSegIm = SurfSegmentation(im, dir, "segmentation.vtk", th, bl, smth, ds);
    } else{
        mitk::CastToMitkImage(im, labelSegIm);
    }

    std::unique_ptr<CemrgScar3D> scar(new CemrgScar3D());
    scar->SetMinStep(-1);
    scar->SetMaxStep(3);
    scar->SetMethodType(2);
    scar->SetScarSegImage(labelSegIm);

    MITK_INFO << "Projection of image labels onto surface";
    surface = scar->Scar3D(dir.toStdString(), labelSegIm);
    surface->GetVtkPolyData()->GetCellData()->GetScalars()->SetName("elemTag");
    surfLoaded=true;

    QString outputPath = dir + "/" + outName;
    mitk::IOUtil::Save(surface, outputPath.toStdString());
    MITK_INFO << ("Saved output shell to " + outputPath).toStdString();
}

void CemrgAtrialTools::ProjectTagsOnExistingSurface(ImageType::Pointer im, QString dir, QString outName, QString existingSurfaceName){
    mitk::Image::Pointer labelSegIm = mitk::Image::New();
    mitk::CastToMitkImage(im, labelSegIm);


    std::unique_ptr<CemrgScar3D> scar(new CemrgScar3D());
    scar->SetMinStep(-1);
    scar->SetMaxStep(3);
    scar->SetMethodType(4); // 2=max, 4=mode
    scar->SetScarSegImage(labelSegIm);

    MITK_INFO << "Projection of image labels onto surface";
    surface = scar->Scar3D(dir.toStdString(), labelSegIm, existingSurfaceName.toStdString());
    surface->GetVtkPolyData()->GetCellData()->GetScalars()->SetName("elemTag");
    surfLoaded=true;

    QString outputPath = dir + "/" + outName;
    mitk::IOUtil::Save(surface, outputPath.toStdString());
    MITK_INFO << ("Saved output shell to " + outputPath).toStdString();
}

void CemrgAtrialTools::ClipMitralValveAuto(QString dir, QString mvNameExt, QString outName, bool insideout){
    MITK_INFO << "[ClipMitralValveAuto]";
    // Make vtk of prodMVI
    MITK_INFO << "Loading Mitral Valve Image (mvi)";
    QString mviPath = dir + "/" + mvNameExt;
    mitk::Image::Pointer mvi = mitk::IOUtil::Load<mitk::Image>(mviPath.toStdString());
    mitk::MorphologicalOperations::Dilate(mvi, 2, mitk::MorphologicalOperations::StructuralElementType::Ball);

    MITK_INFO << "Extract surface from mvi";
    mitk::Surface::Pointer clipperSurf = CemrgCommonUtils::ExtractSurfaceFromSegmentation(mvi, 0.5, 0, 10);

    QFileInfo fi(mviPath);
    QString mvSurfPath = dir + "/" + fi.baseName() + ".vtk";
    mitk::IOUtil::Save(clipperSurf, mvSurfPath.toStdString());

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
    if(insideout){
        mvclipper->InsideOutOn();
    } else{
        mvclipper->InsideOutOff();
    }
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
        QString outPath = dir + "/" + outName;
        MITK_INFO << ("Saving to file: " + outPath).toStdString();
        mitk::IOUtil::Save(surface, outPath.toStdString());
    }
}

void CemrgAtrialTools::ProjectShellScalars(QString dir, QString scalarsShellPath, QString outputShellPath){
    QString prodPathOut = dir + "/";
    QFileInfo fi(outputShellPath);

    mitk::Surface::Pointer _scalarsShell = mitk::IOUtil::Load<mitk::Surface>(scalarsShellPath.toStdString());
    mitk::Surface::Pointer _outputShell = mitk::IOUtil::Load<mitk::Surface>(outputShellPath.toStdString());

    std::unique_ptr<CemrgScarAdvanced> csadv = std::unique_ptr<CemrgScarAdvanced>(new CemrgScarAdvanced());
    csadv->SetOutputPath(prodPathOut.toStdString());
    csadv->SetWeightedCorridorBool(false);
    csadv->SetSourceAndTarget(_outputShell->GetVtkPolyData(), _scalarsShell->GetVtkPolyData());
    csadv->TransformSource2Target(fi.baseName());
}

void CemrgAtrialTools::SetSurfaceLabels(QString correctLabels, QString naiveLabels){
    std::ifstream fileCorrect, fileNaive;
    fileCorrect.open(correctLabels.toStdString());
    fileNaive.open(naiveLabels.toStdString());

    std::string line, nline;

    while(std::getline(fileCorrect, line)){
        std::getline(fileNaive, nline);

        QString qline =  QString::fromStdString(line);
        QString qlineNaive = QString::fromStdString(nline);
        int correctlabel = qline.toInt();
        int naivelabel = qlineNaive.toInt();

        if(correctlabel == laap){
            naivelaap = naivelabel;
        } else if(correctlabel == lspv){
            naivelspv = naivelabel;
        } else if(correctlabel == lipv){
            naivelipv = naivelabel;
        } else if(correctlabel == rspv){
            naiverspv = naivelabel;
        } else if(correctlabel == ripv){
            naiveripv = naivelabel;
        }
    }

    if(HasSurface()){
        MITK_INFO << "Attempting to correct naive labels for default ones";

        vtkFloatArray *cellScalars = vtkFloatArray::SafeDownCast(surface->GetVtkPolyData()->GetCellData()->GetScalars());
        for (vtkIdType ix = 0; ix < surface->GetVtkPolyData()->GetNumberOfCells(); ix++) {
            double testScalar = cellScalars->GetTuple1(ix);

            if(testScalar == naivelaap){
                cellScalars->SetTuple1(ix, laap);
            } else if(testScalar == naivelspv){
                cellScalars->SetTuple1(ix, lspv);
            } else if(testScalar == naivelipv){
                cellScalars->SetTuple1(ix, lipv);
            } else if(testScalar == naiverspv){
                cellScalars->SetTuple1(ix, rspv);
            } else if(testScalar == naiveripv){
                cellScalars->SetTuple1(ix, ripv);
            }
        }
    } else{
        MITK_INFO << "Surface not loaded. Naive labelss associated to default labels, but not changed";
    }
}

void CemrgAtrialTools::ExtractLabelFromShell(QString dir, int label, QString outName){
    // make sure surface is loaded
    mitk::Surface::Pointer outputsurf = mitk::Surface::New();
    outputsurf->SetVtkPolyData(surface->GetVtkPolyData());
    CemrgCommonUtils::SetCellDataToPointData(outputsurf);

    vtkSmartPointer<vtkPolyDataConnectivityFilter> cf = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    cf->SetInputData(outputsurf->GetVtkPolyData());
    cf->ScalarConnectivityOn();
    cf->FullScalarConnectivityOn();
    cf->SetScalarRange(GetNaiveLabel(label), GetNaiveLabel(label));
    cf->Update();
    cf->SetExtractionModeToLargestRegion();

    QString path = dir + "/" + outName + ".vtk";
    outputsurf->SetVtkPolyData(cf->GetOutput());
    mitk::IOUtil::Save(outputsurf, path.toStdString());
}

void CemrgAtrialTools::FindVeinLandmarks(ImageType::Pointer im, vtkSmartPointer<vtkPolyData> pd, int nveins, QString outName){
    vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
    for (int i=0; i<pd->GetNumberOfPoints(); i++) {
        double* point = pd->GetPoint(i);
        point[0] = -point[0];
        point[1] = -point[1];
        pd->GetPoints()->SetPoint(i, point);
    }//_for
    pointLocator->SetDataSet(pd);
    pointLocator->BuildLocator();

    vtkSmartPointer<vtkIdList> pickedSeedIds = vtkSmartPointer<vtkIdList>::New();
    pickedSeedIds->Initialize();

    IteratorType itLMK(im, im->GetRequestedRegion());
    std::vector<std::vector<double>> veinsCentre;
    QString fcentresPath = directory + "/" + outName + "Centres.txt";
    std::ofstream fcentres(fcentresPath.toStdString());

    for (int j=0; j<nveins; j++) {
        int ctrVeinsVoxels = 0;
        std::vector<double> veinLandmark(3, 0.0);
        for (itLMK.GoToBegin(); !itLMK.IsAtEnd(); ++itLMK) {
            if ((int)itLMK.Get() == (j+1)) {
                ImageType::PointType point;
                im->TransformIndexToPhysicalPoint(itLMK.GetIndex(), point);
                veinLandmark[0] += point[0];
                veinLandmark[1] += point[1];
                veinLandmark[2] += point[2];
                ctrVeinsVoxels++;
            }
        }//_for
        veinLandmark[0] /= ctrVeinsVoxels;
        veinLandmark[1] /= ctrVeinsVoxels;
        veinLandmark[2] /= ctrVeinsVoxels;
        veinsCentre.push_back(veinLandmark);

        fcentres << veinLandmark[0] << "," << veinLandmark[1] << "," << veinLandmark[2] << std::endl;
    }//_nveins
    fcentres.close();

    for (int j=0; j<nveins; j++) {
        double veinLandmark[3];
        veinLandmark[0] = veinsCentre.at(j)[0];
        veinLandmark[1] = veinsCentre.at(j)[1];
        veinLandmark[2] = veinsCentre.at(j)[2];

        vtkIdType id = pointLocator->FindClosestPoint(veinLandmark);
        pickedSeedIds->InsertNextId(id);
    }//_nveins

    std::vector<int> pickedSeedLabels;
    for (int j=0; j<nveins; j++){
        pickedSeedLabels.push_back(21);
    }

    QString fidsPath = directory + "/" + outName + "Ids.txt";
    QString flabelsPath = directory + "/" + outName + "Labels.txt";
    std::ofstream fids(fidsPath.toStdString());
    std::ofstream flabels(flabelsPath.toStdString());

    for (int ix = 0; ix < nveins; ix++) {
        vtkIdType vId = ix;
        fids << pickedSeedIds->GetId(vId) << std::endl;
        flabels << pickedSeedLabels.at(ix) << std::endl;
    }
    fids.close();
    flabels.close();

}

bool CemrgAtrialTools::CheckLabelConnectivity(mitk::Surface::Pointer externalSurface, QStringList labelsToCheck, std::vector<int> &labelsVector){
    std::vector<int> v;
    bool foundBodyLabel=false;
    for (int ix = 0; ix < labelsToCheck.size(); ix++) {
        bool ok;
        int labelAtThis = labelsToCheck.at(ix).toInt(&ok);
        if(ok){
            v.push_back(labelAtThis);
            if (labelAtThis==abody){
                foundBodyLabel = true;
            }
        }
    }
    if (!foundBodyLabel){
        v.push_back(abody);
    }

    int totalWrongLabels=0;

    for (unsigned int jx = 0; jx < v.size(); jx++) {
        double currentLabel = (double) v.at(jx);
        vtkSmartPointer<vtkConnectivityFilter> cf = GetLabelConnectivity(externalSurface, currentLabel);
        int currentNumRegions = cf->GetNumberOfExtractedRegions();
        if(currentNumRegions>1){
            totalWrongLabels++;
            labelsVector.push_back(currentLabel);
            std::cout << "label incorrect: " << currentLabel << " num regions: " << currentNumRegions<< '\n';
        }
    }
    std::cout << "Number of wrong labels: " << totalWrongLabels << '\n';
    return (totalWrongLabels>0);
}

void CemrgAtrialTools::FixSingleLabelConnectivityInSurface(mitk::Surface::Pointer externalSurface, int wrongLabel, QString outpath){
    double currentLabel = (double) wrongLabel;
    bool colours = true;
    vtkSmartPointer<vtkConnectivityFilter> cf = GetLabelConnectivity(externalSurface, currentLabel, colours);
    int NumRegions = cf->GetNumberOfExtractedRegions();
    cf->Update();

    // Get max area
    vtkSmartPointer<vtkPolyDataConnectivityFilter> pdcf = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    pdcf->SetInputData(cf->GetPolyDataOutput());
    pdcf->ScalarConnectivityOn();
    pdcf->FullScalarConnectivityOn();
    int maxObjectSize = -1, maxObjectIndex=-1;
    QStringList fnames;
    for (int ix = 0; ix < NumRegions; ix++) {
        pdcf->SetScalarRange(ix, ix);
        pdcf->Update();
        pdcf->SetExtractionModeToLargestRegion();

        vtkSmartPointer<vtkMassProperties> mp = vtkSmartPointer<vtkMassProperties>::New();
        mp->SetInputConnection(pdcf->GetOutputPort());

        if(mp->GetSurfaceArea() > maxObjectSize){
            maxObjectSize = mp->GetSurfaceArea();
            maxObjectIndex = ix;
        }
    }

    vtkFloatArray *cellScalars = vtkFloatArray::SafeDownCast(externalSurface->GetVtkPolyData()->GetCellData()->GetScalars());

    vtkNew<vtkCellLocator> cellLocator;
    cellLocator->SetDataSet(externalSurface->GetVtkPolyData());
    cellLocator->BuildLocator();

    double testScalar;

    for (int jx = 0; jx < NumRegions; jx++) {
        if(jx!=maxObjectIndex){
            pdcf->SetScalarRange(jx, jx);
            pdcf->Update();
            pdcf->SetExtractionModeToLargestRegion();

            std::vector<double> labelsAroundPatch;

            vtkSmartPointer<vtkIdList> cellsAroundPatch = vtkSmartPointer<vtkIdList>::New();
            vtkSmartPointer<vtkIdList> patchCells = vtkSmartPointer<vtkIdList>::New();
            cellsAroundPatch->Initialize();
            patchCells->Initialize();

            std::cout << "Number of cells at region " << jx << ": " << pdcf->GetOutput()->GetNumberOfCells() << '\n';
            for (vtkIdType queryId=0; queryId < pdcf->GetOutput()->GetNumberOfCells(); queryId++) {
                MITK_INFO(debugMessages) << "\t get cell from piece";
                vtkSmartPointer<vtkIdList> ptsInCell = vtkSmartPointer<vtkIdList>::New();
                pdcf->GetOutput()->GetCellPoints(queryId, ptsInCell);
                vtkIdType numPtsInCell = ptsInCell->GetNumberOfIds();

                MITK_INFO(debugMessages) << "\t compute cog of piece";
                double numPts=0;
                double cog[3] = {0, 0, 0};
                for (vtkIdType ptId = 0; ptId < numPtsInCell; ++ptId) {
                    double cP[3];
                    vtkIdType thisPtId = ptsInCell->GetId(ptId);
                    pdcf->GetOutput()->GetPoint(thisPtId, cP);

                    cog[0] += cP[0];
                    cog[1] += cP[1];
                    cog[2] += cP[2];

                    numPts++;
                }
                cog[0] /= numPts;
                cog[1] /= numPts;
                cog[2] /= numPts;

                MITK_INFO(debugMessages) << "\t query cog on cellLocator";
                double closest[3];   // the coordinates of the closest point
                double dist2closest; // the squared distance to the closest point
                vtkIdType closestId; // the cell id of the cell containing the closest point
                int subId;        // this is rarely used (in triangle strips only, I believe)
                cellLocator->FindClosestPoint(cog, closest, closestId, subId, dist2closest);

                if(debugMessages){
                    std::cout << "\t Closest point: ";
                    std::cout << "\t Point: (" ", " << closest[0] << ", " << closest[1] << ", " << closest[2] << ")" << '\n';
                    std::cout << "\t Cell ID: " << closestId;
                    std::cout << "\t Sub ID: " << subId << '\n';
                }

                vtkSmartPointer<vtkIdList> cellPtIds = vtkSmartPointer<vtkIdList>::New();
                externalSurface->GetVtkPolyData()->GetCellPoints(closestId, cellPtIds);
                vtkSmartPointer<vtkIdList> neighbours = vtkSmartPointer<vtkIdList>::New();
                neighbours->Initialize();

                MITK_INFO(debugMessages) << "\t look for surrounding labels to closestId";
                for (vtkIdType i = 0; i < cellPtIds->GetNumberOfIds(); i++){
                    vtkSmartPointer<vtkIdList> ptsIdList = vtkSmartPointer<vtkIdList>::New();
                    ptsIdList->InsertNextId(cellPtIds->GetId(i)); // add one of the edge points
                    if (i + 1 == cellPtIds->GetNumberOfIds()){ // add the other edge point
                        ptsIdList->InsertNextId(cellPtIds->GetId(0));
                    } else{
                        ptsIdList->InsertNextId(cellPtIds->GetId(i + 1));
                    }

                    // get the neighbours of the cell
                    vtkSmartPointer<vtkIdList> neighborCellIds = vtkSmartPointer<vtkIdList>::New();
                    externalSurface->GetVtkPolyData()->GetCellNeighbors(closestId, ptsIdList, neighborCellIds);

                    for (vtkIdType j = 0; j < neighborCellIds->GetNumberOfIds(); j++){
                        neighbours->InsertNextId(neighborCellIds->GetId(j));
                    }
                }
                MITK_INFO(debugMessages) << "\t Got neighbours to closest cell";

                MITK_INFO(debugMessages) << ("\t Number of neighbours: " + QString::number(neighbours->GetNumberOfIds())).toStdString();
                for (vtkIdType cId = 0; cId < neighbours->GetNumberOfIds() ; cId++) {
                    testScalar = cellScalars->GetTuple1(neighbours->GetId(cId));

                    if(testScalar == (double) wrongLabel){
                        patchCells->InsertNextId(neighbours->GetId(cId));
                    } else{
                        cellsAroundPatch->InsertNextId(neighbours->GetId(cId));
                        labelsAroundPatch.push_back(testScalar);
                    }
                }
            }

            // sort labels and make them unique
            std::vector<double> uniqueLabels(labelsAroundPatch);
            std::sort(uniqueLabels.begin(), uniqueLabels.end());
            int numUniqueLabels = std::unique(uniqueLabels.begin(), uniqueLabels.end()) - uniqueLabels.begin();

            int indexOfMax=0;
            if(numUniqueLabels>1){
                std::vector<int> countLabels(numUniqueLabels, 0);
                for (int lx = 0; lx < numUniqueLabels; lx++) {
                    double thisLabel = uniqueLabels.at(lx);
                    for (unsigned int qx = 0; qx < labelsAroundPatch.size(); qx++) {
                        countLabels.at(lx) += (labelsAroundPatch.at(qx) == thisLabel) ? 1 : 0;
                    }
                }
                indexOfMax = std::distance(countLabels.begin(),std::max_element(countLabels.begin(), countLabels.end()));
            }

            // correct labels in cellScalars
            double correctedLabel = uniqueLabels.at(indexOfMax);
            for (vtkIdType cId = 0; cId < patchCells->GetNumberOfIds() ; cId++) {
                cellScalars->SetTuple1(patchCells->GetId(cId), correctedLabel);
            }
        }
    }

    externalSurface->GetVtkPolyData()->GetCellData()->SetScalars(cellScalars);
    if (!outpath.isEmpty()){
        mitk::IOUtil::Save(externalSurface, outpath.toStdString());
    }
}

//helper functions
vtkSmartPointer<vtkConnectivityFilter> CemrgAtrialTools::GetLabelConnectivity(mitk::Surface::Pointer externalSurface, double label, bool colourRegions){
    // threshold at label
    vtkSmartPointer<vtkThreshold> thres = vtkSmartPointer<vtkThreshold>::New();
    thres->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    thres->SetLowerThreshold(label);
    thres->SetUpperThreshold(label);
    thres->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, vtkDataSetAttributes::SCALARS);
    thres->SetInputData(externalSurface->GetVtkPolyData());
    thres->Update();

    // convert threshold output to polydata
    vtkSmartPointer<vtkGeometryFilter> gf = vtkSmartPointer<vtkGeometryFilter>::New();
    gf->SetInputConnection(thres->GetOutputPort());
    gf->Update();

    // connectivity filter
    vtkSmartPointer<vtkConnectivityFilter> cf = vtkSmartPointer<vtkConnectivityFilter>::New();
    cf->SetInputConnection(gf->GetOutputPort());
    cf->Update();
    if(colourRegions){
        cf->SetExtractionModeToAllRegions();
        cf->ColorRegionsOn();
    } else{
        cf->SetExtractionModeToLargestRegion();
    }
    cf->Update();

    return cf;
}

ImageType::Pointer CemrgAtrialTools::ExtractLabel(QString tag, ImageType::Pointer im, uint16_t label, uint16_t filterRadius, int maxNumObjects){
    MITK_INFO << ("Thresholding " + tag + " from clean segmentation").toStdString();
    ThresholdType::Pointer thresVeins = ThresholdImageFilter(im, label);

    ImageType::Pointer auxIm;
    if(filterRadius > 0){
        // morphological opening step
        MITK_INFO << "Morphological opening - remove speckle noise";
        ImFilterType::Pointer imopen = ImOpenFilter(thresVeins->GetOutput(), filterRadius);
        auxIm = imopen->GetOutput();

    } else{
        MITK_INFO << "Morphological opening - ignored";
        auxIm = thresVeins->GetOutput();
    }

    if (maxNumObjects>0){
        ConnectedComponentImageFilterType::Pointer conn2 = ConnectedComponentImageFilterType::New();
        conn2->SetInput(auxIm);
        conn2->Update();

        LabelShapeKeepNObjImgFilterType::Pointer keepObjs2 = LabelShapeKeepNObjImgFilterType::New();
        keepObjs2->SetInput(conn2->GetOutput());
        keepObjs2->SetBackgroundValue(0);
        keepObjs2->SetNumberOfObjects(maxNumObjects);
        keepObjs2->SetAttribute(LabelShapeKeepNObjImgFilterType::LabelObjectType::NUMBER_OF_PIXELS);
        keepObjs2->Update();

        return keepObjs2->GetOutput();
    } else{
        return auxIm;
    }
}

ImageType::Pointer CemrgAtrialTools::AddImage(ImageType::Pointer im1, ImageType::Pointer im2){
    AddFilterType::Pointer sum = AddFilterType::New();
    sum->SetInput1(im1);
    sum->SetInput2(im2);
    sum->Update();

    return sum->GetOutput();
}

ImageType::Pointer CemrgAtrialTools::ThresholdImage(ImageType::Pointer input, uint16_t lowerThres, uint16_t upperThres){
    uint16_t upper = (upperThres < lowerThres) ? lowerThres : upperThres;

    ThresholdType::Pointer thresholdOutput = ThresholdType::New();
    thresholdOutput->SetInput(input);
    thresholdOutput->SetLowerThreshold(lowerThres);
    thresholdOutput->SetUpperThreshold(upper);
    thresholdOutput->SetInsideValue(1);
    thresholdOutput->SetOutsideValue(0);
    thresholdOutput->Update();

    return thresholdOutput->GetOutput();
}

ImageType::Pointer CemrgAtrialTools::ImOpen(ImageType::Pointer input, uint16_t radius){
    StrElType structuringElement;
    structuringElement.SetRadius(static_cast<unsigned long>(radius));
    structuringElement.CreateStructuringElement();

    ImFilterType::Pointer imOpenFilter = ImFilterType::New();
    imOpenFilter->SetInput(input);
    imOpenFilter->SetKernel(structuringElement);
    imOpenFilter->Update();

    return imOpenFilter->GetOutput();
}

ThresholdType::Pointer CemrgAtrialTools::ThresholdImageFilter(ImageType::Pointer input, uint16_t thresholdVal){
    ThresholdType::Pointer thresholdOutput = ThresholdType::New();
    thresholdOutput->SetInput(input);
    thresholdOutput->SetLowerThreshold(thresholdVal);
    thresholdOutput->SetUpperThreshold(thresholdVal);
    thresholdOutput->SetInsideValue(1);
    thresholdOutput->SetOutsideValue(0);
    thresholdOutput->Update();

    return thresholdOutput;
}

ImFilterType::Pointer CemrgAtrialTools::ImOpenFilter(ImageType::Pointer input, uint16_t radius){
    StrElType structuringElement;
    structuringElement.SetRadius(static_cast<unsigned long>(radius));
    structuringElement.CreateStructuringElement();

    ImFilterType::Pointer imOpenFilter = ImFilterType::New();
    imOpenFilter->SetInput(input);
    imOpenFilter->SetKernel(structuringElement);
    imOpenFilter->Update();

    return imOpenFilter;
}

mitk::Image::Pointer CemrgAtrialTools::ImErode(ImageType::Pointer input, int vxls){
    StrElType binaryBall;
    binaryBall.SetRadius(vxls);
    binaryBall.CreateStructuringElement();
    ErosionFilterType::Pointer erosionFilter = ErosionFilterType::New();
    erosionFilter->SetInput(input);
    erosionFilter->SetKernel(binaryBall);
    erosionFilter->UpdateLargestPossibleRegion();

    mitk::Image::Pointer roiImage = mitk::ImportItkImage(erosionFilter->GetOutput())->Clone();

    return roiImage;
}

ShortImageType::Pointer CemrgAtrialTools::Uint16ToShort(ImageType::Pointer im){
    CastUint16ToShortFilterType::Pointer u16ToShort = CastUint16ToShortFilterType::New();
    u16ToShort->SetInput(im);
    return u16ToShort->GetOutput();
}

void CemrgAtrialTools::SaveImageToDisk(ImageType::Pointer im, QString dir, QString imName){
    QString outputPath = dir + "/" + imName;
    outputPath += (!imName.contains(".nii")) ? ".nii" : "";
    MITK_INFO << ("Saving image " + imName + " to: " + dir).toStdString();

    mitk::Image::Pointer outputImg = mitk::Image::New();
    mitk::CastToMitkImage(im, outputImg);
    mitk::IOUtil::Save(outputImg, outputPath.toStdString());
}

int CemrgAtrialTools::GetNaiveLabel(int l){
    int res=-1;
    if(l == laap){
        res = naivelaap;
    } else if(l == lspv){
        res = naivelspv;
    } else if(l == lipv){
        res = naivelipv;
    } else if(l == rspv){
        res = naiverspv;
    } else if(l == ripv){
        res = naiveripv;
    }
    return res;
}

int CemrgAtrialTools::IsLabel(int l){
    int res=-1;
    if(l == laap){
        res = laap;
    } else if(l == lspv){
        res = lspv;
    } else if(l == lipv){
        res = lipv;
    } else if(l == rspv){
        res = rspv;
    } else if(l == ripv){
        res = ripv;
    }
    return res;
}
