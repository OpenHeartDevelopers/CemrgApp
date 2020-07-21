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
 * Simple Common Utilities
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef CemrgCommonUtils_h
#define CemrgCommonUtils_h

#include <MitkCemrgAppModuleExports.h>
#include <mitkImage.h>
#include <mitkBoundingObject.h>
#include <mitkDataNode.h>
#include <mitkDataStorage.h>

class MITKCEMRGAPPMODULE_EXPORT CemrgCommonUtils {

public:

    //Cropping Utils
    static mitk::Image::Pointer CropImage();
    static void SetImageToCut(mitk::Image::Pointer imageToCut);
    static void SetCuttingCube(mitk::BoundingObject::Pointer cuttingCube);
    static void SetImageNode(mitk::DataNode::Pointer imageNode);
    static void SetCuttingNode(mitk::DataNode::Pointer cuttingNode);
    static mitk::DataNode::Pointer GetImageNode();
    static mitk::DataNode::Pointer GetCuttingNode();

    //Sampling Utils
    static mitk::Image::Pointer Downsample(mitk::Image::Pointer image, int factor);
    static mitk::Image::Pointer IsoImageResampling(mitk::Image::Pointer image, bool reorientToRAI=true);
    static mitk::Image::Pointer IsoImageResampling(QString imPath, bool reorientToRAI=true);

    //Mesh Utils
    static mitk::Surface::Pointer LoadVTKMesh(std::string path);
    static QString m3dlibParamFileGenerator(QString dir, QString filename="param-template.par", QString thicknessCalc="0");

    //Generic
    static mitk::DataNode::Pointer AddToStorage(
            mitk::BaseData* data, std::string nodeName, mitk::DataStorage::Pointer ds);

private:

    //Cropping Utils
    static mitk::Image::Pointer imageToCut;
    static mitk::BoundingObject::Pointer cuttingCube;
    static mitk::DataNode::Pointer imageNode;
    static mitk::DataNode::Pointer cuttingNode;
};

#endif // CemrgCommonUtils_h
