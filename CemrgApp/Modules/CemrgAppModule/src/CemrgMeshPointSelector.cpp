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
 * Surface Mesh Points Selector
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
 * =========================================================================*/

// Qmitk
#include <mitkProgressBar.h>
#include <mitkIOUtil.h>
#include <mitkImageCast.h>
#include <mitkITKImageImport.h>
#include <mitkMorphologicalOperations.h>


// VTK
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

// ITK
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

// Qt
#include <QtDebug>
#include <QString>
#include <QFileInfo>
#include <QProcess>
#include <QMessageBox>
#include <numeric>

#include "CemrgAtrialTools.h"
#include "CemrgMeshPointSelector.h"

CemrgMeshPointSelector::CemrgMeshPointSelector() {

    lineSeeds = vtkSmartPointer<vtkPolyData>::New();
    lineSeeds->Initialize();
    lineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());

    pointsData = std::vector<PointLabelData>(); // Unified data structure
}

bool CemrgMeshPointSelector::IsEmpty() {
    int lastLabelIndex = GetLastLabelIndex();
    return (lastLabelIndex == -1);
}

void CemrgMeshPointSelector::SetAvailableLabels(QStringList& names, std::vector<int>& labels) {
    pointsData.clear();
    for (auto ix = 0; ix < names.size(); ix++) { 
        PointLabelData pd(names.at(ix), labels.at(ix));
        pointsData.push_back(pd);
    }
}

bool CemrgMeshPointSelector::AllLabelsSet() {
    bool result = true;
    if (IsEmpty())     {
        return false;
    }

    for (const auto& point : pointsData) {
        if (!point.labelSet) {        
            result = false;
            break;
        }
    }
    return result;
}

int CemrgMeshPointSelector::AddPickedIdToLabel(int pickedSeedId, int label) {
    int whichIndex = -1;
    // look in std::vector<PointLabelData> pointsData
    for (int ix = 0; ix < static_cast<int>(pointsData.size()); ix++) {
        if (pointsData[ix].pointLabel == label) {
            if (!pointsData[ix].labelSet) {
                pointsData[ix].vtkId = pickedSeedId;
                pointsData[ix].labelSet = true;
                whichIndex = ix;
            } else {
                MITK_INFO << "Label " << label << " already set!";
            }
            break;
        }
    }
    return whichIndex;
}

void CemrgMeshPointSelector::AddPointFromSurface(mitk::Surface::Pointer surface, int pickedSeedId, int label) {
    int whichIndex = AddPickedIdToLabel(pickedSeedId, label);
    if (whichIndex == -1) {
        MITK_INFO << "Label " << label << " not found or already set!";
        return;
    }

    double *point = surface->GetVtkPolyData()->GetPoint(pickedSeedId);
    pointsData[whichIndex].coordinates = {point[0], point[1], point[2]};
    lineSeeds->GetPoints()->InsertNextPoint(point);
    lineSeeds->Modified();

    pointsData[whichIndex].index = GetLastLabelIndex() + 1;
}

void CemrgMeshPointSelector::AddPointFromUnstructuredGrid(mitk::UnstructuredGrid::Pointer grid, int pickedSeedId, int label) {
    int whichIndex = AddPickedIdToLabel(pickedSeedId, label);
    if (whichIndex == -1) {
        MITK_INFO << "Label " << label << " not found or already set!";
        return;
    }

    double *point = grid->GetVtkUnstructuredGrid()->GetPoint(pickedSeedId);
    pointsData[whichIndex].coordinates = {point[0], point[1], point[2]};
    lineSeeds->GetPoints()->InsertNextPoint(point);
    lineSeeds->Modified();

    pointsData[whichIndex].index = GetLastLabelIndex() + 1;
}

void CemrgMeshPointSelector::UpdateVtkIdFromSurface(mitk::Surface::Pointer surface, QStringList names) {
    for (int ix = 0; ix < static_cast<int>(pointsData.size()); ix++) {
        if (pointsData[ix].labelSet && names.contains(pointsData[ix].pointName)) {
            vtkIdType vtkId = surface->GetVtkPolyData()->FindPoint(pointsData[ix].coordinates.data());
            pointsData[ix].vtkId = vtkId;
        }
    }
}

void CemrgMeshPointSelector::UpdateVtkIdFromUnstructuredGrid(mitk::UnstructuredGrid::Pointer grid, QStringList names) {
    for (int ix = 0; ix < static_cast<int>(pointsData.size()); ix++) {
        if (pointsData[ix].labelSet && names.contains(pointsData[ix].pointName)) {
            double *pt = pointsData[ix].coordinates.data();
            vtkIdType vtkId = grid->GetVtkUnstructuredGrid()->FindPoint(pointsData[ix].coordinates.data());

            std::cout << "Point: (" << pt[0] << " " << pt[1] << " " << pt[2] << ") Old ID:" << pointsData[ix].vtkId <<", New ID: " << vtkId << std::endl;
            pointsData[ix].vtkId = vtkId;
        }
    }
}

int CemrgMeshPointSelector::CleanupLastPoint() {
    if (pointsData.empty()) {
        return -1;
    }

    // Remove last point from VTK objects
    vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints> points = lineSeeds->GetPoints();
    for (int i = 0; i < points->GetNumberOfPoints() - 1; i++) {
        newPoints->InsertNextPoint(points->GetPoint(i));
    }
    lineSeeds->SetPoints(newPoints);

    return UnsetLastPoint();
}

int CemrgMeshPointSelector::UnsetLastPoint() {
    int lastIndex = GetLastLabelIndex();
    if (lastIndex < 0) {
        return -1;
    }

    PointLabelData pd = pointsData[lastIndex];
    pd.labelSet = false;
    pd.vtkId = -1;
    pd.coordinates = {0.0, 0.0, 0.0};
    pd.index = -1;
    pointsData[lastIndex] = pd;

    return pd.pointLabel;
}

void CemrgMeshPointSelector::Clear() {
    lineSeeds->Reset();
    pointsData.clear();
}

std::string CemrgMeshPointSelector::ToString() {
    std::string res = "";
    for (int ix = 0; ix < static_cast<int>(pointsData.size()); ix++) {
        res += "Point " + std::to_string(ix) + ":\n";
        res += "Name: " + pointsData[ix].pointName.toStdString() + '\n';
        res += "Label: " + std::to_string(pointsData[ix].pointLabel) + '\n';
        res += "Label set: " + std::to_string(pointsData[ix].labelSet) + '\n';
        res += "VTK ID: " + std::to_string(pointsData[ix].vtkId) + '\n';
        res += "Coordinates: (" + std::to_string(pointsData[ix].coordinates[0]) + ", " + std::to_string(pointsData[ix].coordinates[1]) + ", " + std::to_string(pointsData[ix].coordinates[2]) + ")\n";
    }
    return res;
}

PointLabelData CemrgMeshPointSelector::GetPointData(QString name) {
    PointLabelData pd;
    for (const auto& point : pointsData) {
        if (point.pointName == name) {
             pd = point;
            break;
        }
    }
    return pd;
}

int CemrgMeshPointSelector::FindLabel(QString name) {
    PointLabelData pd = GetPointData(name);
    if (pd.pointName == "") {
        return -1;
    }

    return pd.pointLabel;
}

int CemrgMeshPointSelector::FindIndex(QString name) {
    PointLabelData pd = GetPointData(name);
    if (pd.pointName == "") {
        return -1;
    }

    return pd.vtkId;
}

std::vector<double> CemrgMeshPointSelector::FindPoint(QString name) {
    PointLabelData pd = GetPointData(name);
    if (pd.pointName == "") {
        return {0.0, 0.0, 0.0};
    }

    std::vector<double> point = {pd.coordinates[0], pd.coordinates[1], pd.coordinates[2]};

    return point;
}

std::string CemrgMeshPointSelector::PrintManyVtx(QStringList names) {
    std::string res = "";
    std::string aux = "";
    int count = 0;
    for (auto ix = 0; ix < names.size(); ix++) {
        int vtkid = FindIndex(names[ix]);
        if (vtkid != -1){
            aux += std::to_string(vtkid) + '\n';
            count++;
        }
    }
    if (count > 0) {
        res = std::to_string(count) + "\nextra\n" + aux;
    }
    
    return res;
}

std::string CemrgMeshPointSelector::PrintManyCoordTxt(QStringList names) {
    std::string res = "";

    for (auto ix = 0; ix < names.size(); ix++) {
        int auxIx = FindIndex(names[ix]);
        if (auxIx != -1) {
            std::vector<double> point = FindPoint(names[ix]);
            res += std::to_string(point[0]) + " " + std::to_string(point[1]) + " " + std::to_string(point[2]) + '\n';
        }
    }
    return res;
}

void CemrgMeshPointSelector::SaveToFile(QString path, QStringList names, QString type ) {
    // show error if type is not "vtx" or "coord"
    if (type != "vtx" && type != "coord") {
        MITK_WARN << "Type must be 'vtx' or 'coord'";
        return;
    }
    std::string toPrint;
    if (type == "vtx")     {
        toPrint = PrintManyVtx(names);
    }
    else if (type == "coord")     {
        toPrint = PrintManyCoordTxt(names);
    }

    std::ofstream file;
    file.open(path.toStdString());
    file << toPrint;
    file.close();
}

int CemrgMeshPointSelector::GetLastLabelIndex() {
    // get the las index of pointsData (ix, such that poointsData[ix].index is maximum)
    int lastLabelIndex = -1;
    std::vector<int> pointIndices;

    for (const auto& point : pointsData) {
        if (point.index > lastLabelIndex) {
            lastLabelIndex = point.index;
        }
    }

    return lastLabelIndex;
}
