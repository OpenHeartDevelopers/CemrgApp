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
 * Scar Advanced Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * jose.solislemus@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// Qmitk
#include <mitkSurface.h>
#include <mitkIOUtil.h>
#include <mitkImageCast.h>
#include <mitkImagePixelReadAccessor.h>

// VTK
#include <vtkClipPolyData.h>
#include <vtkImplicitBoolean.h>
#include <vtkPlane.h>
#include <vtkSphere.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkIterativeClosestPointTransform.h>
#include <vtkLandmarkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkMath.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkPolyDataNormals.h>
#include <vtkMassProperties.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>

// ITK
#include <itkPoint.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>

// Qt
#include <QtDebug>
#include <QMessageBox>

// C++ Standard
#include <numeric>
#include <string>
#include <sstream>


#include "CemrgScarAdvanced.h"

CemrgScarAdvanced::CemrgScarAdvanced() {

    _SourcePolyData = vtkSmartPointer<vtkPolyData>::New();
    _fileOutName = "encirclements.csv";
    _leftrightpre = "";
    _weightedcorridor = true;
    _neighbourhood_size = 3;
    _fill_threshold = 0.5;
    _max_scalar = -1;
    _run_count = 0;
    fi1_largestSurfaceArea = -1;
    fi2_percentage = -1;
    fi2_largestSurfaceArea = -1;
    fi2_corridorSurfaceArea = -1;
    fi2_connectedAreasTotal = -1;
    fi3_preScarScoreSimple = -1;
    fi3_postScarScoreSimple = -1;
    fi3_totalPoints = -1;
    fi3_emptyPoints = -1;
    fi3_healthy = -1;
    fi3_preScar = -1;
    fi3_postScar = -1;
    fi3_overlapScar = -1;
    _debugScarAdvanced = false;
}

QString CemrgScarAdvanced::GetOutputSufix() {

    QString teststr = QString::fromStdString(_prefix);
    if (teststr.compare("pre", Qt::CaseInsensitive))
        return "Pre";
    else {
        if (teststr.compare("post", Qt::CaseInsensitive))
            return "Post";
        else
            return "";
    }
}

void CemrgScarAdvanced::ResetValues() {

    this->_pointidarray.clear();
    this->_paths.clear();
    this->_shortestPaths.clear();
    this->_pathMappers.clear();
    this->_run_count++;
}

// Helper functions
void CemrgScarAdvanced::PushBackOnPointIDArray(int pointID) {

    _pointidarray.push_back(pointID);
}

std::string CemrgScarAdvanced::ThresholdedShell(double thresho) {

    vtkSmartPointer<vtkIntArray> exploration_values = vtkSmartPointer<vtkIntArray>::New();
    vtkFloatArray* scalars = vtkFloatArray::SafeDownCast(_SourcePolyData->GetPointData()->GetScalars());
    vtkSmartPointer<vtkPolyData> temp = vtkSmartPointer<vtkPolyData>::New();
    temp->DeepCopy(_SourcePolyData);

    for (int i = 0; i < _SourcePolyData->GetNumberOfPoints(); i++) {
        // exploration_values->InsertNextTuple1(0);
        if (scalars->GetTuple1(i) >= thresho)
            exploration_values->InsertNextTuple1(1);
        else
            exploration_values->InsertNextTuple1(0);
    }

    temp->GetPointData()->SetScalars(exploration_values);

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName((this->PathAndPrefix() + "_ShellThreshold.vtk").c_str());
    writer->SetInputData(temp);
    writer->Update();

    return (this->PathAndPrefix() + "_ShellThreshold.vtk");
}

std::string CemrgScarAdvanced::ScarOverlap(vtkSmartPointer<vtkPolyData> prepd, double prethresh, vtkSmartPointer<vtkPolyData> postpd, double postthresh) {

    vtkSmartPointer<vtkIntArray> exploration_values = vtkSmartPointer<vtkIntArray>::New();
    vtkFloatArray* scalars_pre = vtkFloatArray::SafeDownCast(prepd->GetPointData()->GetScalars());
    vtkFloatArray* scalars_post = vtkFloatArray::SafeDownCast(postpd->GetPointData()->GetScalars());
    vtkSmartPointer<vtkPolyData> temp = vtkSmartPointer<vtkPolyData>::New();
    temp->DeepCopy(prepd);

    int valueassigned = 0;
    fi3_totalPoints = (double)prepd->GetNumberOfPoints();
    fi3_emptyPoints = 0.0;
    fi3_healthy = 0.0;
    fi3_preScar = 0.0;
    fi3_postScar = 0.0;
    fi3_overlapScar = 0.0;
    for (int i = 0; i < prepd->GetNumberOfPoints(); i++) {
        if (scalars_post->GetTuple1(i) == 0) { // Veins were clipped here, no value
            valueassigned = -1;
            fi3_emptyPoints++;
        } else {
            if (scalars_pre->GetTuple1(i) >= prethresh) {
                valueassigned += 1;
            }
            if (scalars_post->GetTuple1(i) >= postthresh) {
                valueassigned += 2;
            }

            if (valueassigned == 0) {
                fi3_healthy++;
            } else if (valueassigned == 1) {
                fi3_preScar++;
            } else if (valueassigned == 2) {
                fi3_postScar++;
            } else if (valueassigned == 3) {
                fi3_overlapScar++;
            }
        }

        exploration_values->InsertNextTuple1(valueassigned);
        valueassigned = 0;
    }

    temp->GetPointData()->SetScalars(exploration_values);

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName((GetOutputPath() + "ScarOverlap.vtk").c_str());
    writer->SetInputData(temp);
    writer->Update();

    return (GetOutputPath() + "ScarOverlap.vtk");
}

std::string CemrgScarAdvanced::num2str(double num, int precision) {

    std::stringstream stream;
    stream << std::fixed << std::setprecision(precision) << num;
    return stream.str();
}

void CemrgScarAdvanced::SaveStrToFile(std::string path2file, std::string filename, std::string text) {
    MITK_INFO << "[AdvancedScar] Saving to file: " + path2file + filename;
    std::ofstream outst;
    std::stringstream ss;

    ss << (path2file + filename);
    outst.open(ss.str().c_str(), std::ios_base::out);

    outst << text;
    outst.close();
}

std::string CemrgScarAdvanced::PrintThresholdResults(double mean, double stdv, double val) {

    std::string out = "RESULTS:\t" + _prefix + "-ablation\n\n" +
        "MEAN: " + num2str(mean, 2) +
        "\t STDev: " + num2str(stdv, 2) +
        "\n\n Value chosen: " + num2str(val, 1) + "\n"
        "Threshold value: " + num2str(this->_fill_threshold, 2) + "\n\n" +
        "Surface Area of Ablation: " + num2str(this->fi1_largestSurfaceArea, 2) + " mm^2 \n" +
        "Scar Score: " + num2str(this->fi1_scarScore, 2) + "%";

    SaveStrToFile(GetOutputPath(), _prefix + "_" + GetSurfaceAreaFilename(), out);

    return out;
}

std::string CemrgScarAdvanced::PrintAblationGapsResults(double mean, double stdv, double val) {

    std::string out = "RESULTS:\t" + _prefix + "-ablation\n\n" +
        "MEAN: " + num2str(mean, 2) +
        "\t STDev: " + num2str(stdv, 2) +
        "\n\n Value chosen: " + num2str(val, 1) + "\n"
        "Threshold value: " + num2str(this->_fill_threshold, 2) + "\n\n" +
        "Num connected sections: " + num2str(this->fi2_connectedAreasTotal, 0) + "\n" +
        "Percentage of scar in corridor : " + num2str(this->fi2_percentage, 2) + "%\t\n" +
        "Largest thresholded area in corridor: " + num2str(this->fi2_largestSurfaceArea, 0) + "mm^2 \n" +
        "Area of corridor: " + num2str(this->fi2_corridorSurfaceArea, 0) + " mm^2 \n";
    std::string strisweighted = (_weightedcorridor) ? "Weighted path" : "Geodesic";
    out = out + "\nShortest path calculation: " + strisweighted + "\n";

    SaveStrToFile(PathAndPrefix() + "_", GetGapsFilename(), out);

    return out;
}
std::string CemrgScarAdvanced::PrintScarOverlapResults(double valpre, double valpost) {

    std::string out = "RESULTS:\n\n Simple Scar Score in PRE-ablation: ";
    if (valpre != valpost) {
        out += num2str(this->fi3_preScarScoreSimple, 2) + "% \t(" +
            "Threshold value: " + num2str(valpre, 1) + ")\n" +
            "Simple Scar Score in POST-ablation: " +
            num2str(this->fi3_postScarScoreSimple, 2) + "% \t(" +
            "Threshold value: " + num2str(valpost, 1) + ")\n\n";
    } else {
        out += num2str(this->fi3_preScarScoreSimple, 2) + "% \n" +
            "Simple Scar Score in POST-ablation: " +
            num2str(this->fi3_postScarScoreSimple, 2) + "% \n " +
            "(Both at threshold value: " + num2str(valpre, 1) + ")" + "\n\n";
    }

    double total = fi3_totalPoints - fi3_emptyPoints;
    if (total > 0) {
        out += "\nOVERLAP RESULTS:\n\nHEALTHY %  : " +
            num2str(100 * (fi3_healthy / total)) +
            "\nPRE-SCAR % : " + num2str(100 * (fi3_preScar / total)) +
            "\nPOST-SCAR %: " + num2str(100 * (fi3_postScar / total)) +
            "\nOVERLAP %  : " + num2str(100 * (fi3_overlapScar / total));
    } else {
        MITK_WARN << ("Points: " + QString::number(total)).toStdString();
    }
    SaveStrToFile(GetOutputPath(), GetComparisonFilename(), out);
    return out;
}

std::vector<vtkSmartPointer<vtkActor> > CemrgScarAdvanced::GetPathsMappersAndActors() {

    for (unsigned int i = 0; i < this->_paths.size(); i++) {
        vtkSmartPointer<vtkActor> pathActor = vtkSmartPointer<vtkActor>::New();
        pathActor->SetMapper(this->_pathMappers[i]);
        pathActor->GetProperty()->SetColor(0, 1, 1); // Red
        pathActor->GetProperty()->SetLineWidth(5);
        this->_actors.push_back(pathActor);
    }

    return this->_actors;
}

// F&I T1
void CemrgScarAdvanced::GetSurfaceAreaFromThreshold(double thres, double maxscalar) {

    vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter =
        vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    //connectivityFilter->SetOutputPointsPrecision(outputPointsPrecision);
    connectivityFilter->ScalarConnectivityOn();
    connectivityFilter->SetScalarRange(thres, maxscalar);
    connectivityFilter->SetInputData(this->_SourcePolyData);

    connectivityFilter->Update();
    connectivityFilter->SetExtractionModeToLargestRegion();

    vtkSmartPointer<vtkMassProperties> mpwhole = vtkSmartPointer<vtkMassProperties>::New();
    mpwhole->SetInputData(this->_SourcePolyData);
    MITK_INFO << "MASS PROPERTIES (whole):";
    MITK_INFO << mpwhole->GetSurfaceArea();

    vtkSmartPointer<vtkMassProperties> mp = vtkSmartPointer<vtkMassProperties>::New();
    mp->SetInputConnection(connectivityFilter->GetOutputPort());
    MITK_INFO << "MASS PROPERTIES (threshold):";
    MITK_INFO << mp->GetSurfaceArea();

    this->fi1_largestSurfaceArea = mp->GetSurfaceArea();

}

void CemrgScarAdvanced::ScarScore(double thres) {

    vtkFloatArray* scalars = vtkFloatArray::SafeDownCast(_SourcePolyData->GetPointData()->GetScalars());
    int ctr1 = 0, ctr2 = 0;

    for (int i = 0; i < scalars->GetNumberOfTuples(); i++) {
        double value = scalars->GetValue(i);
        if (value == 0) {
            ctr1++;
            continue;
        }//_if
        if (value > thres) ctr2++;
    }
    double percentage = (ctr2 * 100.0) / (scalars->GetNumberOfTuples() - ctr1);

    fi1_scarScore = percentage;

    if (QString::fromStdString(_prefix).compare("pre", Qt::CaseInsensitive) == 0) {
        this->fi3_preScarScoreSimple = percentage;
        MITK_INFO << "PRE SCAR SCORE (" + _prefix + "): " + num2str(percentage, 2);
    } else {
        this->fi3_postScarScoreSimple = percentage;
        MITK_INFO << "POST SCAR SCORE (" + _prefix + "): " + num2str(percentage, 2);
    }
}

// F&I T2
void CemrgScarAdvanced::ExtractCorridorData(
    std::vector<vtkSmartPointer<vtkDijkstraGraphGeodesicPath> > allShortestPaths) {

    double xyz[3];
    typedef std::map<vtkIdType, int>::iterator it_type;
    std::map<vtkIdType, int> vertex_ids;
    std::vector<std::pair<int, int> > pointNeighbours;
    std::vector<int> pointIDsInCorridor;

    int count = 0;
    std::ofstream out;
    xyz[0] = 1e-10; xyz[1] = 1e-10; xyz[2] = 1e-10;

    std::stringstream ss;

    ss << this->_fileOutName;

    out.open(ss.str().c_str(), std::ios_base::app);
    out << "MainVertexSeq,VertexID,X,Y,Z,VertexDepth,MeshScalar" << std::endl;
    // the recursive order - how many levels deep around a point do you want to explore?
    // default is 3 levels deep, meaning neighbours neighbours neighbour.
    int order = _neighbourhood_size;

    // Bring all the scalars to an array
    vtkSmartPointer<vtkFloatArray> scalars = vtkFloatArray::SafeDownCast(_SourcePolyData->GetPointData()->GetScalars());

    // this will indicate what is vertices are in the exploration corridor
    vtkSmartPointer<vtkIntArray> exploration_corridor = vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray> exploration_scalars = vtkSmartPointer<vtkIntArray>::New();
    for (int i = 0; i < _SourcePolyData->GetNumberOfPoints(); i++) {
        exploration_corridor->InsertNextTuple1(0);
        exploration_scalars->InsertNextTuple1(0);
    }

    // collect all vertex ids lying in shortest path
    for (unsigned int i = 0; i < allShortestPaths.size(); i++) {
        // getting vertex id for each shortest path
        vtkIdList* vertices_in_shortest_path = allShortestPaths[i]->GetIdList();

        for (int j = 0; j < vertices_in_shortest_path->GetNumberOfIds(); j++) {
            // map avoids duplicates
            vertex_ids.insert(std::make_pair(vertices_in_shortest_path->GetId(j), -1));
            // only using keys, no associated value always -2
        }
    }

    MITK_INFO << ("[INFO] There were a total of " + QString::number(vertex_ids.size()) + " vertices in the shortest path you have selected").toStdString();

    for (it_type iterator = vertex_ids.begin(); iterator != vertex_ids.end(); ++iterator) {

        double scalar = -1;

        exploration_corridor->SetTuple1(iterator->first, 1);
        exploration_scalars->SetTuple1(iterator->first, 1);
        if (iterator->first > 0 && iterator->first < _SourcePolyData->GetNumberOfPoints()) {
            _SourcePolyData->GetPoint(iterator->first, xyz);
            scalar = scalars->GetTuple1(iterator->first);
        }
        out << count << "," << iterator->first << ","
            << xyz[0] << "," << xyz[1] << "," << xyz[2]
            << "," << 0 << "," << scalar << std::endl;
        GetNeighboursAroundPoint2(iterator->first, pointNeighbours, order);			// the key is the

        for (unsigned int j = 0; j < pointNeighbours.size(); j++) {

            pointIDsInCorridor.push_back(pointNeighbours[j].first);
            int pointNeighborID = pointNeighbours[j].first;
            int pointNeighborOrder = pointNeighbours[j].second;
            scalar = -1;
            double thresscalar = 0;

            // simple sanity check
            if (pointNeighborID > 0 && pointNeighborID < _SourcePolyData->GetNumberOfPoints()) {
                scalar = scalars->GetTuple1(pointNeighborID);
                _SourcePolyData->GetPoint(pointNeighborID, xyz);
            }

            thresscalar = scalar;
            out << count << "," << pointNeighborID << ","
                << xyz[0] << "," << xyz[1] << "," << xyz[2] << ","
                << pointNeighborOrder << "," << scalar << std::endl;

            exploration_corridor->SetTuple1(pointNeighborID, 1);
            exploration_scalars->SetTuple1(pointNeighborID, thresscalar);
        }

        pointNeighbours.clear();
        count++;
    }
    this->_corridoridarray = pointIDsInCorridor;

    vtkSmartPointer<vtkPolyData> temp = vtkSmartPointer<vtkPolyData>::New();
    temp->DeepCopy(_SourcePolyData);
    temp->GetPointData()->SetScalars(exploration_corridor);
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName((this->PathAndPrefix() + "exploration_corridor.vtk").c_str());
    writer->SetInputData(temp);
    writer->Update();
    MITK_INFO << "Saved Corridor";

    vtkSmartPointer<vtkPolyData> temp2 = vtkSmartPointer<vtkPolyData>::New();
    temp2->DeepCopy(_SourcePolyData);
    temp2->GetPointData()->SetScalars(exploration_scalars);
    vtkSmartPointer<vtkPolyDataWriter> writer2 =
        vtkSmartPointer<vtkPolyDataWriter>::New();
    writer2->SetFileName((this->PathAndPrefix() + "exploration_scalars.vtk").c_str());
    writer2->SetInputData(temp2);
    writer2->Update();
    MITK_INFO << "Saved scalars";

    vtkSmartPointer<vtkPointDataToCellData> p2c = vtkSmartPointer<vtkPointDataToCellData>::New();
    p2c->SetInputData(temp2);
    p2c->PassPointDataOn();
    p2c->Update();

    vtkSmartPointer<vtkThreshold> threshold = vtkSmartPointer<vtkThreshold>::New();
    threshold->ThresholdByUpper(_fill_threshold);
    threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, vtkDataSetAttributes::SCALARS);
    threshold->SetInputData(p2c->GetPolyDataOutput());
    threshold->Update();

    vtkSmartPointer<vtkConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkConnectivityFilter>::New();
    connectivityFilter->SetInputConnection(threshold->GetOutputPort());
    connectivityFilter->Update();
    MITK_INFO << "Normal connectivity filter: ";
    MITK_INFO << connectivityFilter->GetNumberOfExtractedRegions();
    connectivityFilter->SetExtractionModeToLargestRegion();
    connectivityFilter->Update();
    fi2_connectedAreasTotal = connectivityFilter->GetNumberOfExtractedRegions();

    vtkSmartPointer<vtkPolyDataConnectivityFilter> cf = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    //cf->SetOutputPointsPrecision(outputPointsPrecision);
    cf->SetInputData(temp2);
    cf->ScalarConnectivityOn();
    cf->FullScalarConnectivityOn();
    cf->SetScalarRange(_fill_threshold, _max_scalar);
    cf->Update();
    cf->SetExtractionModeToLargestRegion();

    vtkSmartPointer<vtkPolyDataConnectivityFilter> cf2 = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    //cf->SetOutputPointsPrecision(outputPointsPrecision);
    cf2->SetInputData(temp);
    cf2->ScalarConnectivityOn();
    cf2->FullScalarConnectivityOn();
    cf2->SetScalarRange(1, 1);
    cf2->Update();
    cf2->SetExtractionModeToLargestRegion();

    vtkSmartPointer<vtkMassProperties> mp = vtkSmartPointer<vtkMassProperties>::New();
    mp->SetInputConnection(cf->GetOutputPort());
    MITK_INFO << "SURFACE AREA IN CORRIDOR (threshold):";
    MITK_INFO << mp->GetSurfaceArea();
    fi2_largestSurfaceArea = mp->GetSurfaceArea();

    vtkSmartPointer<vtkMassProperties> mp2 = vtkSmartPointer<vtkMassProperties>::New();
    mp2->SetInputConnection(cf2->GetOutputPort());
    MITK_INFO << "SURFACE AREA IN CORRIDOR (full):";
    MITK_INFO << mp2->GetSurfaceArea();
    fi2_corridorSurfaceArea = mp2->GetSurfaceArea();

    vtkSmartPointer<vtkPolyDataWriter> writercf = vtkSmartPointer<vtkPolyDataWriter>::New();
    writercf->SetFileName((this->PathAndPrefix() + "exploration_connectivity.vtk").c_str());
    writercf->SetInputData(cf->GetOutput());
    writercf->Write();

    out.close();
}

void CemrgScarAdvanced::NeighbourhoodFillingPercentage(std::vector<int> points) {

    double fillingcounter = 0;
    double total = points.size();

    // bring al lthe scalars to an array
    vtkFloatArray* scalars = vtkFloatArray::SafeDownCast(_SourcePolyData->GetPointData()->GetScalars());
    MITK_INFO << ("[INFO] Exploring the predeteremined neighbourhood at threshold = " +
        QString::number(_fill_threshold)).toStdString();
    MITK_INFO << ("[INFO] Number of points: " + QString::number(total)).toStdString();
    for (unsigned int i = 0; i < points.size(); i++) {
        if (scalars->GetTuple1(points[i]) > _fill_threshold)
            fillingcounter++;
    }

    double percentage_in_neighbourhood = 100 * (fillingcounter / total);
    this->fi2_percentage = percentage_in_neighbourhood;
    MITK_INFO << ("[INFO] % scar in this neighbourhood = "
        + QString::number(percentage_in_neighbourhood) + ", threshold satisfy? "
        + (percentage_in_neighbourhood > _neighbourhood_size ? "Yes" : "No")).toStdString();
}

int CemrgScarAdvanced::RecursivePointNeighbours(vtkIdType pointId, int order) {

    // double s;
    if (order == 0)
        return 0;
    else {
        if (!InsertPointIntoVisitedList2(pointId, order))
            return 0;			// already visited, no need to look further down this route
        else {
            //vtkIdList* pointList = cell->GetPointIds();
            vtkIdList* pointList = vtkIdList::New();
            // get all neighbouring points of this point
            GetConnectedVertices(_SourcePolyData, pointId, pointList);

            // running through each neighbouring point
            for (vtkIdType e = 0; e < pointList->GetNumberOfIds(); e++) {
                RecursivePointNeighbours(pointList->GetId(e), order - 1);
            }
            return 1;		// keep recursing .. 0 will stop recursing .. returning just a dummy value
        }
    }
}

bool CemrgScarAdvanced::InsertPointIntoVisitedList2(vtkIdType id, int order) {

    for (unsigned int i = 0; i < _visited_point_list.size(); i++) {
        if (_visited_point_list[i].first == id)
            return false;
    }
    _visited_point_list.push_back(std::make_pair(id, order));
    return true;
}

void CemrgScarAdvanced::GetNeighboursAroundPoint2(
    int pointID, std::vector<std::pair<int, int> >& pointNeighbourAndOrder, int max_order) {

    _visited_point_list.clear();
    RecursivePointNeighbours(pointID, max_order);

    for (unsigned int i = 0; i < _visited_point_list.size(); i++) {
        //pointNeighbours.push_back(_visited_point_list[i]);
        pointNeighbourAndOrder.push_back(std::make_pair(_visited_point_list[i].first, _visited_point_list[i].second));
    }
    MITK_INFO(IsDebug()) << ("[INFO] This point has (recursive order n = " +
        QString::number(max_order) + ") = " +
        QString::number(pointNeighbourAndOrder.size()) + " neighbours").toStdString();
}

void CemrgScarAdvanced::GetConnectedVertices(
    vtkSmartPointer<vtkPolyData> mesh, int seed, vtkSmartPointer<vtkIdList> connectedVertices) {

    //Get N-order neighbours of a vertex
    //get all cells that vertex 'seed' is a part of
    vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
    mesh->GetPointCells(seed, cellIdList);

    //loop through all the cells that use the seed point
    for (vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++) {
        vtkCell* cell = mesh->GetCell(cellIdList->GetId(i));

        //if the cell doesn't have any edges, it is a line
        if (cell->GetNumberOfEdges() <= 0) {
            continue;
        }

        // Going through the edges of each cell, remember that an edge is
        // made up of two vertices
        for (vtkIdType e = 0; e < cell->GetNumberOfEdges(); e++) {
            vtkCell* edge = cell->GetEdge(e);

            vtkIdList* pointIdList = edge->GetPointIds();

            if (pointIdList->GetId(0) == seed || pointIdList->GetId(1) == seed) {
                if (pointIdList->GetId(0) == seed) {
                    connectedVertices->InsertNextId(pointIdList->GetId(1));
                } else {
                    connectedVertices->InsertNextId(pointIdList->GetId(0));
                }
            }
        }
    }
    MITK_INFO(IsDebug()) << ("[INFO] There are " + QString::number(connectedVertices->GetNumberOfIds())
        + " points connected to point " + QString::number(seed)).toStdString();
}

void CemrgScarAdvanced::getCorridorPoints(
    std::vector<vtkSmartPointer<vtkDijkstraGraphGeodesicPath> > allShortestPaths) {

    typedef std::map<vtkIdType, int>::iterator it_type;
    std::map<vtkIdType, int> vertex_ids;
    std::vector<std::pair<int, int> > pointNeighbours;
    std::vector<int> pointIDsInCorridor;
    int order = _neighbourhood_size;

    for (unsigned int i = 0; i < allShortestPaths.size(); i++) {
        vtkIdList* vertices_in_shortest_path = allShortestPaths[i]->GetIdList();

        for (int j = 0; j < vertices_in_shortest_path->GetNumberOfIds(); j++)
            vertex_ids.insert(std::make_pair(vertices_in_shortest_path->GetId(j), -1));
    }

    for (it_type iterator = vertex_ids.begin(); iterator != vertex_ids.end(); ++iterator) {
        GetNeighboursAroundPoint2(iterator->first, pointNeighbours, order);

        for (unsigned int j = 0; j < pointNeighbours.size(); j++)
            pointIDsInCorridor.push_back(pointNeighbours[j].first);

        pointNeighbours.clear();
    }
    this->_corridoridarray = pointIDsInCorridor;
}

void CemrgScarAdvanced::CorridorFromPointList(std::vector<int> points, bool circleToStart) {

    vtkSmartPointer<vtkPolyData> poly_data = this->GetSourcePolyData();
    this->_pointidarray = points;

    int lim = this->_pointidarray.size();
    for (int i = 0; i < lim; i++) {
        vtkSmartPointer<vtkDijkstraGraphGeodesicPath> dijkstra = vtkSmartPointer<vtkDijkstraGraphGeodesicPath>::New();
        dijkstra->SetInputData(poly_data);

        if (this->IsWeighted()) {
            dijkstra->UseScalarWeightsOn();
        }

        dijkstra->Update();

        if (i < lim - 1) {
            dijkstra->SetStartVertex(this->_pointidarray[i]);
            dijkstra->SetEndVertex(this->_pointidarray[i + 1]);
        } else if (circleToStart) {
            dijkstra->SetStartVertex(this->_pointidarray[i]);
            dijkstra->SetEndVertex(this->_pointidarray[0]);
        }
        dijkstra->Update();
        vtkSmartPointer<vtkPolyDataMapper> pathMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        pathMapper->SetInputConnection(dijkstra->GetOutputPort());

        this->_shortestPaths.push_back(dijkstra);
        this->_pathMappers.push_back(pathMapper);
        this->_paths.push_back(dijkstra->GetOutput());
    }

    // compute percentage encirlcement
    this->ExtractCorridorData(this->_shortestPaths);
    this->NeighbourhoodFillingPercentage(this->_corridoridarray);

    this->_pointidarray.clear();
    this->_corridoridarray.clear();
}

void CemrgScarAdvanced::SetSourceAndTarget(vtkSmartPointer<vtkPolyData> sc, vtkSmartPointer<vtkPolyData> tg) {
    _source = sc;
    _target = tg;
}

void CemrgScarAdvanced::TransformSource2Target(QString outName) {

    // Copy scalar values from target to source
    vtkSmartPointer<vtkPolyData> Output_Poly = vtkSmartPointer<vtkPolyData>::New();
    Output_Poly->DeepCopy(_source);

    vtkSmartPointer<vtkPointLocator> Target_Poly_PointLocator = vtkSmartPointer<vtkPointLocator>::New();
    Target_Poly_PointLocator->SetDataSet(_target);
    Target_Poly_PointLocator->AutomaticOn();
    Target_Poly_PointLocator->BuildLocator();

    vtkSmartPointer<vtkFloatArray> Target_Poly_Scalar = vtkFloatArray::SafeDownCast(_target->GetPointData()->GetScalars());
    vtkSmartPointer<vtkFloatArray> Output_Poly_Scalar = vtkSmartPointer<vtkFloatArray>::New();
    Output_Poly_Scalar->SetNumberOfComponents(1);

    double pStart[3];
    double _mapping_default_value = 0;
    for (vtkIdType i = 0; i < _source->GetNumberOfPoints(); ++i) {
        _source->GetPoint(i, pStart);
        vtkIdType id_on_target = Target_Poly_PointLocator->FindClosestPoint(pStart);

        float mapped_value = 0;
        if (id_on_target > 0) {
            mapped_value = Target_Poly_Scalar->GetTuple1(id_on_target);
        } else {
            mapped_value = _mapping_default_value;
        }
        Output_Poly_Scalar->InsertNextTuple1(mapped_value);
    }

    Output_Poly->GetPointData()->SetScalars(Output_Poly_Scalar);

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName((GetOutputPath()+outName.toStdString()+".vtk").c_str());
    writer->SetInputData(Output_Poly);
    writer->Write();
}
