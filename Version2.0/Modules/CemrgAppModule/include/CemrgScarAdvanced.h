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

#ifndef CemrgScarAdvanced_h
#define CemrgScarAdvanced_h

// VTK
#include <vtkAppendFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkConnectivityFilter.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkLookupTable.h>
#include <vtkCellLocator.h>
#include "vtkPolyData.h"
#include <vtkCellArray.h>
#include "vtkGenericCell.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include <vtkCellData.h>
#include <vtkPointLocator.h>
#include <vtkMath.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataReader.h>
#include <vtkPlanes.h>
#include <vtkPointData.h>
#include <vtkThreshold.h>
#include <vtkThresholdPoints.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkDijkstraGraphGeodesicPath.h>
#include <vtkTriangle.h>
#include <vtkCellDataToPointData.h>
#include <vtkPointDataToCellData.h>
#include <vtkSelectPolyData.h>
#include <vtkCamera.h>
#include <vtkImageActor.h>
#include <vtkStructuredPoints.h>
#include <vtkTubeFilter.h>
#include <vtkPolyLine.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkLineSource.h>
#include <vtkCallbackCommand.h>
#include <vtkCellPicker.h>
#include <vtkActor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <string>
#include <sstream>
#include <vtkPointPicker.h>

// Qmitk
#include <mitkImage.h>
#include <mitkPointSet.h>
#include <vtkFloatArray.h>
#include <MitkCemrgAppModuleExports.h>

class MITKCEMRGAPPMODULE_EXPORT CemrgScarAdvanced {

public:

    int _neighbourhood_size;
    int _run_count;
    double _fill_threshold;
    double _max_scalar;
    bool _weightedcorridor;
    std::string _fileOutName;
    std::string _outPath;
    std::string _prefix;
    std::string _leftrightpre;

    double fandi1_largestSurfaceArea;
    double fandi1_scarScore;
    double fandi2_percentage;
    double fandi2_largestSurfaceArea;
    double fandi2_corridorSurfaceArea;
    int fandi2_connectedAreasTotal;
    double fandi3_preSurfacePercentage;
    double fandi3_postSurfacePercentage;
    double fandi3_preScarScore;
    double fandi3_postScarScore;

    std::vector<std::pair<int, int> > _visited_point_list; // stores the neighbours around a point
    std::vector<vtkSmartPointer<vtkPolyData> > _paths; // container to store shortest paths between points
    std::vector<vtkSmartPointer<vtkPolyDataMapper> > _pathMappers;
    std::vector<vtkSmartPointer<vtkActor> > _actors;

    vtkSmartPointer<vtkCellPicker> _cell_picker;
    vtkSmartPointer<vtkPointPicker> _point_picker;
    vtkSmartPointer<vtkPolyData> _SourcePolyData;

    vtkSmartPointer<vtkPolyData> _source;
    vtkSmartPointer<vtkPolyData> _target;

    std::vector<vtkSmartPointer<vtkDijkstraGraphGeodesicPath> > _shortestPaths;
    std::vector<int> _pointidarray;
    std::vector<int> _corridoridarray;

    // Getters and setters
    bool IsWeighted();
    bool PreScoresExist();
    bool PostScoresExist();
    void SetWeightedCorridorBool(bool isw);
    void SetInputData(vtkSmartPointer<vtkPolyData> mesh);
    void SetNeighbourhoodSize(int s);
    void SetFillThreshold(double s);
    void SetMaxScalar(double s);
    void SetOutputFileName(std::string filename);
    void SetOutputPrefix(std::string prefixname);
    void SetLeftRightPrefix(std::string lrpre);
    void ClearLeftRightPrefix();
    QString GetOutputSufix();
    void SetOutputPath(std::string pathname);
    void GetConnectedVertices(
            vtkSmartPointer<vtkPolyData> mesh, int seed, vtkSmartPointer<vtkIdList> connectedVertices);
    vtkSmartPointer<vtkPolyData> GetSourcePolyData();
    int GetThresholdValue();
    void ResetValues();

    // Helper functions
    void PushBackOnPointIDArray(int pointID);
    bool isPointIDArrayEmpty();
    std::string ThresholdedShell(double thresho);
    std::string ScarOverlap(vtkSmartPointer<vtkPolyData> prepd, double prethresh, vtkSmartPointer<vtkPolyData> postpd, double posttresh);
    std::string PathAndPrefix();
    std::string GetPrefix();
    std::vector<vtkSmartPointer<vtkActor> > GetPathsMappersAndActors();
    std::string PrintAblationGapsResults(double mean, double stdv, double val);
    std::string PrintThresholdResults(double mean, double stdv, double val);
    std::string PrintScarOverlapResults(double valpre, double valpost);
    std::string num2str(double num, int precision);

    // F&I T1
    //void GetSurfaceAreaFromThreshold();
    void GetSurfaceAreaFromThreshold(double thres, double maxscalar);
    void ScarScore(double thres);

    // F&I T2
    void ExtractCorridorData(std::vector<vtkSmartPointer<vtkDijkstraGraphGeodesicPath>> allShortestPaths);
    void NeighbourhoodFillingPercentage(std::vector<int> points);
    int RecursivePointNeighbours(vtkIdType pointId, int order);
    void GetNeighboursAroundPoint2(int pointID, std::vector<std::pair<int, int>>& pointNeighbourAndOrder, int max_order);
    void getCorridorPoints(std::vector<vtkSmartPointer<vtkDijkstraGraphGeodesicPath>> allShortestPaths);
    bool InsertPointIntoVisitedList2(vtkIdType id, int order);
    void CorridorFromPointList(std::vector<int> points);

    // F&I T3
    void SetSourceAndTarget(vtkSmartPointer<vtkPolyData> sc, vtkSmartPointer<vtkPolyData> tg);
    void TransformSource2Target();

    CemrgScarAdvanced();
};
#endif // CemrgScarAdvanced_h
