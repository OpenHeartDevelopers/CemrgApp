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

// Qmitk
#include <mitkImage.h>
#include <mitkPointSet.h>
#include <vtkFloatArray.h>
#include <MitkCemrgAppModuleExports.h>

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
#include <vtkPointPicker.h>

// C++ Standard
#include <string>
#include <sstream>

class MITKCEMRGAPPMODULE_EXPORT CemrgScarAdvanced {

public:

    bool _debugScarAdvanced;

    int _neighbourhood_size;
    int _run_count;
    double _fill_threshold;
    double _max_scalar;
    bool _weightedcorridor;
    std::string _fileOutName;
    std::string _outPath;
    std::string _prefix;
    std::string _leftrightpre;

    double fi1_largestSurfaceArea, fi1_scarScore;
    double fi2_percentage, fi2_largestSurfaceArea, fi2_corridorSurfaceArea;
    int fi2_connectedAreasTotal;
    double fi3_preScarScoreSimple, fi3_postScarScoreSimple;
    double fi3_totalPoints, fi3_emptyPoints, fi3_healthy, fi3_preScar, fi3_postScar, fi3_overlapScar;
    std::string fi1_fname, fi2_fname, fi3_fname;

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
    inline bool IsDebug() { return _debugScarAdvanced; };

    inline bool IsWeighted() { return _weightedcorridor; };
    inline bool PreScoresExist() { return (fi3_preScarScoreSimple >= 0); };
    inline bool PostScoresExist() { return (fi3_postScarScoreSimple >= 0); };
    inline int GetThresholdValue() { return _fill_threshold; };
    inline vtkSmartPointer<vtkPolyData> GetSourcePolyData() { return _SourcePolyData; };

    QString GetOutputSufix();
    void GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int seed, vtkSmartPointer<vtkIdList> connectedVertices);

    inline void SetDebug(bool db) { _debugScarAdvanced = db; };
    inline void SetDebugOn() { SetDebug(true); };
    inline void SetDebugOff() { SetDebug(false); };

    inline void SetWeightedCorridorBool(bool isw) { _weightedcorridor = isw; };
    inline void SetWeightedCorridorOn() { SetWeightedCorridorBool(true); };
    inline void SetWeightedCorridorOff() { SetWeightedCorridorBool(false); };

    inline void SetNeighbourhoodSize(int s) { _neighbourhood_size = s; };
    inline void SetFillThreshold(double s) { _fill_threshold = s; };
    inline void SetMaxScalar(double s) { _max_scalar = s; };
    inline void SetInputData(vtkSmartPointer<vtkPolyData> inputmesh) { _SourcePolyData->DeepCopy(inputmesh); };
    inline void SetOutputFileName(std::string filename) { _fileOutName = filename; };
    inline void SetOutputPath(std::string pathname) { _outPath = pathname; };
    inline void SetLeftRightPrefix(std::string lrpre) { _leftrightpre = lrpre; };
    inline void SetOutputPrefix(std::string prefixname) { _prefix = _leftrightpre + prefixname; };

    inline void SetSurfaceAreaFilename(std::string fn1) { fi1_fname = fn1; };
    inline void SetGapsFilename(std::string fn2) { fi2_fname = fn2; };
    inline void SetComparisonFilename(std::string fn3) { fi3_fname = fn3; };

    inline std::string GetOutputPath() { return _outPath; };
    inline std::string GetSurfaceAreaFilename() { return fi1_fname; };
    inline std::string GetGapsFilename() { return fi2_fname; };
    inline std::string GetComparisonFilename() { return fi3_fname; };
    inline std::string GetPrefix() { return (_leftrightpre + _prefix); };
    inline std::string PathAndPrefix() { return (GetOutputPath() + GetPrefix()); };

    inline bool isPointIDArrayEmpty() { return (_pointidarray.empty()); };

    inline void ClearLeftRightPrefix() { _leftrightpre = ""; };

    void ResetValues();

    // Helper functions
    void SaveStrToFile(std::string path2file, std::string filename, std::string text);
    void PushBackOnPointIDArray(int pointID);
    std::string ThresholdedShell(double thresho);
    std::string ScarOverlap(vtkSmartPointer<vtkPolyData> prepd, double prethresh, vtkSmartPointer<vtkPolyData> postpd, double posttresh);
    std::vector<vtkSmartPointer<vtkActor> > GetPathsMappersAndActors();
    std::string PrintAblationGapsResults(double mean, double stdv, double val);
    std::string PrintThresholdResults(double mean, double stdv, double val);
    std::string PrintScarOverlapResults(double valpre, double valpost);
    std::string num2str(double num, int precision = 2);

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
    void CorridorFromPointList(std::vector<int> points, bool circleToStart = true);

    // F&I T3
    void SetSourceAndTarget(vtkSmartPointer<vtkPolyData> sc, vtkSmartPointer<vtkPolyData> tg);
    void TransformSource2Target();

    CemrgScarAdvanced();
};
#endif // CemrgScarAdvanced_h
