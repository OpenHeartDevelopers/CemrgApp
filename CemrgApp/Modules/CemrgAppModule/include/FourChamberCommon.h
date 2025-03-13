#ifndef FOURCHAMBERCOMMON_H
#define FOURCHAMBERCOMMON_H

#include <QString>
#include <QFileInfo>
#include <QStringList>
#include <QJsonObject>
#include <QJsonArray>
#include <QMap>
#include <QFile>

#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

#include "CemrgCommonUtils.h"

enum ManualPoints { CYLINDERS, SLICERS, VALVE_PLAINS };


enum LabelsType {
    BACKGROUND = 0,
    LV_BP = 1,
    LV_myo = 2,
    RV_BP = 3,
    LA_BP = 4,
    RA_BP = 5,
    Ao_BP = 6,
    PArt_BP = 7,
    LPV1 = 8,  // lspv
    LPV2 = 9,  // lipv
    RPV1 = 10, // rspv
    RPV2 = 11, // ripv
    LAA = 12,
    SVC = 13,
    IVC = 14,
    LV_neck = 101,
    RV_myo = 103,
    LA_myo = 104,
    RA_myo = 105,
    Ao_wall = 106,
    PArt_wall = 107,
    MV = 201,
    TV = 202,
    AV = 203,
    PV = 204,
    plane_LPV1 = 205,
    plane_LPV2 = 206,
    plane_RPV1 = 207,
    plane_RPV2 = 208,
    plane_LAA = 209,
    plane_SVC = 210,
    plane_IVC = 211,
    LPV1_ring = 221,
    LPV2_ring = 222,
    RPV1_ring = 223,
    RPV2_ring = 224,
    LAA_ring = 225,
    SVC_ring = 226,
    IVC_ring = 227
};

class SegmentationLabels {
    private:
        std::unordered_map<LabelsType, unsigned int> labelMap;

        void UpdateLabels(const SegmentationLabels& other) {
            for (const auto& pair : other.labelMap) {
                labelMap[pair.first] = pair.second; // Update the labelMap
            }
        }

    public:
        SegmentationLabels() {
            labelMap = {
                {BACKGROUND, 0},
                {LV_BP,  1},
                {LV_myo,  2},
                {RV_BP,  3},
                {LA_BP,  4},
                {RA_BP,  5},
                {Ao_BP,  6},
                {PArt_BP,  7},
                {LPV1,  8},
                {LPV2,  9},
                {RPV1,  10},
                {RPV2,  11},
                {LAA,  12},
                {SVC,  13},
                {IVC,  14},
                {LV_neck,  101},
                {RV_myo,  103},
                {LA_myo,  104},
                {RA_myo,  105},
                {Ao_wall,  106},
                {PArt_wall,  107},
                {MV,  201},
                {TV,  202},
                {AV,  203},
                {PV,  204},
                {plane_LPV1,  205},
                {plane_LPV2,  206},
                {plane_RPV1,  207},
                {plane_RPV2,  208},
                {plane_LAA,  209},
                {plane_SVC,  210},
                {plane_IVC,  211},
                {LPV1_ring,  221},
                {LPV2_ring,  222},
                {RPV1_ring,  223},
                {RPV2_ring,  224},
                {LAA_ring,  225},
                {SVC_ring,  226},
                {IVC_ring,  227}
            };
        }

        SegmentationLabels(const SegmentationLabels& other) {
            Update(other);
        }

        unsigned int Get(LabelsType labelType) const {
            auto it = labelMap.find(labelType);
            if (it != labelMap.end()) {
                return it->second;
            }
            return 0; // Default to BACKGROUND
        }

        unsigned int GetLabelFromString(const std::string &labelNameString) const {
            for (const auto &pair : labelMap) {
                if (LabelName(pair.first) == labelNameString) {
                    return pair.second;
                }
            }
            return 0; // Default to BACKGROUND
        }

        void Set(LabelsType labelType, unsigned int value) {
            labelMap[labelType] = value;
        }

        void SetLabelFromString(const std::string &labelNameString, unsigned int newTag) {
            for (auto &pair : labelMap) {
                if (LabelName(pair.first) == labelNameString) {
                    pair.second = newTag;
                    return; // Exit the loop once the label is found and updated
                }
            }
        }

        void Update(const SegmentationLabels& other) {
            std::unordered_map<LabelsType, unsigned int> otherMap;
            other.GetMap(otherMap); // Retrieve the labelMap from the other instance
            labelMap = otherMap;    // Synchronize the maps
        }

        void GetMap(std::unordered_map<LabelsType, unsigned int>& map) const {
            map = labelMap; // Fill the provided map with the unsigned internal labelMap
        }

        std::string LabelName(LabelsType labelType) const {
            switch (labelType) {
                case BACKGROUND: return "BACKGROUND"; break; 
                case LV_BP: return "LV_BP"; break; 
                case LV_myo: return "LV_myo"; break; 
                case RV_BP: return "RV_BP"; break; 
                case LA_BP: return "LA_BP"; break; 
                case RA_BP: return "RA_BP"; break; 
                case Ao_BP: return "Ao_BP"; break; 
                case PArt_BP: return "PArt_BP"; break; 
                case LPV1: return "LPV1"; break; 
                case LPV2: return "LPV2"; break; 
                case RPV1: return "RPV1"; break; 
                case RPV2: return "RPV2"; break; 
                case LAA: return "LAA"; break; 
                case SVC: return "SVC"; break; 
                case IVC: return "IVC"; break; 
                case LV_neck: return "LV_neck"; break; 
                case RV_myo: return "RV_myo"; break; 
                case LA_myo: return "LA_myo"; break; 
                case RA_myo: return "RA_myo"; break; 
                case Ao_wall: return "Ao_wall"; break; 
                case PArt_wall: return "PArt_wall"; break; 
                case MV: return "MV"; break; 
                case TV: return "TV"; break; 
                case AV: return "AV"; break; 
                case PV: return "PV"; break; 
                case plane_LPV1: return "plane_LPV1"; break; 
                case plane_LPV2: return "plane_LPV2"; break; 
                case plane_RPV1: return "plane_RPV1"; break; 
                case plane_RPV2: return "plane_RPV2"; break; 
                case plane_LAA: return "plane_LAA"; break; 
                case plane_SVC: return "plane_SVC"; break; 
                case plane_IVC: return "plane_IVC"; break; 
                case LPV1_ring: return "LPV1_ring"; break; 
                case LPV2_ring: return "LPV2_ring"; break; 
                case RPV1_ring: return "RPV1_ring"; break; 
                case RPV2_ring: return "RPV2_ring"; break; 
                case LAA_ring: return "LAA_ring"; break; 
                case SVC_ring: return "SVC_ring"; break; 
                case IVC_ring: return "IVC_ring"; break; 
                default: return "BACKGROUND";
            }
        }


        bool LabelExists(unsigned int label, LabelsType& ltt) const {
            ltt = BACKGROUND;
            for (const auto &pair : labelMap) {
                if (pair.second == label) {
                    ltt = pair.first;
                    return true;
                }
            }
            return false;
        }

        unsigned int GetLargestLabel() const {
            unsigned int largestLabel = 0;
            for (const auto &pair : labelMap) {
                if (pair.second > largestLabel) {
                    largestLabel = pair.second;
                }
            }
            return largestLabel;
        }

        unsigned int GenerateNewLabel() const {
            return GetLargestLabel() + 1;
        }

        void LoadJsonObject(QJsonObject json) {
            QStringList keys = json.keys();
            MITK_INFO << "Loading json object.";
            for (const auto &key : keys) {
                unsigned int value = json[key].toInt();
                std::cout << "Key: " << key.toStdString() << " Value: " << value << '\n';
                SetLabelFromString(key.toStdString(), json[key].toInt());
            }
        }

        // functions to convert to json
        void LabelInfoLists(QStringList& labelsKeys, QStringList& labelsValues, QStringList& labelsTypes) {
            for (const auto &pair : labelMap) {
                labelsKeys << QString::fromStdString(LabelName(pair.first));
                labelsValues << QString::number(pair.second);
                labelsTypes << "int";
            }
        }

        QStringList LabelNames() const {
            QStringList labelNames;
            for (const auto &pair : labelMap) {
                labelNames << QString::fromStdString(LabelName(pair.first));
            }
            return labelNames;
        }

        QStringList LabelTags() const {
            QStringList labelTags;
            for (const auto &pair : labelMap) {
                labelTags << QString::number(pair.second);
            }
            return labelTags;
        }

        SegmentationLabels& operator=(const SegmentationLabels& other) {
            if (this != &other) {
                Update(other);
            }
            return *this;
        }

        void toString() { 
            for (const auto &pair : labelMap) {
                std::cout << "Label: " << LabelName(pair.first) << " Tag: " << pair.second << '\n';
            }
        }

        std::vector<int> GetMeshingLabels() { 
            std::vector<int> meshingLabels = {
                LV_myo, RV_myo, LA_myo, RA_myo, Ao_wall, PArt_wall, MV, TV, AV, PV,
                plane_LPV1, plane_LPV2, plane_RPV1, plane_RPV2, plane_LAA, plane_SVC, plane_IVC,
                LAA_ring, SVC_ring, IVC_ring, LPV1_ring, LPV2_ring, RPV1_ring, RPV2_ring
            };

            return meshingLabels;
        }

        QString ExtractMeshingLabels() {
            std::vector<int> meshingLabels = GetMeshingLabels();
            QString meshingLabelsString = "";
            for (const auto &label : meshingLabels) {
                meshingLabelsString += QString::number(Get((LabelsType)label)) + ",";
            }
            meshingLabelsString.chop(1); // Remove the last comma

            return meshingLabelsString;
        }
};

enum MeshLabelsType {
    LV_mesh = 1,
    RV_mesh = 2,
    LA_mesh = 3,
    RA_mesh = 4,
    Ao_mesh = 5,
    PArt_mesh = 6,
    MV_mesh = 7,
    TV_mesh = 8,
    AV_mesh = 9,
    PV_mesh = 10,
    LSPV_mesh = 11,
    LIPV_mesh = 12,
    RSPV_mesh = 13,
    RIPV_mesh = 14,
    LAA_mesh = 15,
    SVC_mesh = 16,
    IVC_mesh = 17,
    LAA_ring_mesh = 18,
    SVC_ring_mesh = 19,
    IVC_ring_mesh = 20,
    LSPV_ring_mesh = 21,
    LIPV_ring_mesh = 22,
    RSPV_ring_mesh = 23,
    RIPV_ring_mesh = 24
};

class MeshingLabels { 
    private:
        std::unordered_map<MeshLabelsType, unsigned int> labelMap;

        void UpdateLabels(const MeshingLabels& other) {
            for (const auto& pair : other.labelMap) {
                labelMap[pair.first] = pair.second; // Update the labelMap
            }
        }
    
    public: 
        MeshingLabels() {
            labelMap = {
                {LV_mesh, 1},
                {RV_mesh, 2},
                {LA_mesh, 3},
                {RA_mesh, 4},
                {Ao_mesh, 5},
                {PArt_mesh, 6},
                {MV_mesh, 7},
                {TV_mesh, 8},
                {AV_mesh, 9},
                {PV_mesh, 10},
                {LSPV_mesh, 11},
                {LIPV_mesh, 12},
                {RSPV_mesh, 13},
                {RIPV_mesh, 14},
                {LAA_mesh, 15},
                {SVC_mesh, 16},
                {IVC_mesh, 17},
                {LAA_ring_mesh, 18},
                {SVC_ring_mesh, 19},
                {IVC_ring_mesh, 20},
                {LSPV_ring_mesh, 21},
                {LIPV_ring_mesh, 22},
                {RSPV_ring_mesh, 23},
                {RIPV_ring_mesh, 24}
            };
        }

        MeshingLabels(const MeshingLabels& other) {
            Update(other);
        }

        unsigned int Get(MeshLabelsType labelType) const {
            auto it = labelMap.find(labelType);
            if (it != labelMap.end()) {
                return it->second;
            }
            return 0; // Default to 0
        }

        void Set(MeshLabelsType labelType, unsigned int value) {
            labelMap[labelType] = value;
        }

        void Update(const MeshingLabels& other) {
            std::unordered_map<MeshLabelsType, unsigned int> otherMap;
            other.GetMap(otherMap); // Retrieve the labelMap from the other instance
            labelMap = otherMap;    // Synchronize the maps
        }

        void GetMap(std::unordered_map<MeshLabelsType, unsigned int>& map) const {
            map = labelMap; // Fill the provided map with the unsigned internal labelMap
        }

        unsigned int GetLargestLabel() const {
            unsigned int largestLabel = 0;
            for (const auto &pair : labelMap) {
                if (pair.second > largestLabel) {
                    largestLabel = pair.second;
                }
            }
            return largestLabel;
        }

        unsigned int GenerateNewLabel() const {
            return GetLargestLabel() + 1;
        }

        std::vector<int> GetLabelsVector() {
            std::vector<int> labelsVector = {LV_mesh, RV_mesh, LA_mesh, RA_mesh, Ao_mesh, PArt_mesh, MV_mesh, TV_mesh, AV_mesh, 
                                            PV_mesh, LSPV_mesh, LIPV_mesh, RSPV_mesh, RIPV_mesh, LAA_mesh, SVC_mesh, IVC_mesh, 
                                            LAA_ring_mesh, SVC_ring_mesh, IVC_ring_mesh, LSPV_ring_mesh, LIPV_ring_mesh, RSPV_ring_mesh, RIPV_ring_mesh};
            
            return labelsVector;
        }

        std::vector<int> SetLabelCorrespondence(SegmentationLabels& segLabelsObject) {
            std::vector<int> segLabels = segLabelsObject.GetMeshingLabels();
            std::vector<int> meshLabels = GetLabelsVector();

            int numLabels = segLabels.size();

            std::vector<int> correspondenceVector(numLabels * 2);
            std::vector<int> correspondenceVectorAux;

            for (int ix = 0; ix < numLabels; ix++) {
                int sourceId = segLabels.at(ix);
                int targetId = meshLabels.at(ix);
                
                unsigned int sourceLabel = segLabelsObject.Get((LabelsType)sourceId);
                int targetLabel = Get((MeshLabelsType)targetId);

                auto it = std::find(segLabels.begin(), segLabels.end(), targetLabel);
                if (it != segLabels.end()) { // target exists in segLabels
                    std::cout << "Target label " << targetLabel << " already exists in segmentation labels.\n";
                    int auxLabel = GenerateNewLabel();
                    
                    std::cout << "Source: " << sourceLabel << " aux Target: " << auxLabel << ". Then, Target: " << targetLabel << '\n';
                    correspondenceVector.at(2*ix) = sourceLabel;
                    correspondenceVector.at(2*ix + 1) = auxLabel;

                    correspondenceVectorAux.push_back(auxLabel);
                    correspondenceVectorAux.push_back(targetLabel);
                } else {
                    std::cout << "Source: " << sourceLabel << " Target: " << targetLabel << '\n';

                    correspondenceVector.at(2*ix) = sourceLabel;
                    correspondenceVector.at(2*ix + 1) = targetLabel;

                }
            }

            // append auxCorr to correspondence
            correspondenceVector.insert(correspondenceVector.end(), correspondenceVectorAux.begin(), correspondenceVectorAux.end());
            
            return correspondenceVector;
        }

        void RelabelMesh(QString elemFile, SegmentationLabels& segLabelsObject) {
            std::vector<int> correspondence = SetLabelCorrespondence(segLabelsObject);
            unsigned int numCorrespondence = correspondence.size()/2;
            QFileInfo fi(elemFile);
            QString oldElemFile = fi.absolutePath() + "/" + fi.baseName() + "_old_labels.elem";
            QFile::copy(elemFile, oldElemFile);

            int nElem;
            std::ifstream elemFileRead;
            elemFileRead.open(oldElemFile.toStdString());
            if (!elemFileRead.is_open()) {
                MITK_INFO << ("Error: Failed to open input file:" + oldElemFile).toStdString();
                return;
            }

            elemFileRead >> nElem;
            std::vector<int> values(4 * nElem); 
            std::vector<int> labels(nElem);
            MITK_INFO << "Reading elem file.";
            for (int ix = 0; ix < nElem; ix++) {
                std::string type;
                elemFileRead >> type;
                elemFileRead >> values.at(4*ix + 0);
                elemFileRead >> values.at(4*ix + 1);
                elemFileRead >> values.at(4*ix + 2);
                if (type.compare("Tt") == 0) {
                    elemFileRead >> values.at(4*ix + 3);
                } else {
                    values.at(4*ix + 3) = -1;
                }
                
                elemFileRead >> labels.at(ix); 
            }

            elemFileRead.close();
            
            std::vector<int> newLabels(nElem);
            for (int ix = 0; ix < nElem; ix++) {
                newLabels.at(ix) = labels.at(ix);
            }

            for (unsigned int ix = 0; ix < numCorrespondence; ix++) {
                int source = correspondence.at(2*ix);
                int target = correspondence.at((2*ix) + 1);
                MITK_INFO << ("Updating label: " + QString::number(source) + " to " + QString::number(target)).toStdString();

                for (int jx = 0; jx < nElem; jx++) {
                    if (newLabels.at(jx) == source) {
                        newLabels.at(jx) = target;
                    }
                }
            }

            MITK_INFO << "Writing new elem file.";
            std::ofstream elemFileWrite(elemFile.toStdString());
            elemFileWrite << nElem << '\n';

            for (int ix = 0; ix < nElem; ix++) {
                bool isTetra = values.at(4 * ix + 3) == -1;

                std::string type = (isTetra) ? "Tr" : "Tt";
                elemFileWrite << type << " " << values.at(4*ix + 0) << " " << values.at(4*ix + 1) << " " << values.at(4*ix + 2);
                if (isTetra) {
                    elemFileWrite << " ";
                } else {
                    elemFileWrite << " " << values.at(4*ix + 3) << " ";
                }

                elemFileWrite << newLabels.at(ix) << '\n';
            }

            elemFileWrite.close();
        }

        QString GetMeshingLabelString(MeshLabelsType labelType) {
            QString res;
            switch (labelType) {
            case LV_mesh : res = "LV" ; break; 
            case RV_mesh : res = "RV" ; break;
            case LA_mesh : res = "LA" ; break;
            case RA_mesh : res = "RA" ; break;
            case Ao_mesh : res = "Ao" ; break;  
            case PArt_mesh : res = "PArt" ; break;
            case MV_mesh : res = "MV" ; break;
            case TV_mesh : res = "TV" ; break;
            case AV_mesh : res = "AV" ; break;
            case PV_mesh : res = "PV" ; break;
            case LSPV_mesh : res = "LSPV" ; break;
            case LIPV_mesh : res = "LIPV" ; break;
            case RSPV_mesh : res = "RSPV" ; break;
            case RIPV_mesh : res = "RIPV" ; break;
            case LAA_mesh : res = "LAA" ; break;
            case SVC_mesh : res = "SVC" ; break;
            case IVC_mesh : res = "IVC" ; break;
            case LAA_ring_mesh : res = "LAA_ring" ; break;
            case SVC_ring_mesh : res = "SVC_ring" ; break;
            case IVC_ring_mesh : res = "IVC_ring" ; break;
            case LSPV_ring_mesh : res = "LSPV_ring" ; break;
            case LIPV_ring_mesh : res = "LIPV_ring" ; break;
            case RSPV_ring_mesh : res = "RSPV_ring" ; break;
            case RIPV_ring_mesh : res = "RIPV_ring" ; break;
            default: res = "UNKNOWN" ; break;
            }
            return res;
        }

        void SaveToJson(QString pathToSave) { 
            // Save to json
            QStringList keys, values, types;
            for (const auto &pair : labelMap) {
                keys << GetMeshingLabelString((MeshLabelsType)pair.first);
                values << QString::number(pair.second);
                types << "int";
            }
            QJsonObject json = CemrgCommonUtils::CreateJSONObject(keys, values, types);
            QFileInfo fi(pathToSave);
            bool success = CemrgCommonUtils::WriteJSONFile(json, fi.absolutePath(), fi.baseName() + ".json");
            MITK_INFO(success) << "Meshing labels saved to " << pathToSave.toStdString();
        }
};

struct FourChamberSubfolders {
    QString SEG, MESH, UVC, UVC_LA, UVC_RA, AFIB, PRESIM, SIMS, PAR;
    FourChamberSubfolders()
    :SEG("segmentations"),
     MESH("meshing"),
     UVC("surfaces_uvc"),
     UVC_LA("surfaces_uvc_LA"),
     UVC_RA("surfaces_uvc_RA"),
     AFIB("atrial_fibres"),
     PRESIM("pre_simulation"),
     SIMS("sims_folder"), 
     PAR("parfiles")
    {}

    QStringList Subdirectories(){
         return (QStringList() << SEG << MESH << UVC << UVC_LA << UVC_RA << AFIB << PRESIM << SIMS << PAR);
    }
};

struct ManualPointsStruct {
    QStringList cylinders, slicers, valve_plains;
    ManualPointsStruct(){
        cylinders << "SVC_1" << "SVC_2" << "SVC_3" << "IVC_1" << "IVC_2" << "IVC_3"  
                 << "Ao_1"  << "Ao_2"  << "Ao_3" << "PArt_1" << "PArt_2" << "PArt_3";
        slicers << "SVC_slicer_1" << "SVC_slicer_2" << "SVC_slicer_3" << "SVC_tip"
                << "IVC_slicer_1" << "IVC_slicer_2" << "IVC_slicer_3" << "IVC_tip"
                << "Ao_tip" << "PArt_tip";
        valve_plains << "Ao_WT_tip" << "PArt_WT_tip";
    }

    QStringList CYLINDERS(){ return cylinders; };
    QStringList SLICERS(){ return slicers; };
    QStringList VALVE_PLAINS() { return valve_plains; };

    QStringList GetPointLabelOptions(ManualPoints mpt){
        QStringList res = QStringList();

        switch (mpt) {
            case ManualPoints::CYLINDERS :
                res = CYLINDERS();
                break;

            case ManualPoints::SLICERS : 
                res = SLICERS();
                break;

            case ManualPoints::VALVE_PLAINS :
                res = VALVE_PLAINS();
                break;
            default:
                break;
            return res;
        }
        return res;
    }

    QString title(ManualPoints mpt) {
        QString res = QString();
        switch (mpt) {
            case ManualPoints::CYLINDERS :
                res = "Cylinders";
                break;

            case ManualPoints::SLICERS : 
                res = "Slicers";
                break;

            case ManualPoints::VALVE_PLAINS :
                res = "Valve Plains";
                break;
            default:
                break;
            return res;
        }
        return res;
    }
};

/// @brief Meshtool3D static libraries parameters structures
struct M3DParameters {
    QString seg_dir, seg_name, out_dir, out_name, working_mesh;

    bool mesh_from_segmentation, boundary_relabeling, verbose, eval_thickness;

    double facet_angle, facet_size, facet_distance;
    double cell_rad_edge_ratio, cell_size;
    double rescaleFactor, abs_tol, rel_tol;
    int itr_max, dimKrilovSp;
    
    bool  out_medit, out_carp, out_carp_binary, out_vtk, out_vtk_binary, out_potential;

    M3DParameters()
    :seg_dir(""),
     seg_name(""),
     out_dir(""),
     out_name(""), 
     working_mesh(""),
     mesh_from_segmentation(true),
     boundary_relabeling(false),
     verbose(true),
     eval_thickness(false),
     facet_angle(30.0),
     facet_size(0.8),
     facet_distance(4),
     cell_rad_edge_ratio(2.0),
     cell_size(0.8),
     rescaleFactor(1000),
     abs_tol(1e-6),
     rel_tol(1e-6),
     itr_max(700),
     dimKrilovSp(500),
     out_medit(false),
     out_carp(true),
     out_carp_binary(false),
     out_vtk(false),
     out_vtk_binary(false),
     out_potential(false)
    {}

    bool CheckFormats(){
        return !(out_medit || out_carp || out_carp_binary ||out_vtk || out_vtk_binary);
    }

    void SetBinaryFormats(){
        out_carp_binary = out_carp;
        out_vtk_binary = out_vtk;
    }

    void KeysAndValues(QStringList& keys, QStringList& values, QStringList& types){
        keys << "seg_dir" << "seg_name" << "mesh_from_segmentation" << "boundary_relabeling" << 
                "facet_angle" << "facet_size" << "facet_distance" << "cell_rad_edge_ratio" << "cell_size" << "rescaleFactor" << 
                "abs_tol" << "rel_tol" << "itr_max" << "dimKrilovSp" << "verbose" << "eval_thickness" << 
                "out_dir" << "out_name" << "out_medit" << "out_carp" << "out_carp_binary" << "out_vtk" << "out_vtk_binary" << "out_potential" << "working_mesh";
        values << seg_dir << seg_name << QString::number(mesh_from_segmentation) << QString::number(boundary_relabeling)
               << QString::number(facet_angle) << QString::number(facet_size) << QString::number(facet_distance) << QString::number(cell_rad_edge_ratio) 
               << QString::number(cell_size) << QString::number(rescaleFactor) << QString::number(abs_tol)
               << QString::number(rel_tol) << QString::number(itr_max) << QString::number(dimKrilovSp) << QString::number(verbose) << QString::number(eval_thickness) 
               << out_dir << out_name
               << QString::number(out_medit) << QString::number(out_carp) << QString::number(out_carp_binary) << QString::number(out_vtk) << QString::number(out_vtk_binary) << QString::number(out_potential) <<
               working_mesh;
        types << "string" << "string" << "int" << "int" << 
                "float" << "float" << "float" << "float" << "float" << "int" << 
                "float" << "float" << "int" << "int" << "int" << "int" << 
                "string" << "string" << "int" << "int" << "int" << "int" << "int" << "int" << "string";
    }

};

struct ParfilesNamesStruct {
    QString dir, dirApexSeptumTemplates, dirEtags;
    QString tagsVentFibres, tagsPresim, tagsLvrv, tagsAtrialFibres, settingsBachmannBundleFec, SettingsAtriaMapSettings;
    QStringList apexSeptumTemplatesList, etagsList;

    ParfilesNamesStruct()
    : dir(""),
      dirApexSeptumTemplates("apex_septum_templates"),
      dirEtags("etags"),
      tagsVentFibres("tags_vent_fibres.json"),
      tagsPresim("tags_presim.json"),
      tagsLvrv("tags_lvrv.json"),
      tagsAtrialFibres("tags_atrial_fibres.json"),
      settingsBachmannBundleFec("bachmann_bundles_fec_settings.json"),
      SettingsAtriaMapSettings("atria_map_settings.json")
    {
        apexSeptumTemplatesList = QStringList() << "la.lvapex.vtx" << "la.rvsept_pt.vtx" << "raa_apex.txt"  "ra.lvapex.vtx" << "ra.rvsept_pt.vtx";
        etagsList = QStringList() << "etags_la.sh" << "etags_ra.sh" << "etags.sh";
    }

    inline void SetDirectory(QString directory) { dir = directory; };
    
    void ResetApexSeptumTemplates() {
        QString apexSeptumTemplates = dir + "/" + dirApexSeptumTemplates;
        foreach (QString file, apexSeptumTemplatesList) {
            QString filePath = apexSeptumTemplates + "/" + file;
            if (QFile::exists(filePath)) {
                QFile::remove(filePath);
            }

            std::ofstream fo(filePath.toStdString());
            if (file.contains("raa_apex")) {
                fo << "replace this text with the 3 coordinates of the RAA_apex";
            }
            else {
                fo << "1\nextra\nsave_vtx_here";
            }
        }    
    }

    QString VentFibres() { return dir + "/" + tagsVentFibres; };
    QString ApexSeptumTemplates() { return dir + "/" + dirApexSeptumTemplates; };
    QString Etags() { return dir + "/" + dirEtags; };
    QString Presim() { return dir + "/" + tagsPresim; };
    QString LVRV() { return dir + "/" + tagsLvrv; };
    QString AtrialFibres() { return dir + "/" + tagsAtrialFibres; };
    QString BachmannBundleFecSettings() { return dir + "/" + settingsBachmannBundleFec; };
    QString AtriaMapSettings() { return dir + "/" + SettingsAtriaMapSettings; };
};

struct VentricularFibresParams
{
    double _alpha_endo, _alpha_epi, _beta_endo, _beta_epi;
    QString _bad_elem, _type, _apex_to_base, _epi, _lv, _rv;
    QString _directory, _meshname, _output;
    int _nonmyo;
    bool _nomyo_set;

    VentricularFibresParams()
    :_alpha_endo(60),
     _alpha_epi(-60),
     _beta_endo(-65), 
     _beta_epi(25), 
     _bad_elem(""), 
     _type("biv"), 
     _apex_to_base("uvc/BiV.sol_apba_lap.dat"), 
     _epi("uvc/BiV.sol_endoepi_lap.dat"),
     _lv("uvc/BiV.sol_lvendo_lap.dat"),
     _rv("uvc/BiV.sol_rvendo_lap.dat"), 
     _meshname("BiV"),
     _nonmyo(-1), 
     _nomyo_set(false)
     {}

    inline void SetDirectory(QString dir) { _directory = dir; };
    inline void SetMeshname(QString name) { _meshname = name; };
    inline void SetOutput(QString out) { _output = out; };
    inline void SetDefaultOutput() { 
        _output = "fibres_bayer_" + a_endo() + "_" + a_epi() + ".lon"; 
    }

    inline void SetBadElem(QString bad) { _bad_elem = bad; };
    inline void SetApexToBase(QString apex) { _apex_to_base = apex; };

    inline void SetEpi(QString epi) { _epi = epi; };
    inline void SetLV(QString lv) { _lv = lv; };
    inline void SetRV(QString rv) { _rv = rv; };

    inline void SetNonmyo(int nonmyo) { _nonmyo = nonmyo; _nomyo_set = true; };
    inline void SetType(QString type) { _type = type; };

    inline void SetAlphaEndo(double alpha) { _alpha_endo = alpha; };
    inline void SetAlphaEpi(double alpha) { _alpha_epi = alpha; };
    inline void SetBetaEndo(double beta) { _beta_endo = beta; };
    inline void SetBetaEpi(double beta) { _beta_epi = beta; };

    inline QString a_endo() { return QString::number(_alpha_endo); }; 
    inline QString a_epi() { return QString::number(_alpha_epi); }; 
    inline QString b_endo() { return QString::number(_beta_endo); }; 
    inline QString b_epi() { return QString::number(_beta_epi); };
    inline QString nonmyo() { return QString::number(_nonmyo); };
    inline QString type() { return _type; };

    inline QString output() { return path(_output); };
    inline QString bad_elem() { return path(_bad_elem); };
    inline QString apex_to_base() { return path(_apex_to_base); };
    inline QString meshname() { return path(_meshname); };
    inline QString epi() { return path(_epi); };
    inline QString lv() { return path(_lv); };
    inline QString rv() { return path(_rv); };

    inline QString directory() { return _directory; };
    inline QString path(QString fname = "") { return fname.isEmpty() ? "" : _directory + "/" + fname; };

    void KeysAndValues(QStringList& keys, QStringList& values, QStringList& types) {
        keys << "alpha_endo" << "alpha_epi" << "beta_endo" << "beta_epi" << "bad_elem" << "type" << "apex_to_base" 
            << "epi" << "lv" << "rv" << "nonmyo" << "directory" << "meshname" << "output";
        values << a_endo() << a_epi() << b_endo() << b_epi() << bad_elem() << type() << apex_to_base() 
            << epi() << lv() << rv() << nonmyo() << directory() << meshname() << output();
        types << "float" << "float" << "float" << "float" << "string" << "string" << "string" 
            << "string" << "string" << "string" << "int" << "string" << "string" << "string";
    }
};

enum PointsNamesType {
    SVC1, SVC2, SVC3, IVC1, IVC2, IVC3, Ao1, Ao2, Ao3, PArt1, PArt2, PArt3,
    SVC_SLICER1, SVC_SLICER2, SVC_SLICER3, SVC_TIP, IVC_SLICER1, IVC_SLICER2, IVC_SLICER3, IVC_TIP, Ao_TIP, PArt_TIP,
    Ao_WT_TIP, PArt_WT_TIP
};

class BasePointsType {

public:
    BasePointsType(int numEnumValues)
        : points(3 * numEnumValues), pointsSet(false) {}
    
    BasePointsType(ManualPoints mpt)
        : points(3 * GetNumEnumValues(mpt)), pointsSet(false) {}

    double GetPointAt(PointsNamesType ptType, unsigned int index) {
        return points.at((static_cast<int>(ptType) - offset) * 3 + index);
    }

    void SetPointAt(PointsNamesType ptType, unsigned int index, double value) {
        points.at((static_cast<int>(ptType) - offset) * 3 + index) = value;
    }

    void SetPoint(PointsNamesType ptType, std::vector<double> value) {
        for (int ix = 0; ix < 3; ix++) {
            SetPointAt(ptType, ix, value.at(ix));
        }
    }

    virtual PointsNamesType FromKey(QString key) = 0;

    void SetPoint(QJsonObject json, QString key) {
        if (json[key].isUndefined()) {
            MITK_WARN << ("Undefined key: " + key + '\n').toStdString();
            return;
        }

        for (int ix = 0; ix < 3; ix++) {
            SetPointAt(FromKey(key), ix, json[key].toArray().at(ix).toDouble());
        }
    }

    void SetPointsFromJson(QJsonObject json) {
        QStringList keys = json.keys();
        foreach (QString key, keys) {
            SetPoint(json, key);

        }
    }

    inline void PointSet(bool value) { pointsSet = value; };
    inline void PointSetOn() { PointSet(true); };
    inline void PointSetOff() { PointSet(false); };
    inline bool IsPointSet() { return pointsSet; };

protected:
    std::vector<double> points;
    bool pointsSet;
    int offset;

private: 
    int GetNumEnumValues(ManualPoints mpt) {
        switch (mpt) {
            case ManualPoints::CYLINDERS:
                return static_cast<int>(PArt3) - static_cast<int>(SVC1) + 1; // Or any calculation that fits your needs
            case ManualPoints::SLICERS:
                return static_cast<int>(PArt_TIP) - static_cast<int>(SVC_SLICER1) + 1; // Calculate for SlicersPointsType
            case ManualPoints::VALVE_PLAINS:
                return static_cast<int>(PArt_WT_TIP) - static_cast<int>(Ao_WT_TIP) + 1; // Calculate for ValvePlainsPointsType
        }
    }
};

class CylinderPointsType : public BasePointsType {
public:
    CylinderPointsType()
        : BasePointsType(static_cast<int>(PArt3) - static_cast<int>(SVC1) + 1) {
            offset = static_cast<int>(SVC1);
        } 

    PointsNamesType FromKey(QString key) override {
        if (key == "SVC_1") return SVC1;
        if (key == "SVC_2") return SVC2;
        if (key == "SVC_3") return SVC3;
        if (key == "IVC_1") return IVC1;
        if (key == "IVC_2") return IVC2;
        if (key == "IVC_3") return IVC3;
        if (key == "Ao_1") return Ao1;
        if (key == "Ao_2") return Ao2;
        if (key == "Ao_3") return Ao3;
        if (key == "PArt_1") return PArt1;
        if (key == "PArt_2") return PArt2;
        if (key == "PArt_3") return PArt3;
        return SVC1; // Default value, you can adjust this as needed
    }
}; 

class SlicersPointsType : public BasePointsType {
public:
    SlicersPointsType()
        : BasePointsType(static_cast<int>(PArt_TIP) - static_cast<int>(SVC_SLICER1) + 1) {
            offset = static_cast<int>(SVC_SLICER1);
        } 

    PointsNamesType FromKey(QString key) override {
        if (key == "SVC_slicer_1") return SVC_SLICER1;
        if (key == "SVC_slicer_2") return SVC_SLICER2;
        if (key == "SVC_slicer_3") return SVC_SLICER3;
        if (key == "SVC_tip") return SVC_TIP;
        if (key == "IVC_slicer_1") return IVC_SLICER1;
        if (key == "IVC_slicer_2") return IVC_SLICER2;
        if (key == "IVC_slicer_3") return IVC_SLICER3;
        if (key == "IVC_tip") return IVC_TIP;
        if (key == "Ao_tip") return Ao_TIP;
        if (key == "PArt_tip") return PArt_TIP;
        return SVC_SLICER1; // Default value, you can adjust this as needed
    }
};

class ValvePlainsPointsType : public BasePointsType {
public:
    ValvePlainsPointsType() 
        : BasePointsType(static_cast<int>(PArt_WT_TIP) - static_cast<int>(Ao_WT_TIP) + 1) {
            offset = static_cast<int>(Ao_WT_TIP);
        } // Total number of enum values for ValvePlains

    PointsNamesType FromKey(QString key) override {
        if (key == "Ao_WT_TIP") return Ao_WT_TIP;
        if (key == "PArt_WT_TIP") return PArt_WT_TIP;
        return Ao_WT_TIP; // Default value, you can adjust this as needed
    }
};

struct ThicknessInfo { 
    QJsonObject json;
    double scale_factor;
    ThicknessInfo() {
        QStringList keys1, values1, keys2, values2;
        QStringList types;
        keys1 << "scale_factor" << "valve_WT_multiplier" << "valve_WT_svc_ivc_multiplier"
             << "ring_thickness_multiplier"
             << "LV_neck_WT_multiplier" << "RV_WT_multiplier" << "LA_WT_multiplier" 
             << "RA_WT_multiplier" << "Ao_WT_multiplier" << "PArt_WT_multiplier";
        keys2 << "valve_WT" << "valve_WT_svc_ivc"
             << "ring_thickness"
             << "LV_neck_WT" << "RV_WT" << "LA_WT" << "RA_WT" << "Ao_WT" << "PArt_WT";
        values1 << "2.50978817388" << "4" << "4"
               << "4"
               << "2.00" << "3.50" << "2.00" << "2.00" << "2.00" << "2.00";

        scale_factor = 2.50978817388;

        QStringList keys = keys1;
        QStringList values = values1;

        for (int ix = 0; ix < keys2.size(); ix++) {
            QString valuePrint = QString::number(scale_factor * values1.at(ix + 1).toDouble());
            values2 << valuePrint;
            std::cout << values2.at(ix).toStdString() << ": " << valuePrint.toStdString() << '\n';
        }
        keys.append(keys2);
        values.append(values2);

        for (int ix = 0; ix < keys.size(); ix++){
            types << "float";
        }
        
        json = CemrgCommonUtils::CreateJSONObject(keys, values, types);
    }
    // start here: do everything through json
    QStringList GetKeys() {
        return json.keys();
    }

    void SetKey(QString key, QString value) {
        json[key] = value;
        if (key.contains("multiplier")) {
            json[key.left(key.size() - 11)] = value.toDouble()*scale_factor;
        }
    }

    QJsonObject GetJson() {
        return json;
    }

    void Print() {
        QStringList keys = json.keys();
        for (int ix = 0; ix < keys.size(); ix++) {
            std::cout << keys.at(ix).toStdString() << ": " << json[keys.at(ix)].toDouble() << '\n';
        }
    }

    void SetKey(QString key, double value) {
        if (key == "scale_factor") {
            scale_factor = value;
        } else if (key.contains("multiplier")) {
            json[key.left(key.size() - 11)] = value*scale_factor;
        }
        json[key] = value;
    }

    void Update() {
        QStringList keys = json.keys();
        for (int ix = 0; ix < keys.size(); ix++) {
            if (keys.at(ix).contains("multiplier")) {
                json[keys.at(ix).left(keys.at(ix).size() - 11)] = json[keys.at(ix)].toDouble()*scale_factor;
            }
        }
    }

};

struct HeartLabels {
    QJsonObject json;
    HeartLabels() {
        QStringList keys, values, types;
        keys <<"LV_BP_label" <<"LV_myo_label" <<"RV_BP_label" <<"LA_BP_label" <<"RA_BP_label" <<"Ao_BP_label" <<"PArt_BP_label";
        keys <<"LPV1_label" << "LPV2_label" << "RPV1_label" << "RPV2_label" << "LAA_label" << "SVC_label" << "IVC_label";
        keys <<"LV_neck_label" << "RV_myo_label" << "LA_myo_label" << "RA_myo_label" << "Ao_wall_label" << "PArt_wall_label";
        keys <<"MV_label" << "TV_label" << "AV_label" << "PV_label";
        keys <<"plane_LPV1_label" <<"plane_LPV2_label" <<"plane_RPV1_label" <<"plane_RPV2_label" <<"plane_LAA_label" <<"plane_SVC_label" <<"plane_IVC_label";
        keys <<"LPV1_ring_label" << "LPV2_ring_label" << "RPV1_ring_label" << "RPV2_ring_label" << "LAA_ring_label" << "SVC_ring_label" << "IVC_ring_label";
        values << "1" << "2" << "3" << "4" << "5" << "6" << "7" << "8"
               << "9" << "10" << "11" << "12" << "13" << "14" << "101" << "103"
               << "104" << "105" << "106" << "107" << "201" << "202" << "203" << "204"
               << "205" << "206" << "207" << "208" << "209" << "210" << "211" << "221"
               << "222" << "223" << "224" << "225" << "226" << "227";
        for (int ix = 0; ix < keys.size(); ix++) {
            types << "int";
        }

        json = CemrgCommonUtils::CreateJSONObject(keys, values, types);
    }

    int GetLabel(QString key) {
        return json[key].toInt();
    }

    void SetLabel(QString key, int value) {
        json[key] = value;
    }

    QStringList GetKeys() {
        return json.keys();
    }

    void Save(QString dir, QString name) {
        QString path = dir + "/" + name;
        if (QFile::exists(path)) {
            // bool CemrgCommonUtils::ModifyJSONFile(QString dir, QString fname, QString key, QString value, QString type) {
            foreach (QString key, json.keys()) {
                CemrgCommonUtils::ModifyJSONFile(dir, name, key, QString::number(json[key].toDouble()), "int");
            }
        } else {
            CemrgCommonUtils::WriteJSONFile(json, dir, name);
        }
    }

    QJsonObject UniteJson(QJsonObject other) {
        QJsonObject newJson = json;
        for (const QString &key : other.keys()) {
            newJson[key] = other[key];
        }
        return newJson;
    }

    void Print() {
        QStringList keys = json.keys();
        for (int ix = 0; ix < keys.size(); ix++) {
            std::cout << keys.at(ix).toStdString() << ": " << json[keys.at(ix)].toInt() << '\n';
        }
    }
};

enum AtrialLandmarksType {NOT_SET = -1, LA_APEX=11, LA_SEPTUM=13, RA_APEX=15, RA_SEPTUM=17, RAA_APEX=19, ATRIAL_SEPTUM=21};
class PickedPointType { 
    public :
        PickedPointType() {
            seedIds = vtkSmartPointer<vtkIdList>::New();
            seedIds->Initialize();
            
            lineSeeds = vtkSmartPointer<vtkPolyData>::New();
            lineSeeds->Initialize();
            lineSeeds->SetPoints(vtkSmartPointer<vtkPoints>::New());

            seedLabels = std::vector<int>();
            pointNames = QStringList();
        }

        void AddPointFromSurface(mitk::Surface::Pointer surface, int pickedSeedId) { 
            double* point = surface->GetVtkPolyData()->GetPoint(pickedSeedId);
            seedIds->InsertNextId(pickedSeedId);
            lineSeeds->GetPoints()->InsertNextPoint(point);
            lineSeeds->Modified();
        }

        void AddPoint(double* point, int pickedSeedId) {
            seedIds->InsertNextId(pickedSeedId);
            lineSeeds->GetPoints()->InsertNextPoint(point);
            lineSeeds->Modified();
        }
        void PushBackLabel(AtrialLandmarksType label) {
            seedLabels.push_back(static_cast<int>(label));
        }

        void PushBackLabelFromAvailable(int index) {
            seedLabels.push_back(availableLabels.at(index));
            labelSet.at(index) = true;
        }

        bool AllLabelsSet() {
            bool result = true;
            if (IsEmpty()) {
                return false;
            }

            for (unsigned int ix = 0; ix < labelSet.size(); ix++) {
                if (!labelSet.at(ix)) {
                    result = false;
                    break;
                }
            }
            return result;
        }

        void CleanupLastPoint() {
            // Clean up last dropped seed point
            vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
            vtkSmartPointer<vtkPoints> points = lineSeeds->GetPoints();
            for (int i = 0; i < points->GetNumberOfPoints() - 1; i++) {
                newPoints->InsertNextPoint(points->GetPoint(i));
            }
            lineSeeds->SetPoints(newPoints);
            vtkSmartPointer<vtkIdList> newSeedIds = vtkSmartPointer<vtkIdList>::New();
            newSeedIds->Initialize();
            vtkSmartPointer<vtkIdList> roughSeedIds = seedIds;
            for (int i = 0; i < roughSeedIds->GetNumberOfIds() - 1; i++) {
                newSeedIds->InsertNextId(roughSeedIds->GetId(i));
            }
            seedIds = newSeedIds;
        }

        void Clear() {
            seedIds->Reset();
            lineSeeds->Reset();
            seedLabels.clear();
            pointNames.clear();
            availableLabels.clear();
        }

        bool IsEmpty() {
            return seedIds->GetNumberOfIds() == 0;
        }

        int NumIds() {
            return seedIds->GetNumberOfIds();
        }

        vtkSmartPointer<vtkPolyData> GetLineSeeds() {
            return lineSeeds;
        }

        void SetAvailableLabels(QStringList names, std::vector<int> labels) {
            pointNames = names;
            availableLabels = labels;
            for (unsigned int ix = 0; ix < availableLabels.size(); ix++) {
                labelSet.push_back(false);
            }
        }

        std::string ToString() {
            std::string res = "";
            for (unsigned int i = 0; i < seedLabels.size(); i++) {
                res += "Label: " + std::to_string(seedLabels.at(i)) + " Point ID: " + std::to_string(seedIds->GetId(i)) + '\n';
            }

            for (int j = 0; j < lineSeeds->GetPoints()->GetNumberOfPoints(); j++) {
                double* point = lineSeeds->GetPoints()->GetPoint(j);
                res += "Point: (" + std::to_string(point[0]) + ", " + std::to_string(point[1]) + ", " + std::to_string(point[2]) + ")\n";
            }

            for (unsigned int k = 0; k < seedLabels.size(); k++) {
                res += "Label: " + std::to_string(seedLabels.at(k));
            }

                return res;
        }

        std::string PrintVtxFile(QString name) {
            std::string res = ""; 
            
            // look for label with name in available labels
            int thisLabel;
            for (int i = 0; i < pointNames.size(); i++) {
                if (pointNames.at(i) == name) { 
                    thisLabel = availableLabels.at(i);
                    break;
                }
            }

            // look for ID with label in seedLabels
            int thisId;
            for (unsigned int j = 0; j < seedLabels.size(); j++) {
                if (seedLabels.at(j) == thisLabel) {
                    thisId = seedIds->GetId(j);
                    break;
                }
            }

            res = "1\nextra\n" + std::to_string(thisId) + '\n';
            return res;
        }
        
    private :
        std::vector<int> seedLabels;
        vtkSmartPointer<vtkIdList> seedIds;
        vtkSmartPointer<vtkPolyData> lineSeeds;

        QStringList pointNames;
        std::vector<int> availableLabels;
        std::vector<int> labelSet;

        QStringList saveFiles;
};

#endif
