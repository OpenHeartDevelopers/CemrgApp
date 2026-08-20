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
 * NIfTI Header Tools
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

// MITK
#include <mitkLogMacros.h>

// Qt
#include <QString>
#include <QFile>
#include <QFileInfo>

// Standard library
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#include <zlib.h>

#include "CemrgNiftiUtils.h"


namespace {

// NIfTI-1 header layout. Field offsets in bytes, per the nifti1.h specification.
const int NIFTI_HDR_SIZE     = 348;
const int NIFTI_OFF_PIXDIM   = 76;
const int NIFTI_OFF_QFORM    = 252;
const int NIFTI_OFF_SFORM    = 254;
const int NIFTI_OFF_QUATERN  = 256;
const int NIFTI_OFF_QOFFSET  = 268;
const int NIFTI_OFF_SROW     = 280;
const int NIFTI_OFF_MAGIC    = 344;

struct NiftiHeaderFields {
    bool swap = false;          // header is opposite-endian to this machine
    short qformCode = 0;
    short sformCode = 0;
    float pixdim[8] = {0};
    float quatern[3] = {0};     // b, c, d
    float qoffset[3] = {0};
    float srow[3][4] = {{0}};   // srow_x, srow_y, srow_z
};

template <typename T> T SwapBytes(T value) {
    char* p = reinterpret_cast<char*>(&value);
    for (size_t i = 0; i < sizeof(T) / 2; i++)
        std::swap(p[i], p[sizeof(T) - 1 - i]);
    return value;
}

template <typename T> T ReadField(const char* buf, int offset, bool swap) {
    T value;
    memcpy(&value, buf + offset, sizeof(T));
    return swap ? SwapBytes(value) : value;
}

/**
 * @brief Read the leading NIFTI_HDR_SIZE bytes, transparently handling gzip.
 *
 * zlib's gzread reads uncompressed files verbatim, so this covers .nii and .nii.gz alike.
 * Compression is detected from the gzip magic bytes rather than the file extension.
 */
bool ReadNiftiHeaderBytes(QString path, char* buf, bool& isCompressed) {
    isCompressed = false;
    QFile probe(path);
    if (!probe.open(QIODevice::ReadOnly))
        return false;
    char magic[2] = {0, 0};
    isCompressed = (probe.read(magic, 2) == 2) &&
                   (static_cast<unsigned char>(magic[0]) == 0x1f) &&
                   (static_cast<unsigned char>(magic[1]) == 0x8b);
    probe.close();

    gzFile gz = gzopen(path.toStdString().c_str(), "rb");
    if (gz == nullptr)
        return false;
    int bytesRead = gzread(gz, buf, NIFTI_HDR_SIZE);
    gzclose(gz);
    return (bytesRead == NIFTI_HDR_SIZE);
}

bool ParseNiftiHeader(const char* buf, NiftiHeaderFields& out) {
    int sizeofHdr = ReadField<int>(buf, 0, false);
    if (sizeofHdr != NIFTI_HDR_SIZE) {
        if (SwapBytes(sizeofHdr) != NIFTI_HDR_SIZE)
            return false;
        out.swap = true;
    }

    // Only NIfTI-1 is handled; "n+1" is single-file, "ni1" is the .hdr/.img pair.
    if (memcmp(buf + NIFTI_OFF_MAGIC, "n+1\0", 4) != 0 &&
        memcmp(buf + NIFTI_OFF_MAGIC, "ni1\0", 4) != 0)
        return false;

    out.qformCode = ReadField<short>(buf, NIFTI_OFF_QFORM, out.swap);
    out.sformCode = ReadField<short>(buf, NIFTI_OFF_SFORM, out.swap);
    for (int i = 0; i < 8; i++)
        out.pixdim[i] = ReadField<float>(buf, NIFTI_OFF_PIXDIM + 4 * i, out.swap);
    for (int i = 0; i < 3; i++) {
        out.quatern[i] = ReadField<float>(buf, NIFTI_OFF_QUATERN + 4 * i, out.swap);
        out.qoffset[i] = ReadField<float>(buf, NIFTI_OFF_QOFFSET + 4 * i, out.swap);
    }
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 4; c++)
            out.srow[r][c] = ReadField<float>(buf, NIFTI_OFF_SROW + 16 * r + 4 * c, out.swap);
    return true;
}

/**
 * @brief Rebuild the 3x4 qform matrix from the quaternion fields.
 *
 * Port of nifti_quatern_to_mat44. The quaternion encodes a proper rotation only, so the
 * left-handed case is carried by qfac (pixdim[0]) negating the third column.
 */
bool QformToMatrix(const NiftiHeaderFields& h, double mat[3][4]) {
    double b = h.quatern[0], c = h.quatern[1], d = h.quatern[2];
    double a2 = 1.0 - (b * b + c * c + d * d);
    double a;
    if (a2 > 1e-7) {
        a = std::sqrt(a2);
    } else {
        double norm = std::sqrt(b * b + c * c + d * d);
        if (norm < 1e-12)
            return false;               // quaternion is entirely absent
        b /= norm; c /= norm; d /= norm;
        a = 0.0;
    }

    double dx = h.pixdim[1], dy = h.pixdim[2], dz = h.pixdim[3];
    if (dx <= 0.0 || dy <= 0.0 || dz <= 0.0)
        return false;                   // no usable voxel scaling
    double qfac = (h.pixdim[0] < 0.0) ? -1.0 : 1.0;
    double scale[3] = {dx, dy, dz * qfac};

    double rot[3][3] = {
        {a*a + b*b - c*c - d*d, 2.0*(b*c - a*d),       2.0*(b*d + a*c)},
        {2.0*(b*c + a*d),       a*a + c*c - b*b - d*d, 2.0*(c*d - a*b)},
        {2.0*(b*d - a*c),       2.0*(c*d + a*b),       a*a + d*d - c*c - b*b}
    };

    for (int r = 0; r < 3; r++) {
        for (int c2 = 0; c2 < 3; c2++)
            mat[r][c2] = rot[r][c2] * scale[c2];
        mat[r][3] = h.qoffset[r];
    }
    return true;
}

double SformDeterminant(const NiftiHeaderFields& h) {
    const float (*m)[4] = h.srow;
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
         - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
         + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

/**
 * @brief Largest absolute discrepancy between the qform-derived matrix and the sform.
 */
double MaxMatrixDifference(const double qform[3][4], const NiftiHeaderFields& h) {
    double worst = 0.0;
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 4; c++) {
            double q = qform[r][c], s = h.srow[r][c];
            worst = std::max(worst, std::fabs(q - s));
        }
    }
    return worst;
}

// Absolute, in the units of the matrix itself: millimetres for the translation column, and mm-per-voxel for the direction/spacing block.
const double QFORM_SFORM_TOLERANCE = 1e-4;

/**
 * @brief Overwrite qform_code in an uncompressed NIfTI, in place.
 */
bool WriteQformCodeUncompressed(QString path, short newCode, bool swap) {
    QFile file(path);
    if (!file.open(QIODevice::ReadWrite))
        return false;
    if (!file.seek(NIFTI_OFF_QFORM)) {
        file.close();
        return false;
    }
    short value = swap ? SwapBytes(newCode) : newCode;
    bool ok = (file.write(reinterpret_cast<const char*>(&value), sizeof(value)) == sizeof(value));
    file.close();
    return ok;
}

/**
 * @brief Overwrite qform_code in a gzipped NIfTI.
 *
 * A gzip stream cannot be patched in place, so the file is streamed through a temporary and
 * swapped in on success. Chunked to keep peak memory flat regardless of volume size.
 */
bool WriteQformCodeCompressed(QString path, short newCode, bool swap) {
    QString tempPath = path + ".qformfix.tmp";

    gzFile in = gzopen(path.toStdString().c_str(), "rb");
    if (in == nullptr)
        return false;
    gzFile out = gzopen(tempPath.toStdString().c_str(), "wb");
    if (out == nullptr) {
        gzclose(in);
        return false;
    }

    const int BUFFER_SIZE = 1 << 20;
    std::vector<char> buffer(BUFFER_SIZE);
    bool ok = true;
    bool headerPatched = false;

    while (true) {
        int bytesRead = gzread(in, buffer.data(), BUFFER_SIZE);
        if (bytesRead < 0) { ok = false; break; }
        if (bytesRead == 0) break;

        if (!headerPatched) {
            if (bytesRead < NIFTI_HDR_SIZE) { ok = false; break; }
            short value = swap ? SwapBytes(newCode) : newCode;
            memcpy(buffer.data() + NIFTI_OFF_QFORM, &value, sizeof(value));
            headerPatched = true;
        }

        if (gzwrite(out, buffer.data(), bytesRead) != bytesRead) { ok = false; break; }
    }

    gzclose(in);
    gzclose(out);

    if (!ok || !headerPatched) {
        QFile::remove(tempPath);
        return false;
    }

    // Swap the repaired file in. QFile::rename will not overwrite, so clear the target first.
    if (!QFile::remove(path)) {
        QFile::remove(tempPath);
        return false;
    }
    if (!QFile::rename(tempPath, path)) {
        MITK_ERROR << ("[NiftiQform] Repaired file left at: " + tempPath).toStdString();
        return false;
    }
    return true;
}

CemrgNiftiUtils::NiftiQformStatus AssessNiftiQform(
    QString path, QString& message, NiftiHeaderFields& fields, bool& isCompressed) {

    using Status = CemrgNiftiUtils::NiftiQformStatus;

    if (!QFileInfo::exists(path)) {
        message = "File does not exist: " + path;
        return Status::IoError;
    }

    std::vector<char> buffer(NIFTI_HDR_SIZE);
    if (!ReadNiftiHeaderBytes(path, buffer.data(), isCompressed)) {
        message = "Could not read a NIfTI header from: " + path;
        return Status::IoError;
    }

    if (!ParseNiftiHeader(buffer.data(), fields)) {
        message = "Not a NIfTI-1 file: " + path;
        return Status::NotNifti;
    }

    if (fields.qformCode != 0) {
        message = "qform_code is already set (" + QString::number(fields.qformCode) + ").";
        return Status::AlreadyValid;
    }

    if (fields.sformCode == 0) {
        message = "qform_code and sform_code are both 0 - the file declares no usable transform, "
                  "so there is nothing to validate a repair against.";
        return Status::CannotRepair;
    }

    double qform[3][4];
    if (!QformToMatrix(fields, qform)) {
        message = "qform_code is 0 and the quaternion/voxel-size fields are empty or degenerate, "
                  "so no transform can be recovered.";
        return Status::CannotRepair;
    }

    if (std::fabs(SformDeterminant(fields)) < 1e-12) {
        message = "The sform matrix is singular, so it cannot be used to validate the quaternion.";
        return Status::CannotRepair;
    }

    double difference = MaxMatrixDifference(qform, fields);
    if (difference > QFORM_SFORM_TOLERANCE) {
        message = "qform_code is 0 and the stored quaternion does not match the sform "
                  "(max difference " + QString::number(difference, 'e', 3) +
                  "). Enabling it would change the image geometry, so it needs manual review.";
        return Status::CannotRepair;
    }

    message = "qform_code is 0 but the stored quaternion reproduces the sform "
              "(max difference " + QString::number(difference, 'e', 3) +
              "), so enabling it is geometrically lossless.";
    return Status::Repairable;
}

} // namespace

QString CemrgNiftiUtils::NiftiQformStatusToString(NiftiQformStatus status) {
    switch (status) {
        case NiftiQformStatus::NotNifti:     return "NotNifti";
        case NiftiQformStatus::AlreadyValid: return "AlreadyValid";
        case NiftiQformStatus::Repairable:   return "Repairable";
        case NiftiQformStatus::CannotRepair: return "CannotRepair";
        case NiftiQformStatus::Repaired:     return "Repaired";
        case NiftiQformStatus::IoError:      return "IoError";
    }
    return "Unknown";
}

CemrgNiftiUtils::NiftiQformStatus CemrgNiftiUtils::InspectNiftiQform(QString pathToImage, QString& message) {
    NiftiHeaderFields fields;
    bool isCompressed = false;
    return AssessNiftiQform(pathToImage, message, fields, isCompressed);
}

CemrgNiftiUtils::NiftiQformStatus CemrgNiftiUtils::RepairNiftiQform(QString pathToImage, QString& message) {
    NiftiHeaderFields fields;
    bool isCompressed = false;
    NiftiQformStatus status = AssessNiftiQform(pathToImage, message, fields, isCompressed);

    if (status != NiftiQformStatus::Repairable)
        return status;

    // Mirror sform_code so the qform claims exactly the coordinate space the file already
    // declared - enabling the existing transform without asserting anything new about it.
    short newCode = fields.sformCode;
    bool written = isCompressed
        ? WriteQformCodeCompressed(pathToImage, newCode, fields.swap)
        : WriteQformCodeUncompressed(pathToImage, newCode, fields.swap);

    if (!written) {
        message = "Could not write the repaired header to " + pathToImage +
                  " (check file permissions and free disk space).";
        return NiftiQformStatus::IoError;
    }

    message = "Set qform_code to " + QString::number(newCode) +
              " (matching sform_code) in " + QFileInfo(pathToImage).fileName() +
              ". Image geometry is unchanged; the existing quaternion was already correct.";
    return NiftiQformStatus::Repaired;
}
