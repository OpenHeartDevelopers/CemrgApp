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

#ifndef CemrgNiftiUtils_h
#define CemrgNiftiUtils_h

#include <MitkCemrgAppModuleExports.h>
#include <QString>

class MITKCEMRGAPPMODULE_EXPORT CemrgNiftiUtils {

public:

    //Nifti Header Repair Utils
    /**
     * @brief Outcome of inspecting a NIfTI file's qform/sform declaration.
     *
     * ITK (5.2.x) determines image direction from the qform alone. A file with qform_code=0 is
     * rejected at load with "ITK only supports orthonormal direction cosines. No orthonormal
     * definition found!" even when its sform is perfectly well formed. Such a file is repairable
     * when it already carries a quaternion that reproduces its sform.
     */
    enum class NiftiQformStatus {
        NotNifti,       ///< Not a readable NIfTI-1 file.
        AlreadyValid,   ///< qform_code is already non-zero; nothing to do.
        Repairable,     ///< qform_code=0 but the stored quaternion matches the sform.
        CannotRepair,   ///< qform_code=0 and no trustworthy quaternion to enable.
        Repaired,       ///< qform_code was enabled (RepairNiftiQform only).
        IoError         ///< The file could not be read or written.
    };

    static NiftiQformStatus InspectNiftiQform(QString pathToImage, QString& message);
    static NiftiQformStatus RepairNiftiQform(QString pathToImage, QString& message);
    static QString NiftiQformStatusToString(NiftiQformStatus status);
};

#endif // CemrgNiftiUtils_h
