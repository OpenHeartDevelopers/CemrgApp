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
 * Eikonal Activation Simulation (EASI) Plugin for MITK
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef QMITKCEMRGEASIPERSPECTIVE_H_
#define QMITKCEMRGEASIPERSPECTIVE_H_

#include <berryIPerspectiveFactory.h>

class QmitkCemrgEasiPerspective : public QObject, public berry::IPerspectiveFactory {

    Q_OBJECT
    Q_INTERFACES(berry::IPerspectiveFactory)

public:

    QmitkCemrgEasiPerspective();
    QmitkCemrgEasiPerspective(const QmitkCemrgEasiPerspective& other);

    void CreateInitialLayout(berry::IPageLayout::Pointer layout);

};

#endif /* QMITKCEMRGEASIPERSPECTIVE */
