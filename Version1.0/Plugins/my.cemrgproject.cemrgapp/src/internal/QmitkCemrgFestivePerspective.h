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
 * Fetal Motion Corrected Slice to Volume Registration (Festive) Plugin for MITK
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrg.co.uk/
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef QMITKCEMRGFESTIVEPERSPECTIVE_H_
#define QMITKCEMRGFESTIVEPERSPECTIVE_H_

#include <berryIPerspectiveFactory.h>

class QmitkCemrgFestivePerspective : public QObject, public berry::IPerspectiveFactory {

    Q_OBJECT
    Q_INTERFACES(berry::IPerspectiveFactory)

public:

    QmitkCemrgFestivePerspective();
    QmitkCemrgFestivePerspective(const QmitkCemrgFestivePerspective& other);

    void CreateInitialLayout(berry::IPageLayout::Pointer layout);

};

#endif /* QMITKCEMRGFESTIVEPERSPECTIVE */
