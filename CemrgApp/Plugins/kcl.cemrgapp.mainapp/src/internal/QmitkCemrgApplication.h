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
 * CemrgApp Main App
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef QmitkCemrgApplication_h
#define QmitkCemrgApplication_h

#include <berryIApplication.h>

class QmitkCemrgApplication : public QObject, public berry::IApplication {

    Q_OBJECT
    Q_INTERFACES(berry::IApplication)

public:

    QmitkCemrgApplication();
    QmitkCemrgApplication(const QmitkCemrgApplication& other);

    QVariant Start(berry::IApplicationContext* context) override;
    void Stop() override;
};

#endif // QmitkCemrgApplication_h
