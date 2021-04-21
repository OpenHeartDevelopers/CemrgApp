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

#ifndef QMITKCEMRGAPPLICATION_H_
#define QMITKCEMRGAPPLICATION_H_

#include <berryIApplication.h>


class QmitkCemrgApplication : public QObject, public berry::IApplication {

    Q_OBJECT
    Q_INTERFACES(berry::IApplication)

public:

    QmitkCemrgApplication();
    QmitkCemrgApplication(const QmitkCemrgApplication& other);

    QVariant Start(berry::IApplicationContext* context);
    void Stop();
};

#endif /*QMITKCEMRGAPPLICATION_H_*/
