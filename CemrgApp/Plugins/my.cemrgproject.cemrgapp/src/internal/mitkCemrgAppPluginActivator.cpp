/*=========================================================================

 Program:   BlueBerry Platform
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

#include "mitkCemrgAppPluginActivator.h"
#include "QmitkCemrgApplication.h"
#include "QmitkCemrgJBPerspective.h"
#include "QmitkCemrgHCPerspective.h"
#include "QmitkCemrgRRPerspective.h"
#include "QmitkCemrgEasiPerspective.h"
#include "QmitkCemrgFestivePerspective.h"
#include "QmitkCemrgWathcaPerspective.h"
#include <berryLog.h>
#include <mitkVersion.h>
#include <QFileInfo>
#include <QDateTime>


namespace mitk {

CemrgAppPluginActivator* CemrgAppPluginActivator::inst = 0;

CemrgAppPluginActivator::CemrgAppPluginActivator() {

    inst = this;
}

CemrgAppPluginActivator::~CemrgAppPluginActivator() {

}

CemrgAppPluginActivator* CemrgAppPluginActivator::GetDefault() {

    return inst;
}

void CemrgAppPluginActivator::start(ctkPluginContext* context) {

    berry::AbstractUICTKPlugin::start(context);
    this->context = context;

    BERRY_REGISTER_EXTENSION_CLASS(QmitkCemrgApplication, context);
    BERRY_REGISTER_EXTENSION_CLASS(QmitkCemrgJBPerspective, context);
    BERRY_REGISTER_EXTENSION_CLASS(QmitkCemrgHCPerspective, context);
    BERRY_REGISTER_EXTENSION_CLASS(QmitkCemrgRRPerspective, context);
    BERRY_REGISTER_EXTENSION_CLASS(QmitkCemrgEasiPerspective, context);
    BERRY_REGISTER_EXTENSION_CLASS(QmitkCemrgFestivePerspective, context);
    BERRY_REGISTER_EXTENSION_CLASS(QmitkCemrgWathcaPerspective, context);

    QString collectionFile = GetQtHelpCollectionFile();
    //berry::QtAssistantUtil::SetHelpCollectionFile(collectionFile);
    //berry::QtAssistantUtil::SetDefaultHelpUrl("qthelp://my.cemrgproject.cemrgapp/bundle/index.html");
}

ctkPluginContext* CemrgAppPluginActivator::GetPluginContext() const {

    return context;
}

QString CemrgAppPluginActivator::GetQtHelpCollectionFile() const {

    if (!helpCollectionFile.isEmpty()) {

        return helpCollectionFile;
    }

    QString collectionFilename = "CemrgAppQtHelpCollection.qhc";

    QFileInfo collectionFileInfo = context->getDataFile(collectionFilename);
    QFileInfo pluginFileInfo = QFileInfo(QUrl(context->getPlugin()->getLocation()).toLocalFile());
    if (!collectionFileInfo.exists() ||
            pluginFileInfo.lastModified() > collectionFileInfo.lastModified()) {

        // extract the qhc file from the plug-in
        QByteArray content = context->getPlugin()->getResource(collectionFilename);
        if (content.isEmpty()) {
            BERRY_WARN << "Could not get plug-in resource: " << collectionFilename.toStdString();
        } else {
            QFile file(collectionFileInfo.absoluteFilePath());
            file.open(QIODevice::WriteOnly);
            file.write(content);
            file.close();
        }
    }

    if (QFile::exists(collectionFileInfo.absoluteFilePath())) {
        helpCollectionFile = collectionFileInfo.absoluteFilePath();
    }

    return helpCollectionFile;
}

}

#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
#include <QtPlugin>
Q_EXPORT_PLUGIN2(my_cemrgproject_cemrgapp, mitk::CemrgAppPluginActivator)
#endif
