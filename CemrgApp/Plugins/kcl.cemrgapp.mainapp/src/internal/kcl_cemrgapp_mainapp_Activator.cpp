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

#include "kcl_cemrgapp_mainapp_Activator.h"
#include "perspectives/QmitkCemrgJBPerspective.h"
#include "perspectives/QmitkCemrgHCPerspective.h"
#include "perspectives/QmitkCemrgRRPerspective.h"
#include "perspectives/QmitkCemrgEasiPerspective.h"
#include "perspectives/QmitkCemrgPowertransPerspective.h"
#include "perspectives/QmitkCemrgWathcaPerspective.h"
#include "QmitkCemrgApplication.h"
#include "QmitkCemrgAppCommonTools.h"
#include <mitkVersion.h>
#include <berryLog.h>
#include <QFileInfo>
#include <QDateTime>

namespace mitk {

    kcl_cemrgapp_mainapp_Activator* kcl_cemrgapp_mainapp_Activator::inst = nullptr;

    kcl_cemrgapp_mainapp_Activator::kcl_cemrgapp_mainapp_Activator() {
        inst = this;
        context = nullptr;
    }

    kcl_cemrgapp_mainapp_Activator::~kcl_cemrgapp_mainapp_Activator() {
    }

    kcl_cemrgapp_mainapp_Activator* kcl_cemrgapp_mainapp_Activator::GetDefault() {
        return inst;
    }

    void kcl_cemrgapp_mainapp_Activator::start(ctkPluginContext* context) {

        berry::AbstractUICTKPlugin::start(context);
        this->context = context;

        BERRY_REGISTER_EXTENSION_CLASS(QmitkCemrgApplication, context);
        BERRY_REGISTER_EXTENSION_CLASS(QmitkCemrgAppCommonTools, context);
        BERRY_REGISTER_EXTENSION_CLASS(QmitkCemrgJBPerspective, context);
        BERRY_REGISTER_EXTENSION_CLASS(QmitkCemrgHCPerspective, context);
        BERRY_REGISTER_EXTENSION_CLASS(QmitkCemrgRRPerspective, context);
        BERRY_REGISTER_EXTENSION_CLASS(QmitkCemrgEasiPerspective, context);
        BERRY_REGISTER_EXTENSION_CLASS(QmitkCemrgPowertransPerspective, context);
        BERRY_REGISTER_EXTENSION_CLASS(QmitkCemrgWathcaPerspective, context);

        // QString collectionFile = GetQtHelpCollectionFile();
        // berry::QtAssistantUtil::SetHelpCollectionFile(collectionFile);
        // berry::QtAssistantUtil::SetDefaultHelpUrl("qthelp://kcl.cemrgapp.cemrgapp/bundle/index.html");
    }

    ctkPluginContext* kcl_cemrgapp_mainapp_Activator::GetPluginContext() const {
        return context;
    }

    QString kcl_cemrgapp_mainapp_Activator::GetQtHelpCollectionFile() const {

        if (!helpCollectionFile.isEmpty()) {
            return helpCollectionFile;
        }

        QString collectionFilename = "CemrgAppQtHelpCollection.qhc";
        QFileInfo collectionFileInfo = context->getDataFile(collectionFilename);
        QFileInfo pluginFileInfo = QFileInfo(QUrl(context->getPlugin()->getLocation()).toLocalFile());

        if (!collectionFileInfo.exists() || pluginFileInfo.lastModified()
            > collectionFileInfo.lastModified()) {

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
