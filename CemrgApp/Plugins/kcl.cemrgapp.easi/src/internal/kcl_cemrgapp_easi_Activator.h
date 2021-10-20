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
 * http://www.cemrg.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#ifndef kcl_cemrgapp_easi_Activator_h
#define kcl_cemrgapp_easi_Activator_h

#include <ctkPluginActivator.h>

namespace mitk {

    class kcl_cemrgapp_easi_Activator: public QObject, public ctkPluginActivator {

        Q_OBJECT
        Q_PLUGIN_METADATA(IID "kcl_cemrgapp_easi")
        Q_INTERFACES(ctkPluginActivator)

    public:

        void start(ctkPluginContext *context);
        void stop(ctkPluginContext *context);
        static ctkPluginContext* getContext();

    private:

        static ctkPluginContext* pluginContext;

    }; // kcl_cemrgapp_easi_Activator
}

#endif // kcl_cemrgapp_easi_Activator_h
