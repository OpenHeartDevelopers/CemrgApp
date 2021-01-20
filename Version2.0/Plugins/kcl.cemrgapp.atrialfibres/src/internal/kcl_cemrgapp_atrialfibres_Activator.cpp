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
 * Morphological Quantification
 *
 * Cardiac Electromechanics Research Group
 * http://www.cemrgapp.com
 * orod.razeghi@kcl.ac.uk
 *
 * This software is distributed WITHOUT ANY WARRANTY or SUPPORT!
 *
=========================================================================*/

#include "kcl_cemrgapp_atrialfibres_Activator.h"
#include "AtrialFibresView.h"
#include "AtrialFibresClipperView.h"

namespace mitk{

ctkPluginContext* kcl_cemrgapp_atrialfibres_Activator::pluginContext = nullptr;

void kcl_cemrgapp_atrialfibres_Activator::start(ctkPluginContext *context) {

    BERRY_REGISTER_EXTENSION_CLASS(AtrialFibresView, context);
    BERRY_REGISTER_EXTENSION_CLASS(AtrialFibresClipperView, context);
    pluginContext = context;
}

void kcl_cemrgapp_atrialfibres_Activator::stop(ctkPluginContext *context) {

    Q_UNUSED(context)
    pluginContext = nullptr;
}

ctkPluginContext* kcl_cemrgapp_atrialfibres_Activator::getContext() {

    return pluginContext;
}
}
