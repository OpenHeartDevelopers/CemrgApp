#include "mitkPluginActivator.h"
#include "WallThicknessCalculationsView.h"
#include "WallThicknessCalculationsClipperView.h"

namespace mitk {

ctkPluginContext* PluginActivator::pluginContext = nullptr;

void PluginActivator::start(ctkPluginContext* context) {

    BERRY_REGISTER_EXTENSION_CLASS(WallThicknessCalculationsView, context)
    BERRY_REGISTER_EXTENSION_CLASS(WallThicknessCalculationsClipperView, context)
    pluginContext = context;
}

void PluginActivator::stop(ctkPluginContext* context) {

    Q_UNUSED(context)
    pluginContext = nullptr;
}

ctkPluginContext* PluginActivator::getContext() {
    return pluginContext;
}

}

#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
#include <QtPlugin>
Q_EXPORT_PLUGIN2(my.cemrgproject.wathca, mitk::PluginActivator)
#endif
