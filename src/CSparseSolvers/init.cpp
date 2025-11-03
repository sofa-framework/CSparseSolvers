#include <CSparseSolvers/config.h>

#include <sofa/core/ObjectFactory.h>
using sofa::core::ObjectFactory;
#include <sofa/helper/system/PluginManager.h>

namespace csparsesolvers
{
    extern void registerSparseLUSolver(sofa::core::ObjectFactory* factory);
    extern void registerSparseCholeskySolver(sofa::core::ObjectFactory* factory);
}

extern "C" {
    SOFA_CSPARSESOLVERS_API void initExternalModule();
    SOFA_CSPARSESOLVERS_API const char* getModuleName();
    SOFA_CSPARSESOLVERS_API const char* getModuleVersion();
    SOFA_CSPARSESOLVERS_API const char* getModuleLicense();
    SOFA_CSPARSESOLVERS_API const char* getModuleDescription();
    SOFA_CSPARSESOLVERS_API void registerObjects(sofa::core::ObjectFactory* factory);
}

void initExternalModule()
{
    static bool first = true;
    if (first)
    {
        // make sure that this plugin is registered into the PluginManager
        sofa::helper::system::PluginManager::getInstance().registerPlugin(csparsesolvers::MODULE_NAME);

        first = false;
    }
}

const char* getModuleName()
{
    return csparsesolvers::MODULE_NAME;
}

const char* getModuleVersion()
{
    return csparsesolvers::MODULE_VERSION;
}

const char* getModuleLicense()
{
    return "LGPL";
}

const char* getModuleDescription()
{
    return "A set of linear solvers based on the library CSparse";
}

void registerObjects(sofa::core::ObjectFactory* factory)
{
    csparsesolvers::registerSparseLUSolver(factory);
    csparsesolvers::registerSparseCholeskySolver(factory);
}
