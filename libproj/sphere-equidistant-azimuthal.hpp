#ifndef PROJ_EQUIDISTANT_AZ_H
#define PROJ_EQUIDISTANT_AZ_H

#include "proj.hpp"

namespace Projection {

class EquidistantAzimuthal : public Implementation
{
public:
    EquidistantAzimuthal();
    static void forward(const Context* ctx, const Implementation* imp, double lat, double lon, double& x, double& y);
    static void inverse(const Context* ctx, const Implementation* imp, double x, double y, double& lat, double& lon);
};

}

#endif
