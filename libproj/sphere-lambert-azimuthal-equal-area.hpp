#ifndef PROJ_LAMBERT_AZ_AE_H
#define PROJ_LAMBERT_AZ_AE_H

#include "proj.hpp"

namespace Projection {

class LambertAzimuthalEqualArea : public Implementation
{
public:
    LambertAzimuthalEqualArea();
    static void forward(const Context* ctx, const Implementation* imp, double lat, double lon, double& x, double& y);
    static void inverse(const Context* ctx, const Implementation* imp, double x, double y, double& lat, double& lon);
};

}

#endif
