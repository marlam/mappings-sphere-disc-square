#ifndef PROJ_STEREOGRAPHIC_H
#define PROJ_STEREOGRAPHIC_H

#include "proj.hpp"

namespace Projection {

class Stereographic : public Implementation
{
public:
    Stereographic();
    static void forward(const Context* ctx, const Implementation* imp, double lat, double lon, double& x, double& y);
    static void inverse(const Context* ctx, const Implementation* imp, double x, double y, double& lat, double& lon);
};

}

#endif
