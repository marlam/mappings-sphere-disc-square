#ifndef PROJ_EQAR_CONF_TRADEOFF_H
#define PROJ_EQAR_CONF_TRADEOFF_H

#include "proj.hpp"

namespace Projection {

class MixLambertStereographic : public Implementation
{
public:
    double beta;

    MixLambertStereographic(double beta);
    static void forward(const Context* ctx, const Implementation* imp, double lat, double lon, double& x, double& y);
    static void inverse(const Context* ctx, const Implementation* imp, double x, double y, double& lat, double& lon);
};

}

#endif
