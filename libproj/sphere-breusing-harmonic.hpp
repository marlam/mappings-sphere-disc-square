#ifndef PROJ_BREUSING_HARMONIC_H
#define PROJ_BREUSING_HARMONIC_H

#include "proj.hpp"

namespace Projection {

class BreusingHarmonicMean : public Implementation
{
public:
    BreusingHarmonicMean();
    static void forward(const Context* ctx, const Implementation* imp, double lat, double lon, double& x, double& y);
    static void inverse(const Context* ctx, const Implementation* imp, double x, double y, double& lat, double& lon);
};

}

#endif
