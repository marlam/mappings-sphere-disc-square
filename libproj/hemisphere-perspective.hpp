#ifndef PROJ_PERSPECTIVE_H
#define PROJ_PERSPECTIVE_H

#include "proj.hpp"

namespace Projection {

class Perspective : public Implementation
{
private:
    double _P, _lat0, _lon0;
public:
    Perspective(double P, double lat0, double lon0);
    static void forward(const Context* ctx, const Implementation* imp, double lat, double lon, double& x, double& y);
    static void inverse(const Context* ctx, const Implementation* imp, double x, double y, double& lat, double& lon);
};

}

#endif
