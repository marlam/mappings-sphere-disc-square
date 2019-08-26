#include "hemisphere-perspective.hpp"

namespace Projection {

Perspective::Perspective(double P, double lat0, double lon0) :
    _P(P), _lat0(lat0), _lon0(lon0)
{
}

void Perspective::forward(const Context* ctx, const Implementation* imp, double lat, double lon, double& x, double& y)
{
    const Perspective* me = reinterpret_cast<const Perspective*>(imp);

    if (!ctx->is_sphere)
        lat = ctx->geodetic_to_geocentric(lat);

    double P = me->_P;
    double lat0 = me->_lat0;
    double lon0 = me->_lon0;

    double cos_c = sin(lat0) * sin(lat) + cos(lat0) * cos(lat) * cos(lon - lon0);
    if (cos_c < 1.0 / P) {
        x = y = NAN;
        return;
    }
    double kp = (P - 1.0) / (P - cos_c);
    double scale = 1.0 / sqrt((P-1)/(P+1));
    x = scale * kp * cos(lat) * sin(lon - lon0);
    y = scale * kp * (cos(lat0) * sin(lat) - sin(lat0) * cos(lat) * cos(lon - lon0));
}

void Perspective::inverse(const Context* /* ctx */, const Implementation* /* imp */, double /* x */, double /* y */, double& lat, double& lon)
{
    lat = lon = NAN;
}

}
