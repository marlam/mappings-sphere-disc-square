#include "hemisphere-stereographic.hpp"

namespace Projection {

Stereographic::Stereographic()
{
}

void Stereographic::forward(const Context* ctx, const Implementation* /* imp */, double lat, double lon, double& x, double& y)
{
    if (!ctx->is_sphere)
        lat = ctx->geodetic_to_geocentric(lat);

    if (ctx->center == CenterSouthPole)
        lat = -lat;

    double r = tan((M_PI_2 - lat) / 2);
    double phi = lon - M_PI_2;

    x = r * cos(phi);
    y = r * sin(phi);
}

void Stereographic::inverse(const Context* ctx, const Implementation* /* imp */, double x, double y, double& lat, double& lon)
{
    double r = hypot(x, y);
    double phi = atan2(y, x);

    lat = M_PI_2 - 2 * atan(r);
    lon = phi + M_PI_2;

    if (ctx->center == CenterSouthPole)
        lat = -lat;

    if (!ctx->is_sphere)
        lat = ctx->geocentric_to_geodetic(lat);
}

}
