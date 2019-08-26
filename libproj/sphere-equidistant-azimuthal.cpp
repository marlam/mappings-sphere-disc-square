#include "sphere-equidistant-azimuthal.hpp"

namespace Projection {

EquidistantAzimuthal::EquidistantAzimuthal()
{
}

void EquidistantAzimuthal::forward(const Context* ctx, const Implementation* /* imp */, double lat, double lon, double& x, double& y)
{
    if (!ctx->is_sphere)
        lat = ctx->geodetic_to_geocentric(lat);

    if (ctx->center == CenterSouthPole)
        lat = -lat;

    double r = 1 - lat / M_PI_2;
    double phi = lon - M_PI_2;

    if (ctx->layout == LayoutSphere)
        r /= 2;

    x = r * cos(phi);
    y = r * sin(phi);
}

void EquidistantAzimuthal::inverse(const Context* ctx, const Implementation* /* imp */, double x, double y, double& lat, double& lon)
{
    double r = hypot(x, y);
    double phi = atan2(y, x);

    if (ctx->layout == LayoutSphere)
        r *= 2;

    lat = (1 - r) * M_PI_2;
    lon = phi + M_PI_2;

    if (ctx->center == CenterSouthPole)
        lat = -lat;

    if (!ctx->is_sphere)
        lat = ctx->geocentric_to_geodetic(lat);
}

}
